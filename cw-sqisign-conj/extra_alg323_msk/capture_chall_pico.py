#!/usr/bin/env python3
"""
Capture traces of SQIsign Algorithm 3.23 line 12: (c1, c2) <- Msk * [1, achl]^T.

Unlike the Algorithm 3.18 line 2 target, the secret here is fixed: Msk is part of
the secret key and is identical in every signature, while achl = chall_coeff is a
hash output that ships inside the signature. Fixed secret plus known varying
input is what CPA needs, and traces average.

A trace costs one short serial exchange rather than a full Algorithm 3.18 run, so
this captures roughly three orders of magnitude faster than capture_conj_pico.py.

Modes:
  random   fresh uniform challenge coefficient per trace -- the CPA set
  fixed    one coefficient reused -- noise floor, and the fixed group of a TVLA
  tvla     fixed and random interleaved

The secret matrix is generated from --secret-seed and stored in the npz, since
this is a profiling experiment where the answer has to be known to score the
attack. A real attacker has the traces and the coefficients but not the matrix.

Usage:
  ~/.cwvenv/bin/python3 capture_chall_pico.py --mode tvla   --traces 2000
  ~/.cwvenv/bin/python3 capture_chall_pico.py --mode random --traces 20000 --skip-flash
"""

import argparse
import pathlib
import random
import sys
import time

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import sys as _sys, pathlib as _pathlib
# conj_rig and the analysis modules live in ../common, shared by both parts
_sys.path.insert(0, str(_pathlib.Path(__file__).resolve().parent.parent / "common"))

import conj_rig as rig

CHALLENGE_BYTES = 16
CHALLENGE_BITS = 128


def make_secret(seed):
    r = random.Random(seed)
    # keygen reduces every entry mod 2^CHALLENGE_BITS and normalises it positive
    return [r.getrandbits(CHALLENGE_BITS) for _ in range(4)]


def load_secret(target, msk):
    txt = "".join(f"{m:x}\n" for m in msk).encode()
    rig.load_vector(target, txt)
    if rig.ss_cmd(target, "k", b"", 1, timeout_ms=5000) != bytearray([1]):
        raise RuntimeError("loading the secret matrix failed")


def verify(target, msk, a):
    """Check the board's product against the host before trusting a campaign."""
    rig.ss_cmd(target, "f", b"", 1, timeout_ms=5000)
    lines = [l for l in rig.read_outbuf(target).splitlines() if l.strip()]
    c1, c2 = int(lines[0], 16), int(lines[1], 16)
    prof = {k: int(v) for k, v in (l.split("=") for l in lines[2:] if "=" in l)}
    ok = (c1 == msk[0] + msk[1] * a) and (c2 == msk[2] + msk[3] * a)
    return ok, prof


def main():
    p = argparse.ArgumentParser(description="Capture Alg 3.23 line 12 traces")
    p.add_argument("--mode", choices=["random", "fixed", "tvla"], default="random")
    p.add_argument("-n", "--traces", type=int, default=5000)
    p.add_argument("--samples", type=int, default=60000)
    p.add_argument("--timebase", type=int, default=2)
    p.add_argument("--pre-samples", type=int, default=1000)
    p.add_argument("--cha-range", type=int, default=3)
    p.add_argument("--chb-range", type=int, default=8)
    p.add_argument("--skip-flash", action="store_true")
    p.add_argument("--secret-seed", type=int, default=1)
    p.add_argument("--coeff-seed", type=int, default=12345)
    p.add_argument("--crop", action="store_true", default=True)
    p.add_argument("--no-crop", dest="crop", action="store_false")
    p.add_argument("--tail", type=int, default=9000,
                   help="samples to keep past the last trigger pulse, so the "
                        "trailing add stays in the window")
    p.add_argument("--out", default=None)
    p.add_argument("--plot", default=None)
    args = p.parse_args()

    out = pathlib.Path(args.out or f"/tmp/sqisign_chall_{args.mode}.npz")
    plot = pathlib.Path(args.plot or f"/tmp/sqisign_chall_{args.mode}.png")

    if not rig.HEX.exists():
        print(f"ERROR: firmware not found: {rig.HEX}\n"
              f"Build it with: make PLATFORM=CW308_STM32F4 CONJ_TRIG_MODE=12", file=sys.stderr)
        sys.exit(1)

    import chipwhisperer as cw
    from chipwhisperer.capture.api.programmers import STM32FProgrammer

    msk = make_secret(args.secret_seed)
    rnd = random.Random(args.coeff_seed)
    fixed_coeff = rnd.getrandbits(CHALLENGE_BITS)

    scope = cw.scope()
    target = cw.target(scope)
    pico = None
    traces, coeffs, groups = [], [], []
    pulse_bounds = []

    try:
        scope.default_setup()
        rig.set_target_clock(scope)
        if not args.skip_flash:
            print(f"programming {rig.HEX.name} ...", flush=True)
            cw.program_target(scope, STM32FProgrammer, str(rig.HEX))
        else:
            rig.reset_target(scope)
        time.sleep(0.2)
        target.flush()
        print("banner:", target.read(timeout=1000), flush=True)
        rig.sync(target)

        info = rig.board_info(target)
        print(f"board: trig_mode={info['trig_mode']}", flush=True)
        if info["trig_mode"] not in (10, 11, 12, 13, 14):
            print(f"WARNING: CONJ_TRIG_MODE={info['trig_mode']} is not one of the "
                  f"Alg 3.23 line 12 modes (10-14)", flush=True)

        load_secret(target, msk)
        print(f"secret matrix loaded (seed {args.secret_seed})", flush=True)

        # correctness gate: a campaign on a board that computes the wrong thing is
        # worth nothing, so check once against the host before capturing
        a0 = rnd.getrandbits(CHALLENGE_BITS)
        rig.ss_cmd(target, "e", a0.to_bytes(CHALLENGE_BYTES, "little"), 1, timeout_ms=5000)
        ok, prof = verify(target, msk, a0)
        if not ok:
            raise RuntimeError("board product does not match the host")
        print(f"board verified against host. targeted multiply = {prof.get('mul1')} cycles, "
              f"secret_limbs={prof.get('secret_limbs')} known_limbs={prof.get('known_limbs')}",
              flush=True)

        pico = rig.Pico(args.samples, args.timebase, args.pre_samples,
                        args.cha_range, args.chb_range, auto_trigger_ms=0)

        def one_capture(a, tries=5):
            """Arm, run one eval, return the record. Retries a missed trigger.

            A capture can be missed while the scope is still settling after
            RunBlock; the target is only ~20 ms behind the arm here, so that
            window is reachable. The eval is idempotent, so retrying is free.
            """
            payload = bytearray(a.to_bytes(CHALLENGE_BYTES, "little"))
            for attempt in range(tries):
                pico.arm()
                target.flush()
                target.simpleserial_write("e", payload)
                try:
                    # a miss costs this much, so keep it short: the target fires
                    # ~20 ms after arming, so anything beyond that is a miss
                    rec = pico.get(timeout_s=0.5)
                except TimeoutError:
                    target.simpleserial_read("r", 1, timeout=5000)
                    continue
                if target.simpleserial_read("r", 1, timeout=5000) != bytearray([1]):
                    raise RuntimeError("board reported eval failure")
                return rec
            raise TimeoutError(f"no trigger in {tries} attempts")

        for _ in range(3):  # AC coupling settles, and the scope warms up
            one_capture(fixed_coeff)

        t0 = time.time()
        for i in range(args.traces):
            if args.mode == "fixed":
                is_fixed = True
            elif args.mode == "random":
                is_fixed = False
            else:
                is_fixed = (i % 2 == 0)
            a = fixed_coeff if is_fixed else rnd.getrandbits(CHALLENGE_BITS)

            try:
                power, trig, overflow = one_capture(a)
            except TimeoutError as e:
                # The scope stops honouring triggers after some thousands of
                # arms. Keep what was captured rather than losing the run.
                print(f"\n{e} at trace {i}; keeping the {len(traces)} traces so far",
                      flush=True)
                break

            if args.crop:
                pulses = rig.trigger_pulses(trig)
                if not pulses:
                    raise RuntimeError(f"trace {i}: no trigger pulse")
                # Keep a tail past the last pulse. The final operation of
                # ibz_mat_2x2_eval is an add that happens after the last multiply,
                # so cropping at the last falling edge drops it entirely -- which
                # is what left Msk[1][0] with no target intermediate in the data.
                power = power[pulses[0][0]:min(len(power), pulses[-1][1] + args.tail)]

            traces.append(power.astype(np.float32))
            coeffs.append(a)
            groups.append(1 if is_fixed else 0)
            # Trigger pulse boundaries, relative to the crop. With mode 11 these
            # mark each multiply exactly, which beats deriving the operation
            # windows from the cycle profile -- that estimate overlapped the
            # m10 and m11 windows and let one attack pick up the other's leak.
            if pulse_bounds is not None:
                base = pulses[0][0] if args.crop else 0
                pulse_bounds.append([(s0 - base, e0 - base) for s0, e0, _ in pulses])

            if i == 0 or (i + 1) % 100 == 0 or i == args.traces - 1:
                rate = (time.time() - t0) / (i + 1)
                print(f"{i + 1:6d}/{args.traces}  len={len(power)}  "
                      f"{'OVF' if overflow & 1 else '   '}  "
                      f"{rate * 1000:.0f} ms/trace  eta={rate * (args.traces - i - 1) / 60:.1f} min",
                      flush=True)

    finally:
        if pico is not None:
            pico.close()
        try:
            target.dis()
        finally:
            scope.dis()

    if not traces:
        print("no traces captured", file=sys.stderr)
        return 1

    n = min(len(t) for t in traces)
    T = np.vstack([t[:n] for t in traces])

    # challenge coefficients as limbs, which is what the attack conditions on
    C = np.array(coeffs, dtype=object)
    limbs = np.zeros((len(coeffs), CHALLENGE_BITS // 32), dtype=np.uint32)
    for i, a in enumerate(coeffs):
        for k in range(CHALLENGE_BITS // 32):
            limbs[i, k] = (a >> (32 * k)) & 0xFFFFFFFF

    np.savez_compressed(
        out,
        traces=T,
        group=np.array(groups, dtype=np.uint8),
        coeff_limbs=limbs,
        coeff_hex=np.array([f"{a:032x}" for a in coeffs]),
        pulse_bounds=np.array(pulse_bounds[0] if pulse_bounds else [], dtype=np.int64),
        n_pulses=np.int64(len(pulse_bounds[0]) if pulse_bounds else 0),
        secret_hex=np.array([f"{m:032x}" for m in msk]),
        secret_seed=np.int64(args.secret_seed),
        interval_ns=np.float32(pico.interval_ns if pico else 0),
        samples_per_cycle=np.float32(pico.samples_per_cycle if pico else 0),
        trig_mode=np.uint8(info["trig_mode"]),
    )
    print(f"\nsaved {T.shape} to {out}  (fixed={sum(groups)}, random={len(groups) - sum(groups)})",
          flush=True)

    fig, ax = plt.subplots(2, 1, figsize=(14, 7))
    for i in range(min(8, T.shape[0])):
        ax[0].plot(T[i], lw=0.3, alpha=0.6)
    ax[0].set_title(f"Alg 3.23 line 12, first {min(8, T.shape[0])} of {T.shape[0]} traces")
    ax[0].set_ylabel("mV")
    ax[0].grid(alpha=0.3)
    ax[1].plot(T.mean(axis=0), lw=0.6)
    ax[1].set_title(f"mean of {T.shape[0]} traces")
    ax[1].set_xlabel("sample")
    ax[1].set_ylabel("mV")
    ax[1].grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(plot, dpi=130)
    plt.close(fig)
    print(f"saved {plot}", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
