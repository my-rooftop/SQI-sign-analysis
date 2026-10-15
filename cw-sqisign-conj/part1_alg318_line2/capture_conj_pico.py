#!/usr/bin/env python3
"""
Step 2: capture traces of SQIsign Algorithm 3.18 line 2, Conjugate(alpha_sk,chl).

Three capture sets, all written as .npz:

  --mode random    a fresh Alg 3.18 input and board seed per trace. This is the
                   SNR / profiling set: alpha_sk,chl varies, and every trace ships
                   the raw limbs the trace should be leaking.

  --mode fixed     one input and one seed, reused for every trace. Establishes the
                   noise floor, and is the "fixed" group of a TVLA.

  --mode tvla      interleaves fixed and random groups trace by trace, which is
                   what a non-specific TVLA wants: interleaving keeps slow drift
                   (temperature, supply) from showing up as leakage.

Each trace costs one full Algorithm 3.18 run. Line 1 (EquivalentCoprimeIdeal)
dominates that cost and runs before the trigger fires, so the scope waits on the
hardware trigger with no auto-trigger timeout. Budget accordingly: use
check_trigger.py first to measure the per-run wall time.

Usage:
  ~/.cwvenv/bin/python3 capture_conj_pico.py --mode random --traces 200
  ~/.cwvenv/bin/python3 capture_conj_pico.py --mode tvla --traces 200 --skip-flash
"""

import argparse
import json
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

VECTOR = pathlib.Path("/tmp/sqisign_conj_vec.txt")
FIXED_SEED16 = bytes(range(16))


def build_labels(records):
    """Per-coordinate, per-limb labels for SNR/CPA, taken from the pre-Conjugate value.

    ibz_neg walks limb 0 first, so limb index here is also the loop iteration index
    and therefore the time order inside the triggered window.
    """
    labels = {}
    nlimbs = max(r["pre_limbs"][c]["nlimbs"] for r in records for c in range(4))
    for c in range(4):
        vals = [r["pre"]["coord"][c] for r in records]
        labels[f"c{c}_sign"] = np.array([1 if v < 0 else 0 for v in vals], dtype=np.uint8)
        labels[f"c{c}_magnitude_bits"] = np.array([abs(v).bit_length() for v in vals], dtype=np.uint16)
        labels[f"c{c}_bitlen"] = np.array([r["pre_limbs"][c]["bitlen"] for r in records], dtype=np.uint16)
        labels[f"c{c}_nlimbs"] = np.array([r["pre_limbs"][c]["nlimbs"] for r in records], dtype=np.uint8)
        labels[f"c{c}_signext_limbs"] = np.array(
            [rig.sign_ext_limbs(r["pre_limbs"][c]["limbs"]) for r in records], dtype=np.uint8)
        for k in range(nlimbs):
            lv = np.array([
                (r["pre_limbs"][c]["limbs"][k] if k < r["pre_limbs"][c]["nlimbs"] else 0)
                for r in records], dtype=np.uint32)
            labels[f"c{c}_limb{k}"] = lv
            # value loaded by ibz_neg
            labels[f"c{c}_limb{k}_hw"] = np.array([rig.hw(v) for v in lv], dtype=np.uint8)
            # value it stores back: ~limb + carry, carry ignored here (it is 1 only
            # until the first non-zero limb); the complement is what hits the bus
            labels[f"c{c}_limb{k}_neg_hw"] = np.array(
                [rig.hw(~int(v)) for v in lv], dtype=np.uint8)
    return labels, nlimbs


def main():
    p = argparse.ArgumentParser(description="Capture Alg 3.18 line 2 Conjugate traces")
    p.add_argument("--mode", choices=["random", "fixed", "tvla"], default="random")
    p.add_argument("-n", "--traces", type=int, default=100,
                   help="total traces; in tvla mode this is split evenly between groups")
    p.add_argument("--samples", type=int, default=200000)
    p.add_argument("--timebase", type=int, default=2)
    p.add_argument("--pre-samples", type=int, default=2000)
    p.add_argument("--cha-range", type=int, default=3)
    p.add_argument("--chb-range", type=int, default=8)
    p.add_argument("--skip-flash", action="store_true")
    p.add_argument("--timeout", type=float, default=900.0)
    p.add_argument("--crop", action="store_true",
                   help="store only the triggered window instead of the whole record")
    p.add_argument("--seed-prefix", default="rnd",
                   help="prefix for the per-trace vector seed; change it to capture a "
                        "set disjoint from an earlier run instead of regenerating the "
                        "same vectors")
    p.add_argument("--out", default=None)
    p.add_argument("--plot", default=None)
    args = p.parse_args()

    out = pathlib.Path(args.out or f"/tmp/sqisign_conj_{args.mode}.npz")
    plot = pathlib.Path(args.plot or f"/tmp/sqisign_conj_{args.mode}.png")

    if not rig.HEX.exists():
        print(f"ERROR: firmware not found: {rig.HEX}\n"
              f"Build it with: make PLATFORM=CW308_STM32F4 CONJ_TRIG_MODE=2", file=sys.stderr)
        sys.exit(1)
    rig.ensure_host_tool()

    import chipwhisperer as cw
    from chipwhisperer.capture.api.programmers import STM32FProgrammer

    scope = cw.scope()
    target = cw.target(scope)
    pico = None
    traces, records, groups, pulse_log = [], [], [], []

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
        print(f"board: trig_mode={info['trig_mode']} out_of_place={info['out_of_place']} "
              f"probe_reps={info['probe_reps']}", flush=True)
        if info["trig_mode"] not in (2, 3, 4, 8):
            print(f"NOTE: CONJ_TRIG_MODE={info['trig_mode']} is not one of the line 2 modes "
                  f"(2, 3, 4, 8). The capture will still run.", flush=True)

        pico = rig.Pico(args.samples, args.timebase, args.pre_samples,
                        args.cha_range, args.chb_range, auto_trigger_ms=0)

        # The AC-coupled Ch A input needs a capture or two to settle; discard them.
        fixed_vec = rig.gen_vector(VECTOR, seed="fixed")
        for w in range(2):
            rig.load_vector(target, fixed_vec)
            rig.seed_board(target, FIXED_SEED16)
            pico.arm()
            target.flush()
            target.simpleserial_write("c", bytearray())
            pico.get(timeout_s=args.timeout)
            target.simpleserial_read("r", 1, timeout=int(args.timeout * 1000))
            print(f"warm-up capture {w + 1}/2 done", flush=True)

        t_start = time.time()
        for i in range(args.traces):
            if args.mode == "fixed":
                is_fixed = True
            elif args.mode == "random":
                is_fixed = False
            else:
                is_fixed = (i % 2 == 0)

            if is_fixed:
                rig.load_vector(target, fixed_vec)
                rig.seed_board(target, FIXED_SEED16)
            else:
                vec = rig.gen_vector(VECTOR, seed=f"{args.seed_prefix}{i}")
                rig.load_vector(target, vec)
                srng = random.Random(hash(args.seed_prefix) % 10**6 + i)
                rig.seed_board(target, bytes(srng.getrandbits(8) for _ in range(16)))

            pico.arm()
            target.flush()
            target.simpleserial_write("c", bytearray())
            power, trig, overflow = pico.get(timeout_s=args.timeout)
            ack = target.simpleserial_read("r", 1, timeout=int(args.timeout * 1000))
            if ack != bytearray([1]):
                raise RuntimeError(f"compute failed at trace {i}: {ack}")

            res = rig.parse_result(rig.read_outbuf(target))
            problems = rig.check_conjugate(res)
            if problems:
                raise RuntimeError(f"trace {i}: Conjugate check failed: {problems}")

            pulses = rig.trigger_pulses(trig)
            pulse_log.append(pulses)
            if args.crop:
                if not pulses:
                    raise RuntimeError(f"trace {i}: no trigger pulse found, cannot crop")
                s, e = pulses[0][0], pulses[-1][1]
                power = power[s:e]

            traces.append(power)
            records.append(res)
            groups.append(1 if is_fixed else 0)

            if i == 0 or (i + 1) % 10 == 0 or i == args.traces - 1:
                rate = (time.time() - t_start) / (i + 1)
                eta = rate * (args.traces - i - 1)
                c1 = res["pre_limbs"][1]
                print(f"{i:04d}/{args.traces} {'fixed' if is_fixed else 'rand '} "
                      f"pulses={len(pulses)} len={len(power)} "
                      f"pmV[{power.min():7.1f},{power.max():7.1f}] "
                      f"c1_limb0=0x{c1['limbs'][0]:08x} "
                      f"{'OVF' if overflow & 1 else ''} "
                      f"{rate:.1f}s/trace eta={eta / 60:.1f}min", flush=True)

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

    if args.crop:
        n = min(len(t) for t in traces)
        if max(len(t) for t in traces) != n:
            print(f"cropped lengths differ ({n} .. {max(len(t) for t in traces)}); "
                  f"truncating all to {n}", flush=True)
        trace_arr = np.vstack([t[:n] for t in traces])
    else:
        trace_arr = np.vstack(traces)

    labels, nlimbs = build_labels(records)
    # keep only the serialisable parts of each record
    slim = [{
        "pre": r["pre"], "post": r["post"], "pre_limbs": r["pre_limbs"],
        "norm": r["norm"], "skchall_norm": r["skchall_norm"],
    } for r in records]

    np.savez_compressed(
        out,
        traces=trace_arr,
        group=np.array(groups, dtype=np.uint8),
        interval_ns=np.float32(pico.interval_ns if pico else 0),
        samples_per_cycle=np.float32(pico.samples_per_cycle if pico else 0),
        trig_mode=np.uint8(info["trig_mode"]),
        out_of_place=np.uint8(info["out_of_place"]),
        nlimbs=np.uint8(nlimbs),
        pulses_json=json.dumps(pulse_log),
        records_json=json.dumps(slim),
        **labels,
    )
    print(f"\nsaved {trace_arr.shape} to {out}  "
          f"(fixed={sum(groups)}, random={len(groups) - sum(groups)})", flush=True)

    # ---- overview ----
    fig, ax = plt.subplots(3, 1, figsize=(14, 9))
    for i in range(min(10, len(traces))):
        ax[0].plot(trace_arr[i], lw=0.3, alpha=0.6)
    ax[0].set_title(f"Alg 3.18 line 2 Conjugate: first {min(10, len(traces))} of "
                    f"{len(traces)} traces (mode {args.mode})")
    ax[0].set_ylabel("mV")
    ax[0].grid(alpha=0.3)

    ax[1].plot(trig, color="green", lw=0.5)
    ax[1].set_title("Trigger (Ch B), last trace")
    ax[1].set_ylabel("mV")
    ax[1].grid(alpha=0.3)

    ax[2].plot(trace_arr.mean(axis=0), color="C0", lw=0.6)
    ax[2].set_title(f"Mean of {len(traces)} traces")
    ax[2].set_xlabel("sample")
    ax[2].set_ylabel("mV")
    ax[2].grid(alpha=0.3)

    fig.tight_layout()
    fig.savefig(plot, dpi=130)
    plt.close(fig)
    print(f"saved {plot}", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
