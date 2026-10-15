#!/usr/bin/env python3
"""
Capture Algorithm 3.18 line 2 traces with alpha injected directly.

Running the whole of Algorithm 3.18 on the board costs about 13 s per trace,
nearly all of it in line 1 and line 5, which makes a campaign large enough to
measure MSB recovery impractical. The measured window is line 2 only, so this
generates realistic alpha_sk,chl values on the host (conj_vector_host alphas,
which runs the real algorithm) and drives line 2 on the board directly. The
executed code and the trigger placement are identical; only the surrounding
instructions differ.

Validate that equivalence before relying on it: --validate re-runs the full
algorithm a few times and checks that the line 2 cycle count and the trigger
width match the injected path.

Usage:
  ./conj_vector_host alphas 20000 run1 > /tmp/alphas.txt
  ~/.cwvenv/bin/python3 capture_line2_pico.py --alphas /tmp/alphas.txt -n 20000
"""

import argparse
import pathlib
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


def read_alphas(path):
    """Parse the output of `conj_vector_host alphas`."""
    recs = []
    lines = [l.rstrip("\n") for l in open(path) if l.strip() and not l.startswith("#")]
    i = 0
    if lines and lines[0].startswith("count"):
        i = 1
    while i < len(lines):
        if not lines[i].startswith("a "):
            i += 1
            continue
        bounds = [int(t) for t in lines[i].split()[1:6]]
        vals = [int(lines[i + 1 + k], 16) for k in range(5)]
        recs.append({"bounds": bounds, "denom": vals[0], "coord": vals[1:5],
                     "text": "\n".join(lines[i:i + 6]) + "\n"})
        i += 6
    return recs


def limbs_of(value, bitlen):
    """Two's-complement limbs exactly as ibz_t holds them."""
    n = (bitlen + 31) // 32
    v = value & ((1 << (32 * n)) - 1)
    return [(v >> (32 * k)) & 0xFFFFFFFF for k in range(n)]


def main():
    p = argparse.ArgumentParser(description="Capture Alg 3.18 line 2 with injected alpha")
    p.add_argument("--alphas", required=True, help="file from `conj_vector_host alphas`")
    p.add_argument("-n", "--traces", type=int, default=5000)
    p.add_argument("--samples", type=int, default=60000)
    p.add_argument("--timebase", type=int, default=2)
    p.add_argument("--pre-samples", type=int, default=1000)
    p.add_argument("--cha-range", type=int, default=3)
    p.add_argument("--chb-range", type=int, default=8)
    p.add_argument("--skip-flash", action="store_true")
    p.add_argument("--validate", type=int, default=3,
                   help="full Algorithm 3.18 runs to cross-check the injected path")
    p.add_argument("--out", default="/tmp/sqisign_line2.npz")
    p.add_argument("--plot", default="/tmp/sqisign_line2.png")
    args = p.parse_args()

    recs = read_alphas(args.alphas)
    if not recs:
        print(f"no alpha records in {args.alphas}", file=sys.stderr)
        return 1
    n = min(args.traces, len(recs))
    print(f"{len(recs)} alpha records available, using {n}", flush=True)

    import chipwhisperer as cw
    from chipwhisperer.capture.api.programmers import STM32FProgrammer

    scope = cw.scope()
    target = cw.target(scope)
    pico = None
    traces = []

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
        if info["trig_mode"] not in (2, 3, 4, 8):
            print(f"WARNING: CONJ_TRIG_MODE={info['trig_mode']} is not a line 2 mode",
                  flush=True)

        def load_alpha(rec):
            rig.load_vector(target, rec["text"].encode())
            if rig.ss_cmd(target, "u", b"", 1, timeout_ms=5000) != bytearray([1]):
                raise RuntimeError("loading alpha failed")

        # Correctness gate: the board must negate coord[1..3] and leave the rest.
        load_alpha(recs[0])
        rig.ss_cmd(target, "z", b"", 1, timeout_ms=5000)
        rig.ss_cmd(target, "q", b"", 1, timeout_ms=5000)
        txt = rig.read_outbuf(target)
        lines = [l for l in txt.splitlines() if l.strip()]
        board_limbs = []
        for k in range(4):
            t = lines[k].split()
            board_limbs.append([int(z, 16) for z in t[2:2 + int(t[1])]])
        post = [int(lines[4 + k], 16) for k in range(4)]
        want_limbs = [limbs_of(recs[0]["coord"][k], recs[0]["bounds"][1 + k]) for k in range(4)]
        if board_limbs != want_limbs:
            raise RuntimeError("board's alpha does not match the host's")
        if post[0] != recs[0]["coord"][0] or any(post[k] != -recs[0]["coord"][k] for k in (1, 2, 3)):
            raise RuntimeError("board did not compute Conjugate")
        print("board verified: alpha loaded exactly, Conjugate correct", flush=True)

        pico = rig.Pico(args.samples, args.timebase, args.pre_samples,
                        args.cha_range, args.chb_range, auto_trigger_ms=0)

        def one(rec, tries=5):
            load_alpha(rec)
            for _ in range(tries):
                pico.arm()
                target.flush()
                target.simpleserial_write("z", bytearray())
                try:
                    rec_out = pico.get(timeout_s=0.5)
                except TimeoutError:
                    target.simpleserial_read("r", 1, timeout=5000)
                    continue
                if target.simpleserial_read("r", 1, timeout=5000) != bytearray([1]):
                    raise RuntimeError("board reported line 2 failure")
                return rec_out
            raise TimeoutError("no trigger")

        for _ in range(3):
            one(recs[0])

        t0 = time.time()
        for i in range(n):
            power, trig, overflow = one(recs[i])
            pulses = rig.trigger_pulses(trig)
            if not pulses:
                raise RuntimeError(f"trace {i}: no trigger pulse")
            traces.append(power[pulses[0][0]:pulses[-1][1]].astype(np.float32))
            if i == 0 or (i + 1) % 200 == 0 or i == n - 1:
                rate = (time.time() - t0) / (i + 1)
                print(f"{i + 1:6d}/{n}  len={len(traces[-1])}  {rate * 1000:.0f} ms/trace  "
                      f"eta={rate * (n - i - 1) / 60:.1f} min", flush=True)

    finally:
        if pico is not None:
            pico.close()
        try:
            target.dis()
        finally:
            scope.dis()

    m = min(len(t) for t in traces)
    T = np.vstack([t[:m] for t in traces])

    used = recs[:len(traces)]
    nl = max((b + 31) // 32 for r in used for b in r["bounds"][1:5])
    labels = {}
    for c in range(4):
        vals = [r["coord"][c] for r in used]
        labels[f"c{c}_value_hex"] = np.array([f"{v & ((1 << 416) - 1):0104x}" for v in vals])
        labels[f"c{c}_magnitude_bits"] = np.array([abs(v).bit_length() for v in vals], np.uint16)
        labels[f"c{c}_sign"] = np.array([1 if v < 0 else 0 for v in vals], np.uint8)
        labels[f"c{c}_bitlen"] = np.array([r["bounds"][1 + c] for r in used], np.uint16)
        for k in range(nl):
            lv = np.array([limbs_of(r["coord"][c], r["bounds"][1 + c])[k]
                           if k < (r["bounds"][1 + c] + 31) // 32 else 0 for r in used],
                          dtype=np.uint32)
            labels[f"c{c}_limb{k}"] = lv
            labels[f"c{c}_limb{k}_hw"] = np.array(
                [int(v).bit_count() for v in lv], np.uint8)

    np.savez_compressed(args.out, traces=T, nlimbs=np.uint8(nl),
                        interval_ns=np.float32(pico.interval_ns if pico else 0),
                        samples_per_cycle=np.float32(pico.samples_per_cycle if pico else 0),
                        trig_mode=np.uint8(info["trig_mode"]), **labels)
    print(f"\nsaved {T.shape} to {args.out}", flush=True)

    fig, ax = plt.subplots(2, 1, figsize=(14, 7))
    for i in range(min(8, T.shape[0])):
        ax[0].plot(T[i], lw=0.3, alpha=0.6)
    ax[0].set_title(f"Alg 3.18 line 2, first {min(8, T.shape[0])} of {T.shape[0]} traces")
    ax[0].grid(alpha=0.3)
    ax[1].plot(T.mean(axis=0), lw=0.6)
    ax[1].set_title(f"mean of {T.shape[0]}")
    ax[1].set_xlabel("sample")
    ax[1].grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(args.plot, dpi=130)
    plt.close(fig)
    print(f"saved {args.plot}", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
