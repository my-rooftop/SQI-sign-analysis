#!/usr/bin/env python3
"""
Capture Algorithm 3.18 line 5, the formation of the response quaternion.

quat_lattice_sample_from_ball (lattice.c:159) rejection-samples a vector and then
forms alpha_rsp from it:

    ibz_mat_4x4_eval(&res->coord, &lat.basis, &x);   coord[i] = sum_j b[i][j]*x[j]
    ibz_copy(&res->denom, &lat.denom);
    quat_alg_normalize(res);

Two reasons not to trace the whole subroutine. It costs about 5.6 s per call, so
a campaign large enough to profile would take a day and a half; and the rejection
loop runs a data-dependent number of times, so traces do not align. The formation
runs exactly once after acceptance, so the host replays it with the lattice and
the accepted sample vector the real routine used, and the board traces that
alone -- the same arrangement the line 2 campaign used, and the executed code is
the vendored routine either way.

Unlike Conjugate, alpha_rsp is *formed* by a multiply-accumulate rather than
negated in place, which is the structure that leaks in the intersection
experiment.

Usage:
  ./conj_vector_host responses 20000 run1 > /tmp/responses.txt
  ~/.cwvenv/bin/python3 capture_resp_pico.py --records /tmp/responses.txt -n 20000
  ~/.cwvenv/bin/python3 capture_resp_pico.py --records /tmp/responses.txt \
      --mode tvla -n 1000 --out /tmp/sqisign_resp_tvla.npz
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

N_INPUT = 1 + 16 + 4      # denom, basis[4][4], x[4]
N_TRUTH = 1 + 4           # alpha_rsp denom, coord[0..3]


def read_records(path):
    lines = [l.rstrip("\n") for l in open(path) if l.strip() and not l.startswith("#")]
    i = 1 if lines and lines[0].startswith("count") else 0
    recs = []
    while i < len(lines):
        if not lines[i].startswith("r "):
            i += 1
            continue
        head = lines[i]
        body = lines[i + 1: i + 1 + N_INPUT]
        truth = [int(x, 16) for x in lines[i + 1 + N_INPUT: i + 1 + N_INPUT + N_TRUTH]]
        bounds = [int(t) for t in head.split()[1:]]
        recs.append({
            "text": "\n".join([head] + body) + "\n",
            "radius_bound": bounds[0],
            "bitlens": bounds[1:],
            "denom": truth[0],
            "coord": truth[1:5],
        })
        i += 1 + N_INPUT + N_TRUTH
    return recs


def limbs_of(value, bitlen):
    n = (bitlen + 31) // 32
    v = value & ((1 << (32 * n)) - 1)
    return [(v >> (32 * k)) & 0xFFFFFFFF for k in range(n)]


def main():
    p = argparse.ArgumentParser(description="Capture Alg 3.18 line 5 formation")
    p.add_argument("--records", required=True, help="file from `conj_vector_host responses`")
    p.add_argument("-n", "--traces", type=int, default=20000)
    p.add_argument("--mode", choices=("random", "tvla"), default="random",
                   help="random: every trace a fresh alpha_rsp, for SNR and "
                        "templates. tvla: interleave one fixed input with random "
                        "ones and record the group, for the Welch test.")
    p.add_argument("--seed", type=int, default=0,
                   help="seeds the fixed/random interleave in tvla mode")
    p.add_argument("--samples", type=int, default=60000)
    p.add_argument("--timebase", type=int, default=2)
    p.add_argument("--pre-samples", type=int, default=1000)
    p.add_argument("--cha-range", type=int, default=3)
    p.add_argument("--chb-range", type=int, default=8)
    p.add_argument("--skip-flash", action="store_true")
    p.add_argument("--out", default="/tmp/sqisign_resp.npz")
    p.add_argument("--plot", default="/tmp/sqisign_resp.png")
    args = p.parse_args()

    recs = read_records(args.records)
    if not recs:
        print(f"no records in {args.records}", file=sys.stderr)
        return 1
    n = min(args.traces, len(recs)) if args.mode == "random" else args.traces
    if args.mode == "tvla":
        # recs[0] is the fixed input; the random half draws from the rest, so a
        # short record file is fine here -- it just repeats random inputs.
        rng = np.random.default_rng(args.seed)
        group = rng.integers(0, 2, size=n).astype(np.uint8)
        order = [recs[0] if g else recs[1 + (i % max(1, len(recs) - 1))]
                 for i, g in enumerate(group)]
        print(f"{len(recs)} records available; tvla over {n} traces "
              f"({int(group.sum())} fixed, {n - int(group.sum())} random)", flush=True)
    else:
        group = None
        order = recs[:n]
        print(f"{len(recs)} records available, using {n}", flush=True)

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
        if info["trig_mode"] not in (15, 16, 17, 18):
            print(f"WARNING: CONJ_TRIG_MODE={info['trig_mode']} is not a line 5 "
                  f"formation mode (15-18)", flush=True)

        def load(rec):
            rig.load_vector(target, rec["text"].encode())
            if rig.ss_cmd(target, "b", b"", 1, timeout_ms=5000) != bytearray([1]):
                raise RuntimeError("loading the line 5 input failed")

        # Correctness gate: the board's alpha_rsp must equal the host's.
        load(recs[0])
        if rig.ss_cmd(target, "m", b"", 1, timeout_ms=10000) != bytearray([1]):
            raise RuntimeError("formation failed")
        rig.ss_cmd(target, "o", b"", 1, timeout_ms=5000)
        lines = [l for l in rig.read_outbuf(target).splitlines() if l.strip()]
        board_den = int(lines[4], 16)
        board_coord = [int(lines[5 + k], 16) for k in range(4)]
        if board_den != recs[0]["denom"] or board_coord != recs[0]["coord"]:
            raise RuntimeError("board alpha_rsp does not match the host")
        prof = {k: int(v) for k, v in (l.split("=") for l in lines[9:] if "=" in l)}
        print(f"board verified: alpha_rsp matches the host. "
              f"formation={prof.get('formation')} cycles "
              f"(eval {prof.get('eval')}, normalize {prof.get('normalize')})", flush=True)

        pico = rig.Pico(args.samples, args.timebase, args.pre_samples,
                        args.cha_range, args.chb_range, auto_trigger_ms=0)

        def one(rec, tries=5):
            load(rec)
            for _ in range(tries):
                pico.arm()
                target.flush()
                target.simpleserial_write("m", bytearray())
                try:
                    out = pico.get(timeout_s=1.0)
                except TimeoutError:
                    target.simpleserial_read("r", 1, timeout=5000)
                    continue
                if target.simpleserial_read("r", 1, timeout=5000) != bytearray([1]):
                    raise RuntimeError("board reported a formation failure")
                return out
            raise TimeoutError("no trigger")

        for _ in range(3):
            one(recs[0])

        t0 = time.time()
        for i in range(n):
            try:
                power, trig, overflow = one(order[i])
            except TimeoutError as e:
                print(f"\n{e} at trace {i}; keeping the {len(traces)} captured",
                      flush=True)
                break
            pulses = rig.trigger_pulses(trig)
            if not pulses:
                raise RuntimeError(f"trace {i}: no trigger pulse")
            traces.append(power[pulses[0][0]:pulses[-1][1]].astype(np.float32))
            if i == 0 or (i + 1) % 200 == 0 or i == n - 1:
                rate = (time.time() - t0) / (i + 1)
                print(f"{i + 1:6d}/{n}  len={len(traces[-1])}  pulses={len(pulses)}  "
                      f"{rate * 1000:.0f} ms/trace  eta={rate * (n - i - 1) / 60:.1f} min",
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

    m = min(len(t) for t in traces)
    T = np.vstack([t[:m] for t in traces])
    used = order[:len(traces)]
    extra = {} if group is None else {"group": group[:len(traces)]}

    # alpha_rsp's declared bound comes from the routine's exit, so derive the limb
    # count from the magnitude rather than assuming the line 2 layout.
    nl = max((abs(r["coord"][c]).bit_length() + 32) // 32 for r in used for c in range(4))
    labels = {}
    for c in range(4):
        vals = [r["coord"][c] for r in used]
        labels[f"c{c}_value_hex"] = np.array([f"{v & ((1 << 416) - 1):0104x}" for v in vals])
        labels[f"c{c}_magnitude_bits"] = np.array([abs(v).bit_length() for v in vals], np.uint16)
        labels[f"c{c}_sign"] = np.array([1 if v < 0 else 0 for v in vals], np.uint8)
        labels[f"c{c}_bitlen"] = np.array([abs(v).bit_length() + 1 for v in vals], np.uint16)
        labels[f"c{c}_nlimbs"] = np.full(len(used), nl, np.uint8)
        for k in range(nl):
            lv = np.array([limbs_of(r["coord"][c], 32 * nl)[k] for r in used], dtype=np.uint32)
            labels[f"c{c}_limb{k}"] = lv
            labels[f"c{c}_limb{k}_hw"] = np.array([int(v).bit_count() for v in lv], np.uint8)
        labels[f"c{c}_signext_limbs"] = np.array(
            [nl - (abs(r["coord"][c]).bit_length() + 32) // 32 for r in used], np.uint8)

    np.savez_compressed(args.out, traces=T, nlimbs=np.uint8(nl),
                        interval_ns=np.float32(pico.interval_ns if pico else 0),
                        samples_per_cycle=np.float32(pico.samples_per_cycle if pico else 0),
                        trig_mode=np.uint8(info["trig_mode"]), **extra, **labels)
    print(f"\nsaved {T.shape} to {args.out}", flush=True)

    fig, ax = plt.subplots(2, 1, figsize=(14, 7))
    for i in range(min(8, T.shape[0])):
        ax[0].plot(T[i], lw=0.3, alpha=0.6)
    ax[0].set_title(f"Alg 3.18 line 5 formation ({args.mode}), "
                    f"first {min(8, T.shape[0])} of {T.shape[0]}")
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
