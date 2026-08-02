#!/usr/bin/env python3
import argparse
import json
import pathlib
import random
import subprocess
import time

import chipwhisperer as cw
import numpy as np
from chipwhisperer.capture.api.programmers import STM32FProgrammer


APP = pathlib.Path(__file__).resolve().parent
PKG = APP.parent
FW_DIR = PKG / "firmware"
HOST = PKG / "host" / "response_split_host"
HOST_BUILD = PKG / "host" / "build_host_tool.sh"
HEX = FW_DIR / "simpleserial-sqisign-response-CW308_STM32F4.hex"
VECTOR = pathlib.Path("/tmp/sqisign_response_vec.txt")
CHUNK = 48


def run(cmd, **kwargs):
    print("+", " ".join(str(c) for c in cmd), flush=True)
    return subprocess.run(cmd, check=True, **kwargs)


def ss_cmd(target, cmd, data=b"", pay_len=1, timeout_ms=250000):
    target.simpleserial_write(cmd, bytearray(data))
    return target.simpleserial_read("r", pay_len, timeout=timeout_ms)


def parse_response(text):
    toks = text.split()
    if len(toks) < 7 or toks[0] != "1":
        raise ValueError(f"bad response: {text[:80]!r}")
    vals = [int(x, 16) for x in toks[1:7]]
    return {
        "lattice_content": vals[0],
        "denom": vals[1],
        "coord0": vals[2],
        "coord1": vals[3],
        "coord2": vals[4],
        "coord3": vals[5],
    }


def limbs32(x, n=16):
    y = abs(int(x))
    return [(y >> (32 * i)) & 0xFFFFFFFF for i in range(n)]


def hw32(x):
    return int(x & 0xFFFFFFFF).bit_count()


def read_board_response(target):
    r = ss_cmd(target, "l", b"", 2, timeout_ms=5000)
    outlen = (r[0] << 8) | r[1]
    out = bytearray()
    for idx in range((outlen + CHUNK - 1) // CHUNK):
        r = ss_cmd(target, "g", bytes([(idx >> 8) & 0xFF, idx & 0xFF]), CHUNK, timeout_ms=5000)
        out.extend(r)
    del out[outlen:]
    text = out.decode("ascii", errors="replace")
    return text, parse_response(text)


def setup_vector(target):
    with VECTOR.open("wb") as f:
        run([HOST, "gen"], stdout=f)
    vector = VECTOR.read_bytes()
    r = ss_cmd(target, "x", b"", 1, timeout_ms=5000)
    if r != bytearray([1]):
        raise RuntimeError(f"reset failed: {r}")
    padded = vector
    if len(padded) % CHUNK:
        padded += b"\n" * (CHUNK - (len(padded) % CHUNK))
    for off in range(0, len(padded), CHUNK):
        r = ss_cmd(target, "a", padded[off : off + CHUNK], 1, timeout_ms=5000)
        if r != bytearray([1]):
            raise RuntimeError(f"write failed at {off}: {r}")
    print(f"loaded vector bytes={len(vector)} chunks={len(padded) // CHUNK}", flush=True)


def quick_correlation(traces, labels):
    traces = traces.astype(np.float64)
    labels = labels.astype(np.float64)
    if labels.std() == 0:
        return -1, float("nan")
    labels -= labels.mean()
    traces -= traces.mean(axis=0)
    denom = np.sqrt((labels * labels).sum() * (traces * traces).sum(axis=0))
    denom[denom == 0] = np.inf
    corr = labels @ traces / denom
    idx = int(np.nanargmax(np.abs(corr)))
    return idx, float(corr[idx])


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-n", "--traces", type=int, default=20)
    parser.add_argument("--samples", type=int, default=80000)
    parser.add_argument("--offset", type=int, default=0)
    parser.add_argument("--timeout", type=int, default=180)
    parser.add_argument("--compute-each", action="store_true")
    parser.add_argument("--seed-each", action="store_true")
    parser.add_argument("--fixed-seed", action="store_true")
    parser.add_argument("--calibrate", action="store_true")
    parser.add_argument("--skip-build", action="store_true")
    parser.add_argument("--skip-flash", action="store_true")
    parser.add_argument("--out", default="/tmp/sqisign_response_leakage.npz")
    args = parser.parse_args()

    if not args.skip_build:
        run(["make", "PLATFORM=CW308_STM32F4", "RESP_TRIGGER_MODE=1"], cwd=FW_DIR)
        if not HOST.exists():
            run([HOST_BUILD])

    scope = cw.scope()
    target = cw.target(scope)
    traces = []
    records = []
    try:
        scope.default_setup()
        scope.adc.samples = args.samples
        scope.adc.offset = args.offset
        scope.adc.timeout = args.timeout
        if not args.skip_flash:
            cw.program_target(scope, STM32FProgrammer, str(HEX))
        time.sleep(0.2)
        target.flush()
        print("banner:", target.read(timeout=1000), flush=True)
        if not args.calibrate:
            setup_vector(target)

        if not args.compute_each and not args.calibrate:
            print("computing one response to populate last_resp", flush=True)
            t0 = time.time()
            r = ss_cmd(target, "c", b"", 1, timeout_ms=900000)
            print(f"initial compute: {r} elapsed={time.time() - t0:.3f}s", flush=True)
            if r != bytearray([1]):
                raise RuntimeError(f"initial compute failed: {r}")
            _text, initial = read_board_response(target)
            print(
                f"initial label: denom_bits={abs(initial['denom']).bit_length()} "
                f"coord0_lsw=0x{limbs32(initial['coord0'], 1)[0]:08x}",
                flush=True,
            )

        rng = random.Random(0x515149)
        for i in range(args.traces):
            if args.seed_each:
                seed_rng = random.Random(0 if args.fixed_seed else i)
                seed = bytes(seed_rng.getrandbits(8) for _ in range(16))
                r = ss_cmd(target, "n", seed, 1, timeout_ms=5000)
                if r != bytearray([1]):
                    raise RuntimeError(f"seed failed at trace {i}: {r}")
            if args.calibrate:
                value = rng.getrandbits(32)
                payload = bytes([(value >> 24) & 0xFF, (value >> 16) & 0xFF, (value >> 8) & 0xFF, value & 0xFF])
                r = ss_cmd(target, "q", payload, 1, timeout_ms=5000)
                if r != bytearray([1]):
                    raise RuntimeError(f"calibration set failed at trace {i}: {r}")
            scope.arm()
            target.simpleserial_write("c" if args.compute_each else "s", bytearray())
            timed_out = scope.capture()
            ack = target.simpleserial_read("r", 1, timeout=900000 if args.compute_each else 5000)
            if timed_out:
                raise RuntimeError(f"capture timed out at trace {i}")
            if ack != bytearray([1]):
                raise RuntimeError(f"compute failed at trace {i}: {ack}")
            wave = np.asarray(scope.get_last_trace(), dtype=np.float32)
            text, rec = read_board_response(target)
            rec["raw"] = text
            traces.append(wave)
            records.append(rec)
            print(
                f"{i:04d}: samples={len(wave)} denom_bits={abs(rec['denom']).bit_length()} "
                f"coord0_lsw=0x{limbs32(rec['coord0'], 1)[0]:08x}",
                flush=True,
            )

    finally:
        try:
            target.dis()
        finally:
            scope.dis()

    trace_arr = np.vstack(traces)
    names = ["denom", "coord0", "coord1", "coord2", "coord3"]
    labels = {}
    for name in names:
        vals = np.array([r[name] for r in records], dtype=object)
        labels[f"{name}_sign"] = np.array([1 if int(v) < 0 else 0 for v in vals], dtype=np.uint8)
        labels[f"{name}_bits"] = np.array([abs(int(v)).bit_length() for v in vals], dtype=np.uint16)
        labels[f"{name}_limb_count"] = np.array([(abs(int(v)).bit_length() + 31) // 32 for v in vals], dtype=np.uint8)
        for limb_idx in range(4):
            limb_vals = np.array([limbs32(v, limb_idx + 1)[limb_idx] for v in vals], dtype=np.uint32)
            labels[f"{name}_limb{limb_idx}"] = limb_vals
            labels[f"{name}_limb{limb_idx}_hw"] = np.array([hw32(v) for v in limb_vals], dtype=np.uint8)

    np.savez_compressed(
        args.out,
        traces=trace_arr,
        records_json=json.dumps(records),
        **labels,
    )
    print(f"saved {trace_arr.shape} to {args.out}", flush=True)

    if args.traces >= 5:
        for key in ("coord0_limb0_hw", "coord1_limb0_hw", "coord2_limb0_hw", "coord3_limb0_hw"):
            idx, corr = quick_correlation(trace_arr, labels[key])
            if idx < 0:
                print(f"quick corr {key}: skipped constant label", flush=True)
            else:
                print(f"quick corr {key}: sample={idx} corr={corr:+.4f}", flush=True)


if __name__ == "__main__":
    main()
