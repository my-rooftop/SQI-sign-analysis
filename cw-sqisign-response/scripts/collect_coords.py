#!/usr/bin/env python3
"""
Bulk coordinate-trace collection -- PicoScope 3000E (psospa) + CW-Lite + STM32F415
==================================================================================
Collects a profiling/analysis corpus for recovering coord0..3 of the SQIsign
response quaternion, from the mode-6 per-value serialization pulses.

Storage design (what makes a large corpus affordable):
  - Only the coord0..3 region is kept. The lattice_content pulse is ~105k samples
    (26% of the window) and is not a target, so it is captured but not stored.
  - Samples are stored as RAW int16 ADC codes, not float32 mV. Lossless, half the
    size; the mV scale factor lives in meta.json.
  - The trigger channel is never stored as a waveform -- only its six rise/fall
    sample indices per trace, rebased onto the crop origin.
  => 250,000 samples * 2 B = 0.5 MB/trace.

Ground truth is written incrementally to labels.jsonl (one JSON object per trace),
so an interrupted run keeps everything captured up to that point.

Value source:
  --known (default)  cmd 'k' injects a response with realistic bit-lengths. Fast.
  --real             cmd 'c' runs the true core per trace (~12 s) for a test set.

  ~/.cwvenv/bin/python3 scripts/collect_coords.py --traces 200000 --out ~/sqisign_traces/prof01
"""

import argparse
import ctypes
import json
import pathlib
import random
import subprocess
import sys
import time

import numpy as np

import chipwhisperer as cw
from chipwhisperer.capture.api.programmers import STM32FProgrammer
from picosdk.psospa import psospa as ps
from picosdk.PicoDeviceEnums import picoEnum as enums

APP = pathlib.Path(__file__).resolve().parent
PKG = APP.parent
FW_DIR = PKG / "firmware"
HOST = PKG / "host" / "response_split_host"
HOST_BUILD = PKG / "host" / "build_host_tool.sh"
VECTOR = pathlib.Path("/tmp/sqisign_response_vec.txt")
CHUNK = 48

# coord0 rise -> coord3 fall measured at 226,989 samples; worst-case coordinate
# sizes (328+327+202+201 bit) scale that to ~232.5k. 250k leaves ~7.5% headroom.
CROP_LEN = 250_000
PRE_MARGIN = 1_000          # samples kept before the coord0 rising edge
RANGES_MV = [10, 20, 50, 100, 200, 500, 1000, 2000]
AUTORANGE_HEADROOM = 1.4
WARMUP_CAPTURES = 5
BL = {"lattice_content": (516, 528), "coord0": (316, 328), "coord1": (316, 327),
      "coord2": (188, 202), "coord3": (188, 201)}


def run(cmd, **kw):
    print("+", " ".join(str(c) for c in cmd), flush=True)
    return subprocess.run(cmd, check=True, **kw)


def ss(target, cmd, data=b"", n=1, timeout_ms=8000):
    target.simpleserial_write(cmd, bytearray(data))
    return target.simpleserial_read("r", n, timeout=timeout_ms)


def gen_response(rng):
    """Synthetic response: right bit-lengths, but uniform-random contents, always
    positive, denom=1. A verification against 500 real responses showed this
    matches real coords on bit-length and LS-limb Hamming weight, but NOT on sign
    (real coords are ~50% negative) or denom (~2/3 are denom=2). Prefer --values
    with real precomputed responses; keep this only as a no-host-tool fallback."""
    def val(lo, hi):
        b = rng.randint(lo, hi)
        return (1 << (b - 1)) | rng.getrandbits(b - 1)
    return {"lattice_content": val(*BL["lattice_content"]), "denom": 1,
            "coord0": val(*BL["coord0"]), "coord1": val(*BL["coord1"]),
            "coord2": val(*BL["coord2"]), "coord3": val(*BL["coord3"])}


def load_real_responses(path):
    """Read precomputed real responses from `host/response_split_host dump-responses`.
    One line per response: content denom coord0 coord1 coord2 coord3 (signed hex)."""
    out = []
    with open(path) as f:
        for line in f:
            t = line.split()
            if len(t) != 6:
                continue
            out.append({"lattice_content": int(t[0], 16), "denom": int(t[1], 16),
                        "coord0": int(t[2], 16), "coord1": int(t[3], 16),
                        "coord2": int(t[4], 16), "coord3": int(t[5], 16)})
    return out


def load_known(target, r):
    txt = (f"{r['lattice_content']:x}\n{r['denom']:x}\n{r['coord0']:x}\n"
           f"{r['coord1']:x}\n{r['coord2']:x}\n{r['coord3']:x}\n").encode()
    if ss(target, "x", b"", 1, 5000) != bytearray([1]):
        raise RuntimeError("reset failed")
    padded = txt + b"\n" * ((-len(txt)) % CHUNK)
    for off in range(0, len(padded), CHUNK):
        if ss(target, "a", padded[off:off + CHUNK], 1, 5000) != bytearray([1]):
            raise RuntimeError("append failed")
    if ss(target, "k", b"", 1, 5000) != bytearray([1]):
        raise RuntimeError("set-known failed")


def read_board_response(target):
    r = ss(target, "l", b"", 2, 5000)
    outlen = (r[0] << 8) | r[1]
    out = bytearray()
    for idx in range((outlen + CHUNK - 1) // CHUNK):
        out.extend(ss(target, "g", bytes([(idx >> 8) & 0xFF, idx & 0xFF]), CHUNK, 5000))
    del out[outlen:]
    toks = out.decode("ascii", "replace").split()
    if len(toks) < 7 or toks[0] != "1":
        raise ValueError(f"bad response: {out[:60]!r}")
    v = [int(x, 16) for x in toks[1:7]]
    return {"lattice_content": v[0], "denom": v[1], "coord0": v[2],
            "coord1": v[3], "coord2": v[4], "coord3": v[5]}


def setup_vector(target):
    with VECTOR.open("wb") as f:
        run([str(HOST), "gen"], stdout=f)
    vec = VECTOR.read_bytes()
    if ss(target, "x", b"", 1, 5000) != bytearray([1]):
        raise RuntimeError("reset failed")
    padded = vec + b"\n" * ((-len(vec)) % CHUNK)
    for off in range(0, len(padded), CHUNK):
        if ss(target, "a", padded[off:off + CHUNK], 1, 5000) != bytearray([1]):
            raise RuntimeError("vector append failed")


class Pico3000E:
    A = enums.PICO_CHANNEL["PICO_CHANNEL_A"]
    B = enums.PICO_CHANNEL["PICO_CHANNEL_B"]

    def __init__(self, samples, timebase, pre, cha_mv, chb_mv, resolution=10):
        self.h = ctypes.c_int16()
        self.samples, self.timebase, self.pre = samples, timebase, pre
        self.post = samples - pre
        self.cha_mv, self.chb_mv = cha_mv, chb_mv
        self.res = self._open(resolution)
        self.minADC, self.maxADC = ctypes.c_int16(), ctypes.c_int16()
        self._ok(ps.psospaGetAdcLimits(self.h, self.res,
                                       ctypes.byref(self.minADC), ctypes.byref(self.maxADC)))
        self.set_ranges(cha_mv, chb_mv)
        for ch in (2, 3):
            ps.psospaSetChannelOff(self.h, ch)
        interval = ctypes.c_double(0)
        max_samples = ctypes.c_uint64(0)
        self._ok(ps.psospaGetTimebase(self.h, timebase, samples,
                                      ctypes.byref(interval), ctypes.byref(max_samples), 0))
        self.interval_ns = interval.value
        print(f"PicoScope {self.variant}: {self.res_bits}-bit, {self.interval_ns:.2f} ns/sample, "
              f"{samples:,} samples (~{samples * self.interval_ns / 1e6:.3f} ms)", flush=True)

    def _ok(self, status, what=""):
        if status != 0:
            raise RuntimeError(f"psospa error {status} {what}")

    def _open(self, bits):
        for want in (bits, 8):
            res = enums.PICO_DEVICE_RESOLUTION[f"PICO_DR_{want}BIT"]
            st = ps.psospaOpenUnit(ctypes.byref(self.h), None, res, None)
            if st == 0 and self.h.value > 0:
                self.res_bits = want
                buf = ctypes.create_string_buffer(128)
                rl = ctypes.c_int16()
                ps.psospaGetUnitInfo(self.h, buf, 128, ctypes.byref(rl), 3)
                self.variant = buf.value.decode(errors="replace")
                return res
            if self.h.value > 0:
                ps.psospaCloseUnit(self.h)
            print(f"  {want}-bit rejected (status {st}), falling back", flush=True)
        raise RuntimeError("no psospa device could be opened")

    def set_ranges(self, cha_mv, chb_mv):
        self.cha_mv, self.chb_mv = cha_mv, chb_mv
        full = enums.PICO_BANDWIDTH_LIMITER["PICO_BW_FULL"]
        self._ok(ps.psospaSetChannelOn(self.h, self.A, enums.PICO_COUPLING["PICO_AC"],
                                       -cha_mv * 1000000, cha_mv * 1000000, 0, 0, full), "chA")
        self._ok(ps.psospaSetChannelOn(self.h, self.B, enums.PICO_COUPLING["PICO_DC"],
                                       -chb_mv * 1000000, chb_mv * 1000000, 0, 0, full), "chB")

    def set_trigger(self, threshold_mV=1000, auto_us=5_000_000):
        thr = round(threshold_mV * self.maxADC.value / self.chb_mv)
        self.trig_thr_raw = thr
        self._ok(ps.psospaSetSimpleTrigger(self.h, 1, self.B, thr,
                                           enums.PICO_THRESHOLD_DIRECTION["PICO_RISING"],
                                           0, auto_us), "trigger")

    def arm(self):
        self.bufA = (ctypes.c_int16 * self.samples)()
        self.bufB = (ctypes.c_int16 * self.samples)()
        dtype = enums.PICO_DATA_TYPE["PICO_INT16_T"]
        raw = enums.PICO_RATIO_MODE["PICO_RATIO_MODE_RAW"]
        self._ok(ps.psospaSetDataBuffer(self.h, self.A, ctypes.byref(self.bufA), self.samples,
                                        dtype, 0, raw,
                                        enums.PICO_ACTION["PICO_CLEAR_ALL"] |
                                        enums.PICO_ACTION["PICO_ADD"]), "bufA")
        self._ok(ps.psospaSetDataBuffer(self.h, self.B, ctypes.byref(self.bufB), self.samples,
                                        dtype, 0, raw, enums.PICO_ACTION["PICO_ADD"]), "bufB")
        t_ind = ctypes.c_double(0)
        self._ok(ps.psospaRunBlock(self.h, self.pre, self.post, self.timebase,
                                   ctypes.byref(t_ind), 0, None, None), "runBlock")

    def get_raw(self, timeout_s=30.0):
        """Returns raw int16 ADC codes -- no mV conversion, that is a storage cost."""
        ready = ctypes.c_int16(0)
        t0 = time.time()
        while ready.value == 0:
            ps.psospaIsReady(self.h, ctypes.byref(ready))
            if time.time() - t0 > timeout_s:
                raise TimeoutError("capture did not complete (no trigger?)")
        n = ctypes.c_uint64(self.samples)
        overflow = ctypes.c_int16(0)
        self._ok(ps.psospaGetValues(self.h, 0, ctypes.byref(n), 1,
                                    enums.PICO_RATIO_MODE["PICO_RATIO_MODE_RAW"], 0,
                                    ctypes.byref(overflow)), "getValues")
        return (np.ctypeslib.as_array(self.bufA).copy(),
                np.ctypeslib.as_array(self.bufB).copy(), overflow.value)

    def close(self):
        ps.psospaStop(self.h)
        ps.psospaCloseUnit(self.h)


def pulse_edges_raw(trig_raw, thr_raw):
    hi = trig_raw > thr_raw
    rise = np.flatnonzero((~hi[:-1]) & hi[1:]) + 1
    fall = np.flatnonzero(hi[:-1] & (~hi[1:])) + 1
    return rise, fall


def main():
    p = argparse.ArgumentParser(description="bulk coord0..3 trace collection")
    p.add_argument("-n", "--traces", type=int, default=1000)
    p.add_argument("--out", required=True, help="output directory")
    p.add_argument("--real", action="store_true",
                   help="run the true core per trace (~12 s) instead of cmd 'k'")
    p.add_argument("--values", type=pathlib.Path, default=None,
                   help="file of precomputed REAL responses (host dump-responses); "
                        "injected via 'k' -- real distribution at the fast serialize speed")
    p.add_argument("--samples", type=int, default=400000)
    p.add_argument("--timebase", type=int, default=8000)
    p.add_argument("--pre", type=int, default=2000)
    p.add_argument("--cha-mv", type=int, default=500, help="autorange starts here")
    p.add_argument("--chb-mv", type=int, default=2000)
    p.add_argument("--resolution", type=int, default=10, choices=[8, 10, 12])
    p.add_argument("--crop-len", type=int, default=CROP_LEN)
    p.add_argument("--verify-every", type=int, default=200,
                   help="read the response back from the board every N traces; the "
                        "values are already known in --known mode, so this is a "
                        "periodic integrity check rather than per-trace overhead")
    p.add_argument("--flush-every", type=int, default=500)
    p.add_argument("--max-hours", type=float, default=0.0, help="0 = no limit")
    p.add_argument("--seed", type=lambda s: int(s, 0), default=0x5EED)
    p.add_argument("--platform", default="CW308_STM32F4")
    p.add_argument("--skip-flash", action="store_true")
    a = p.parse_args()

    hexfile = FW_DIR / f"simpleserial-sqisign-response-{a.platform}.hex"
    if not hexfile.exists():
        sys.exit(f"firmware missing: {hexfile}")
    if a.real and not HOST.exists():
        run([str(HOST_BUILD)])

    out = pathlib.Path(a.out).expanduser()
    out.mkdir(parents=True, exist_ok=True)
    L = a.crop_len
    need_gb = a.traces * L * 2 / 1e9
    print(f"output   : {out}\n"
          f"budget   : {a.traces:,} traces x {L:,} samples x 2 B = {need_gb:.1f} GB", flush=True)

    real_vals = None
    if a.values:
        real_vals = load_real_responses(a.values)
        if len(real_vals) < a.traces:
            sys.exit(f"--values has {len(real_vals):,} responses but --traces is "
                     f"{a.traces:,}; generate more with dump-responses")
        print(f"values   : {len(real_vals):,} real responses from {a.values}", flush=True)

    T = np.lib.format.open_memmap(out / "traces.npy", mode="w+",
                                  dtype=np.int16, shape=(a.traces, L))
    labels = (out / "labels.jsonl").open("w", buffering=1)

    scope = cw.scope()
    target = cw.target(scope)
    pico = None
    rng = random.Random(a.seed)
    n_done = 0
    n_bad = 0
    t_start = time.time()
    try:
        scope.default_setup()
        if not a.skip_flash:
            print(f"programming {hexfile.name} ...", flush=True)
            cw.program_target(scope, STM32FProgrammer, str(hexfile))
        scope.io.nrst = "low"
        time.sleep(0.05)
        scope.io.nrst = "high_z"
        time.sleep(0.5)
        target.flush()
        print("banner:", target.read(timeout=1000), flush=True)

        if a.real:
            setup_vector(target)

        pico = Pico3000E(a.samples, a.timebase, a.pre, a.cha_mv, a.chb_mv, a.resolution)
        pico.set_trigger(1000, auto_us=5_000_000)

        def source_new_value(i):
            """Put a fresh response into last_resp; return ground truth."""
            if a.real:
                if ss(target, "c", b"", 1, 900000) != bytearray([1]):
                    raise RuntimeError("compute failed")
                return read_board_response(target)
            r = real_vals[i] if real_vals is not None else gen_response(rng)
            load_known(target, r)
            return r

        def capture():
            pico.arm()
            target.flush()
            target.simpleserial_write("s", bytearray())
            pa, pb, ov = pico.get_raw()
            if target.simpleserial_read("r", 1, timeout=8000) != bytearray([1]):
                raise RuntimeError("serialize ack failed")
            return pa, pb, ov

        source_new_value(0)
        for _ in range(WARMUP_CAPTURES):
            capture()
        print(f"warm-up: {WARMUP_CAPTURES} captures discarded (AC settle)", flush=True)

        # autorange on settled captures, percentile not max
        scale = pico.cha_mv / pico.maxADC.value
        peaks = []
        for _ in range(2):
            pa, _pb, _ov = capture()
            peaks.append(float(np.percentile(np.abs(pa), 99.99)) * scale)
        peak = max(peaks)
        best = next((r for r in RANGES_MV if peak * AUTORANGE_HEADROOM < r), pico.cha_mv)
        if best != pico.cha_mv:
            print(f"autorange: peak(p99.99) {peak:.1f} mV -> Ch A +/-{best} mV", flush=True)
            pico.set_ranges(best, a.chb_mv)
            pico.set_trigger(1000, auto_us=5_000_000)
            for _ in range(2):
                capture()
        else:
            print(f"autorange: keeping Ch A +/-{pico.cha_mv} mV (peak {peak:.1f} mV)", flush=True)

        meta = {"crop_len": L, "pre_margin": PRE_MARGIN, "interval_ns": pico.interval_ns,
                "resolution_bits": pico.res_bits, "cha_mv": pico.cha_mv,
                "chb_mv": pico.chb_mv, "adc_max": pico.maxADC.value,
                "mv_per_code": pico.cha_mv / pico.maxADC.value,
                "capture_samples": a.samples, "timebase": a.timebase, "pre": a.pre,
                "source": ("real_core" if a.real else
                           "real_values" if real_vals is not None else "synthetic"),
                "values_file": str(a.values) if a.values else None, "seed": a.seed,
                "platform": a.platform, "trigger_mode": 6, "target_traces": a.traces}
        (out / "meta.json").write_text(json.dumps(meta, indent=2))
        print(json.dumps(meta), flush=True)

        deadline = t_start + a.max_hours * 3600 if a.max_hours > 0 else None
        for i in range(a.traces):
            if deadline and time.time() > deadline:
                print(f"stopping: --max-hours reached at trace {i}", flush=True)
                break
            truth = source_new_value(i)
            pa, pb, ov = capture()
            rise, fall = pulse_edges_raw(pb, pico.trig_thr_raw)
            if len(rise) < 6 or len(fall) < 6:
                n_bad += 1
                print(f"  trace {i}: only {len(rise)} pulses, skipped", flush=True)
                continue
            origin = int(rise[2]) - PRE_MARGIN            # coord0 rise, minus margin
            if origin < 0 or origin + L > a.samples:
                n_bad += 1
                print(f"  trace {i}: crop {origin}..{origin+L} outside window, skipped", flush=True)
                continue
            if int(fall[5]) > origin + L:
                n_bad += 1
                print(f"  trace {i}: coord3 ends past crop ({fall[5]-origin} > {L}), skipped",
                      flush=True)
                continue
            T[n_done] = pa[origin:origin + L]
            rec = {"i": n_done,
                   "vals": {k: hex(v) for k, v in truth.items()},
                   "rise": [int(x) - origin for x in rise[:6]],
                   "fall": [int(x) - origin for x in fall[:6]],
                   "origin": origin, "ov": int(ov & 1)}
            if a.verify_every and n_done % a.verify_every == 0:
                back = read_board_response(target)
                rec["verified"] = all(back[k] == truth[k] for k in truth)
                if not rec["verified"]:
                    raise RuntimeError(f"read-back mismatch at trace {n_done}")
            labels.write(json.dumps(rec) + "\n")
            n_done += 1

            if n_done % a.flush_every == 0:
                T.flush()
                el = time.time() - t_start
                rate = n_done / el
                gb = n_done * L * 2 / 1e9
                eta = (a.traces - n_done) / rate / 3600
                print(f"{n_done:7,}/{a.traces:,}  {rate:.2f} tr/s  {gb:.1f} GB  "
                      f"elapsed {el/3600:.2f} h  eta {eta:.1f} h  bad={n_bad}", flush=True)
    except KeyboardInterrupt:
        print("\ninterrupted -- flushing what we have", flush=True)
    finally:
        if pico is not None:
            pico.close()
        try:
            target.dis()
        finally:
            scope.dis()
        T.flush()
        labels.close()
        meta_path = out / "meta.json"
        m = json.loads(meta_path.read_text()) if meta_path.exists() else {}
        m["n_collected"] = n_done
        m["n_bad"] = n_bad
        m["elapsed_s"] = time.time() - t_start
        meta_path.write_text(json.dumps(m, indent=2))
        print(f"\ncollected {n_done:,} traces ({n_done*L*2/1e9:.1f} GB), bad={n_bad}, "
              f"{(time.time()-t_start)/3600:.2f} h -> {out}", flush=True)
        if n_done < a.traces:
            print(f"NOTE: traces.npy is sized for {a.traces:,}; rows {n_done:,}.. are unwritten. "
                  f"Use meta.json's n_collected when loading.", flush=True)


if __name__ == "__main__":
    main()
