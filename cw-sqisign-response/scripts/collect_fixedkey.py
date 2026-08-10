#!/usr/bin/env python3
"""
Attack-faithful trace collection: ONE fixed secret key, fresh signature per trace.

This is the manuscript's real scenario and addresses RESEARCH_MEMORY findings #1
and #2. It is a NEW workflow added alongside collect_coords.py (which stays as a
serialization-profile collector); nothing existing is modified.

  Value source : host `dump-fixedkey N` -> one fixed Isk, fresh (Ichl from message,
                 Icom from a fresh commitment) per trace. The board input triple.
  Measured op  : command `c` runs the REAL response core; the trigger wraps
                 sample_finalize (materialization `basis*x` + normalization) --
                 RESP_TRIGGER_MODE=2, the attack target, NOT serialization.
                 (The artificial coord0 probe is now gated to mode 5 only, so the
                 mode-2 window is free of that contamination.)
  Label        : the true response is read back from the board (l/g) after capture.
  Cost         : ~12 s/trace (the core runs on the Cortex-M4); the trigger arrives
                 ~12 s after arming, so the scope waits indefinitely (autoTrigger=0).

  ~/.cwvenv/bin/python3 scripts/collect_fixedkey.py --traces 3000 \
      --out /home/ugrp/sqisign_traces/fixedkey01
"""
import argparse, ctypes, json, pathlib, subprocess, sys, time
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
CHUNK = 48
RANGES_MV = [10, 20, 50, 100, 200, 500, 1000, 2000]
AUTORANGE_HEADROOM = 1.4
WARMUP = 3


def run(cmd, **kw):
    print("+", " ".join(str(c) for c in cmd), flush=True)
    return subprocess.run(cmd, check=True, **kw)


def ss(t, cmd, data=b"", n=1, tmo=8000):
    t.simpleserial_write(cmd, bytearray(data))
    return t.simpleserial_read("r", n, timeout=tmo)


def load_fixedkey_triples(path):
    """Read `dump-fixedkey` output: '# trace i' comments + 54 hex tokens per trace."""
    toks = []
    for line in open(path):
        if line.startswith("#"):
            continue
        toks.extend(line.split())
    if len(toks) % 54:
        raise ValueError(f"{len(toks)} tokens not a multiple of 54")
    return [toks[i:i + 54] for i in range(0, len(toks), 54)]


def send_triple(target, tokens):
    data = ("\n".join(tokens) + "\n").encode()
    if ss(target, "x", b"", 1, 5000) != bytearray([1]):
        raise RuntimeError("reset failed")
    padded = data + b"\n" * ((-len(data)) % CHUNK)
    for off in range(0, len(padded), CHUNK):
        if ss(target, "a", padded[off:off + CHUNK], 1, 5000) != bytearray([1]):
            raise RuntimeError("triple append failed")


def read_board_response(target):
    r = ss(target, "l", b"", 2, 5000)
    outlen = (r[0] << 8) | r[1]
    out = bytearray()
    for idx in range((outlen + CHUNK - 1) // CHUNK):
        out.extend(ss(target, "g", bytes([(idx >> 8) & 0xFF, idx & 0xFF]), CHUNK, 5000))
    del out[outlen:]
    tk = out.decode("ascii", "replace").split()
    if len(tk) < 7 or tk[0] != "1":
        raise ValueError(f"bad response: {out[:60]!r}")
    v = [int(x, 16) for x in tk[1:7]]
    return {"lattice_content": v[0], "denom": v[1], "coord0": v[2],
            "coord1": v[3], "coord2": v[4], "coord3": v[5]}


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
        iv = ctypes.c_double(0); mx = ctypes.c_uint64(0)
        self._ok(ps.psospaGetTimebase(self.h, timebase, samples, ctypes.byref(iv),
                                      ctypes.byref(mx), 0))
        self.interval_ns = iv.value
        print(f"PicoScope {self.variant}: {self.res_bits}-bit, {self.interval_ns:.1f} ns/sample, "
              f"{samples:,} samples (~{samples*self.interval_ns/1e6:.2f} ms window)", flush=True)

    def _ok(self, st, what=""):
        if st != 0:
            raise RuntimeError(f"psospa error {st} {what}")

    def _open(self, bits):
        for want in (bits, 8):
            res = enums.PICO_DEVICE_RESOLUTION[f"PICO_DR_{want}BIT"]
            st = ps.psospaOpenUnit(ctypes.byref(self.h), None, res, None)
            if st == 0 and self.h.value > 0:
                self.res_bits = want
                buf = ctypes.create_string_buffer(128); rl = ctypes.c_int16()
                ps.psospaGetUnitInfo(self.h, buf, 128, ctypes.byref(rl), 3)
                self.variant = buf.value.decode(errors="replace")
                return res
            if self.h.value > 0:
                ps.psospaCloseUnit(self.h)
        raise RuntimeError("no psospa device")

    def set_ranges(self, cha_mv, chb_mv):
        self.cha_mv, self.chb_mv = cha_mv, chb_mv
        full = enums.PICO_BANDWIDTH_LIMITER["PICO_BW_FULL"]
        self._ok(ps.psospaSetChannelOn(self.h, self.A, enums.PICO_COUPLING["PICO_AC"],
                                       -cha_mv*1000000, cha_mv*1000000, 0, 0, full), "chA")
        self._ok(ps.psospaSetChannelOn(self.h, self.B, enums.PICO_COUPLING["PICO_DC"],
                                       -chb_mv*1000000, chb_mv*1000000, 0, 0, full), "chB")

    def set_trigger(self, threshold_mV=1000):
        thr = round(threshold_mV * self.maxADC.value / self.chb_mv)
        self.trig_thr_raw = thr
        # autoTrigger = 0 -> wait forever (finalize trigger arrives ~12 s after arming)
        self._ok(ps.psospaSetSimpleTrigger(self.h, 1, self.B, thr,
                                           enums.PICO_THRESHOLD_DIRECTION["PICO_RISING"], 0, 0),
                 "trigger")

    def arm(self):
        self.bufA = (ctypes.c_int16 * self.samples)()
        self.bufB = (ctypes.c_int16 * self.samples)()
        dt = enums.PICO_DATA_TYPE["PICO_INT16_T"]; raw = enums.PICO_RATIO_MODE["PICO_RATIO_MODE_RAW"]
        self._ok(ps.psospaSetDataBuffer(self.h, self.A, ctypes.byref(self.bufA), self.samples,
                 dt, 0, raw, enums.PICO_ACTION["PICO_CLEAR_ALL"] | enums.PICO_ACTION["PICO_ADD"]))
        self._ok(ps.psospaSetDataBuffer(self.h, self.B, ctypes.byref(self.bufB), self.samples,
                 dt, 0, raw, enums.PICO_ACTION["PICO_ADD"]))
        ti = ctypes.c_double(0)
        self._ok(ps.psospaRunBlock(self.h, self.pre, self.post, self.timebase,
                                   ctypes.byref(ti), 0, None, None), "runBlock")

    def get_raw(self, poll_s=0.05, timeout_s=60.0):
        ready = ctypes.c_int16(0); t0 = time.time()
        while ready.value == 0:
            ps.psospaIsReady(self.h, ctypes.byref(ready))
            if time.time() - t0 > timeout_s:
                raise TimeoutError("no trigger within timeout")
            time.sleep(poll_s)
        n = ctypes.c_uint64(self.samples); ov = ctypes.c_int16(0)
        self._ok(ps.psospaGetValues(self.h, 0, ctypes.byref(n), 1,
                 enums.PICO_RATIO_MODE["PICO_RATIO_MODE_RAW"], 0, ctypes.byref(ov)), "getValues")
        return (np.ctypeslib.as_array(self.bufA).copy(),
                np.ctypeslib.as_array(self.bufB).copy(), ov.value)

    def close(self):
        ps.psospaStop(self.h); ps.psospaCloseUnit(self.h)


def main():
    p = argparse.ArgumentParser(description="fixed-key attack-target collection (mode 2, cmd c)")
    p.add_argument("--traces", type=int, default=3000)
    p.add_argument("--out", required=True)
    p.add_argument("--triples", type=pathlib.Path, default=None,
                   help="fixedkey triples file; generated with host if absent")
    p.add_argument("--samples", type=int, default=2_700_000, help="~21.6 ms at 8 ns")
    p.add_argument("--timebase", type=int, default=8000)
    p.add_argument("--pre", type=int, default=5000)
    p.add_argument("--crop", type=int, default=2_500_000, help="samples stored from trigger rise")
    p.add_argument("--cha-mv", type=int, default=500)
    p.add_argument("--chb-mv", type=int, default=2000)
    p.add_argument("--resolution", type=int, default=10)
    p.add_argument("--platform", default="CW308_STM32F4")
    p.add_argument("--skip-flash", action="store_true")
    p.add_argument("--flush-every", type=int, default=50)
    a = p.parse_args()

    hexfile = FW_DIR / f"simpleserial-sqisign-response-{a.platform}.hex"
    if not hexfile.exists():
        sys.exit(f"firmware missing: {hexfile} (build RESP_TRIGGER_MODE=2)")
    if not HOST.exists():
        run([str(HOST_BUILD)])

    out = pathlib.Path(a.out).expanduser(); out.mkdir(parents=True, exist_ok=True)
    triples_path = a.triples or (out / "triples.txt")
    if not triples_path.exists():
        print(f"generating {a.traces} fixed-key triples ...", flush=True)
        with open(triples_path, "wb") as f:
            run([str(HOST), "dump-fixedkey", str(a.traces)], stdout=f)
    triples = load_fixedkey_triples(triples_path)
    if len(triples) < a.traces:
        sys.exit(f"triples file has {len(triples)} < {a.traces}")
    # confirm the fixed-key property
    sec = [tuple(t[:18]) for t in triples[:min(200, len(triples))]]
    if any(s != sec[0] for s in sec):
        sys.exit("secret ideal is NOT constant across triples -- wrong generator")
    print(f"loaded {len(triples):,} triples; secret ideal constant (fixed key) OK", flush=True)

    L = a.crop
    T = np.lib.format.open_memmap(out / "traces.npy", mode="w+", dtype=np.int16,
                                  shape=(a.traces, L))
    labels = (out / "labels.jsonl").open("w", buffering=1)

    scope = cw.scope(); target = cw.target(scope); pico = None
    n_done = 0; n_bad = 0; t_start = time.time()
    try:
        scope.default_setup()
        if not a.skip_flash:
            print(f"programming {hexfile.name} (mode 2) ...", flush=True)
            cw.program_target(scope, STM32FProgrammer, str(hexfile))
        scope.io.nrst = "low"; time.sleep(0.05); scope.io.nrst = "high_z"; time.sleep(0.5)
        target.flush()
        print("banner:", target.read(timeout=1000), flush=True)

        pico = Pico3000E(a.samples, a.timebase, a.pre, a.cha_mv, a.chb_mv, a.resolution)
        pico.set_trigger(1000)

        def capture_one(tok):
            send_triple(target, tok)
            pico.arm(); target.flush()
            target.simpleserial_write("c", bytearray())     # ~12 s real core
            pa, pb, ov = pico.get_raw()
            if target.simpleserial_read("r", 1, timeout=900000) != bytearray([1]):
                raise RuntimeError("compute ack failed")
            return pa, pb, ov

        # warm-up (AC settle) using the first triple
        for _ in range(WARMUP):
            capture_one(triples[0])
        print(f"warm-up: {WARMUP} captures done", flush=True)

        # autorange on a settled capture
        pa, _pb, _ov = capture_one(triples[0])
        peak = float(np.percentile(np.abs(pa), 99.99)) * (pico.cha_mv / pico.maxADC.value)
        best = next((r for r in RANGES_MV if peak * AUTORANGE_HEADROOM < r), pico.cha_mv)
        if best != pico.cha_mv:
            print(f"autorange: peak {peak:.1f} mV -> Ch A +/-{best} mV", flush=True)
            pico.set_ranges(best, a.chb_mv); pico.set_trigger(1000); capture_one(triples[0])

        meta = {"crop": L, "interval_ns": pico.interval_ns, "resolution_bits": pico.res_bits,
                "cha_mv": pico.cha_mv, "adc_max": pico.maxADC.value,
                "mv_per_code": pico.cha_mv / pico.maxADC.value, "capture_samples": a.samples,
                "timebase": a.timebase, "pre": a.pre, "source": "fixedkey_realcore",
                "trigger_mode": 2, "fixed_key": True, "platform": a.platform,
                "target_traces": a.traces}
        (out / "meta.json").write_text(json.dumps(meta, indent=2))
        print(json.dumps(meta), flush=True)

        for i in range(a.traces):
            pa, pb, ov = capture_one(triples[i])
            hi = pb > pico.trig_thr_raw
            idx = np.flatnonzero(hi)
            if len(idx) == 0:
                n_bad += 1; print(f"  trace {i}: no trigger, skip", flush=True); continue
            rise = int(idx[0]); fall = int(idx[-1])
            origin = max(0, rise - a.pre)
            if origin + L > a.samples:
                origin = a.samples - L
            rec = read_board_response(target)
            T[n_done] = pa[origin:origin + L]
            labels.write(json.dumps({
                "i": n_done, "trace_idx": i,
                "vals": {k: hex(v) for k, v in rec.items()},
                "trig_rise": rise - origin, "trig_fall": fall - origin,
                "trig_width": fall - rise, "ov": int(ov & 1)}) + "\n")
            n_done += 1
            if n_done % a.flush_every == 0:
                T.flush()
                el = time.time() - t_start; rate = n_done / el
                print(f"{n_done:5,}/{a.traces:,}  {rate*3600:.0f} tr/h  "
                      f"{n_done*L*2/1e9:.1f} GB  elapsed {el/3600:.2f} h  "
                      f"eta {(a.traces-n_done)/rate/3600:.1f} h  bad={n_bad}  "
                      f"trig_width={fall-rise}", flush=True)
    except KeyboardInterrupt:
        print("\ninterrupted -- flushing", flush=True)
    finally:
        if pico is not None:
            pico.close()
        try:
            target.dis()
        finally:
            scope.dis()
        T.flush(); labels.close()
        m = json.loads((out / "meta.json").read_text()) if (out / "meta.json").exists() else {}
        m["n_collected"] = n_done; m["n_bad"] = n_bad; m["elapsed_s"] = time.time() - t_start
        (out / "meta.json").write_text(json.dumps(m, indent=2))
        print(f"\ncollected {n_done:,} traces, bad={n_bad}, "
              f"{(time.time()-t_start)/3600:.2f} h -> {out}", flush=True)


if __name__ == "__main__":
    main()
