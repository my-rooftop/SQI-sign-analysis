#!/usr/bin/env python3
"""
Large-scale PROFILING with realistic known values (firmware mode 6 + 'v' command).
No 12 s response core: we synthesize responses whose coordinates match the real
bit-length distribution, load them via 'v', and serialize under 6 per-value pulses.
Ground truth = exactly the values we set -> ideal for training per-limb templates
that then transfer to the real 300-trace set.

  ~/.cwvenv/bin/python3 exp_profile_known.py --traces 5000 --samples 400000 --timebase 3
"""
import argparse, ctypes, json, pathlib, random, sys, time
import numpy as np
import chipwhisperer as cw
from chipwhisperer.capture.api.programmers import STM32FProgrammer
from picosdk.ps3000a import ps3000a as ps
from picosdk.functions import adc2mV, assert_pico_ok

APP = pathlib.Path(__file__).resolve().parent
HEX = APP / "simpleserial-sqisign-response-CW308_STM32F4.hex"
CHUNK = 48
# real lvl1 coordinate bit-length ranges (from RESPONSE_FLOW + measured)
BL = {"content": (516, 528), "coord0": (316, 328), "coord1": (316, 327),
      "coord2": (188, 202), "coord3": (188, 201)}


def ss(t, cmd, data=b"", n=1, tmo=8000):
    t.simpleserial_write(cmd, bytearray(data)); return t.simpleserial_read("r", n, timeout=tmo)


def gen_response(rng):
    def val(lo, hi):
        b = rng.randint(lo, hi); return (1 << (b - 1)) | rng.getrandbits(b - 1)
    return dict(content=val(*BL["content"]), denom=1,
                coord0=val(*BL["coord0"]), coord1=val(*BL["coord1"]),
                coord2=val(*BL["coord2"]), coord3=val(*BL["coord3"]))


def load_known(t, r):
    txt = f"{r['content']:x}\n{r['denom']:x}\n{r['coord0']:x}\n{r['coord1']:x}\n{r['coord2']:x}\n{r['coord3']:x}\n"
    data = txt.encode()
    if ss(t, "x", b"", 1, 5000) != bytearray([1]): raise RuntimeError("reset")
    padded = data + b"\n" * ((-len(data)) % CHUNK)
    for off in range(0, len(padded), CHUNK):
        if ss(t, "a", padded[off:off + CHUNK], 1, 5000) != bytearray([1]): raise RuntimeError("append")
    if ss(t, "k", b"", 1, 5000) != bytearray([1]): raise RuntimeError("set-known")


class Pico:
    def __init__(s, samples, tb, pre, cha=5, chb=7):
        s.h = ctypes.c_int16(); s.samples = samples; s.tb = tb; s.pre = pre; s.post = samples - pre
        s.cha = cha; s.chb = chb
        st = ps.ps3000aOpenUnit(ctypes.byref(s.h), None)
        try: assert_pico_ok(st)
        except Exception:
            if st in (282, 286): assert_pico_ok(ps.ps3000aChangePowerSource(s.h, st))
            else: raise
        assert_pico_ok(ps.ps3000aSetChannel(s.h, 0, 1, 0, cha, 0))
        assert_pico_ok(ps.ps3000aSetChannel(s.h, 1, 1, 1, chb, 0))
        ps.ps3000aSetChannel(s.h, 2, 0, 1, 0, 0); ps.ps3000aSetChannel(s.h, 3, 0, 1, 0, 0)
        assert_pico_ok(ps.ps3000aSetSimpleTrigger(s.h, 1, 1, 1000, 2, 0, 5000))
        iv = ctypes.c_float(); mx = ctypes.c_int16()
        assert_pico_ok(ps.ps3000aGetTimebase2(s.h, tb, samples, ctypes.byref(iv), 1, ctypes.byref(mx), 0))
        s.interval_ns = iv.value; s.maxADC = ctypes.c_int16(); ps.ps3000aMaximumValue(s.h, ctypes.byref(s.maxADC))
        print(f"PicoScope: {s.interval_ns:.2f} ns/sample, {samples} samples "
              f"(~{samples*s.interval_ns/1e6:.2f} ms)", flush=True)

    def arm(s):
        s.bA = (ctypes.c_int16 * s.samples)(); s.bB = (ctypes.c_int16 * s.samples)()
        assert_pico_ok(ps.ps3000aSetDataBuffers(s.h, 0, s.bA, None, s.samples, 0, 0))
        assert_pico_ok(ps.ps3000aSetDataBuffers(s.h, 1, s.bB, None, s.samples, 0, 0))
        assert_pico_ok(ps.ps3000aRunBlock(s.h, s.pre, s.post, s.tb, 1, None, 0, None, None))

    def get(s):
        rdy = ctypes.c_int16(0)
        for _ in range(20_000_000):
            ps.ps3000aIsReady(s.h, ctypes.byref(rdy))
            if rdy.value: break
        ov = ctypes.c_int16(); got = ctypes.c_int32(s.samples)
        assert_pico_ok(ps.ps3000aGetValues(s.h, 0, ctypes.byref(got), 0, 0, 0, ctypes.byref(ov)))
        p = np.array(adc2mV(s.bA, s.cha, s.maxADC), dtype=np.float32)
        tr = np.array(adc2mV(s.bB, s.chb, s.maxADC), dtype=np.float32)
        return p, tr, ov.value

    def close(s): ps.ps3000aStop(s.h); ps.ps3000aCloseUnit(s.h)


def pulse_edges(trig, thr=1000):
    hi = trig > thr
    rise = np.where((~hi[:-1]) & hi[1:])[0] + 1
    fall = np.where(hi[:-1] & (~hi[1:]))[0] + 1
    return rise, fall


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-n", "--traces", type=int, default=5000)
    ap.add_argument("--samples", type=int, default=400000)
    ap.add_argument("--timebase", type=int, default=3)   # ~8 ns
    ap.add_argument("--pre", type=int, default=500)
    ap.add_argument("--skip-flash", action="store_true")
    ap.add_argument("--out", default="/tmp/sqisign_profile_known.npz")
    a = ap.parse_args()
    if not HEX.exists(): sys.exit("firmware missing (mode 6 + v)")

    scope = cw.scope(); target = cw.target(scope); pico = None
    T = None; records = []; edges = []; trig_last = None
    try:
        scope.default_setup()
        if not a.skip_flash:
            print("programming mode-6+v firmware ...", flush=True)
            cw.program_target(scope, STM32FProgrammer, str(HEX))
        time.sleep(0.2); target.flush()
        pico = Pico(a.samples, a.timebase, a.pre)
        T = np.lib.format.open_memmap(a.out.replace(".npz", "_traces.npy"), mode="w+",
                                      dtype=np.float32, shape=(a.traces, a.samples))
        rng = random.Random(0x5EED)
        t0 = time.time()
        for i in range(a.traces):
            r = gen_response(rng)
            load_known(target, r)
            pico.arm(); target.flush()
            target.simpleserial_write("s", bytearray())
            power, trig, ov = pico.get()
            if target.simpleserial_read("r", 1, timeout=8000) != bytearray([1]):
                raise RuntimeError(f"serialize ack {i}")
            T[i] = power; trig_last = trig
            rise, fall = pulse_edges(trig)
            e = np.full(12, -1, np.int64)
            e[:min(6, len(rise))] = rise[:6]; e[6:6 + min(6, len(fall))] = fall[:6]
            edges.append(e)
            records.append({k: hex(v) for k, v in r.items()})
            if (i + 1) % 500 == 0:
                rate = (i + 1) / (time.time() - t0)
                print(f"{i+1:5d}/{a.traces} pulses={len(rise)} rate={rate:.1f}/s "
                      f"eta={ (a.traces-i-1)/rate/60:.1f}min ov={ov&1}", flush=True)
    finally:
        if pico is not None: pico.close()
        try: target.dis()
        finally: scope.dis()
    if T is not None: T.flush()
    np.savez_compressed(a.out, records_json=json.dumps(records),
                        interval_ns=np.float32(pico.interval_ns if pico else 0),
                        edges=np.array(edges) if edges else np.zeros((0, 12)),
                        trig_last=trig_last if trig_last is not None else np.zeros(1),
                        n=a.traces, samples=a.samples)
    print(f"\nsaved {a.out}; traces -> {a.out.replace('.npz','_traces.npy')}", flush=True)


if __name__ == "__main__":
    main()
