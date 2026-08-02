#!/usr/bin/env python3
"""
Stockpile capture of the ATTACK-TARGET function on random inputs.
Firmware RESP_TRIGGER_MODE=2: trigger wraps sample_finalize (coord = basis*x, denom
copy, normalize) -- where coord0..3 are first materialized as integers.

Each trace: fresh RNG seed -> 'c' runs the full response core (~12 s); the trigger
fires DURING finalize (near the end of the compute). PicoScope must wait indefinitely
(autoTrigger=0) because the trigger arrives ~11-12 s after arming. Ground-truth
coord0..3 are read back.

  ~/.cwvenv/bin/python3 exp_target_finalize.py --traces 3000 --samples 1300000 --timebase 4
  ~/.cwvenv/bin/python3 exp_target_finalize.py --probe    # 2 traces, size the window
"""
import argparse, ctypes, json, pathlib, random, subprocess, sys, time
import numpy as np
import chipwhisperer as cw
from chipwhisperer.capture.api.programmers import STM32FProgrammer
from picosdk.ps3000a import ps3000a as ps
from picosdk.functions import adc2mV, assert_pico_ok

APP = pathlib.Path(__file__).resolve().parent
HOST = APP / "response_split_host"
HEX = APP / "simpleserial-sqisign-response-CW308_STM32F4.hex"
VECTOR = pathlib.Path("/tmp/sqisign_response_vec.txt")
CHUNK = 48


def run(cmd, **kw):
    print("+", " ".join(str(c) for c in cmd), flush=True); return subprocess.run(cmd, check=True, **kw)


def ss(t, cmd, data=b"", n=1, tmo=250000):
    t.simpleserial_write(cmd, bytearray(data)); return t.simpleserial_read("r", n, timeout=tmo)


def parse_response(text):
    toks = text.split()
    if len(toks) < 7 or toks[0] != "1":
        raise ValueError(f"bad response: {text[:80]!r}")
    v = [int(x, 16) for x in toks[1:7]]
    return dict(lattice_content=v[0], denom=v[1], coord0=v[2], coord1=v[3], coord2=v[4], coord3=v[5])


def read_board_response(t):
    r = ss(t, "l", b"", 2, 5000); outlen = (r[0] << 8) | r[1]
    out = bytearray()
    for idx in range((outlen + CHUNK - 1) // CHUNK):
        out.extend(ss(t, "g", bytes([(idx >> 8) & 0xFF, idx & 0xFF]), CHUNK, 5000))
    del out[outlen:]
    return parse_response(out.decode("ascii", "replace"))


def setup_vector(t):
    with VECTOR.open("wb") as f:
        run([str(HOST), "gen"], stdout=f)
    vec = VECTOR.read_bytes()
    if ss(t, "x", b"", 1, 5000) != bytearray([1]): raise RuntimeError("reset")
    padded = vec + b"\n" * ((-len(vec)) % CHUNK)
    for off in range(0, len(padded), CHUNK):
        if ss(t, "a", padded[off:off + CHUNK], 1, 5000) != bytearray([1]): raise RuntimeError("append")
    print(f"loaded vector bytes={len(vec)}", flush=True)


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
        # autoTrigger = 0  -> wait indefinitely (finalize trigger arrives ~12 s after arming)
        assert_pico_ok(ps.ps3000aSetSimpleTrigger(s.h, 1, 1, 1000, 2, 0, 0))
        iv = ctypes.c_float(); mx = ctypes.c_int16()
        assert_pico_ok(ps.ps3000aGetTimebase2(s.h, tb, samples, ctypes.byref(iv), 1, ctypes.byref(mx), 0))
        s.interval_ns = iv.value; s.maxADC = ctypes.c_int16(); ps.ps3000aMaximumValue(s.h, ctypes.byref(s.maxADC))
        print(f"PicoScope: {s.interval_ns:.1f} ns/sample, {samples} samples "
              f"(~{samples*s.interval_ns/1e6:.2f} ms window)", flush=True)

    def arm(s):
        s.bA = (ctypes.c_int16 * s.samples)(); s.bB = (ctypes.c_int16 * s.samples)()
        assert_pico_ok(ps.ps3000aSetDataBuffers(s.h, 0, s.bA, None, s.samples, 0, 0))
        assert_pico_ok(ps.ps3000aSetDataBuffers(s.h, 1, s.bB, None, s.samples, 0, 0))
        assert_pico_ok(ps.ps3000aRunBlock(s.h, s.pre, s.post, s.tb, 1, None, 0, None, None))

    def get(s, poll_s=0.02):
        rdy = ctypes.c_int16(0)
        while rdy.value == 0:                       # waits through the ~12 s compute
            ps.ps3000aIsReady(s.h, ctypes.byref(rdy))
            time.sleep(poll_s)
        ov = ctypes.c_int16(); got = ctypes.c_int32(s.samples)
        assert_pico_ok(ps.ps3000aGetValues(s.h, 0, ctypes.byref(got), 0, 0, 0, ctypes.byref(ov)))
        p = np.array(adc2mV(s.bA, s.cha, s.maxADC), dtype=np.float32)
        tr = np.array(adc2mV(s.bB, s.chb, s.maxADC), dtype=np.float32)
        return p, tr, ov.value

    def close(s): ps.ps3000aStop(s.h); ps.ps3000aCloseUnit(s.h)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-n", "--traces", type=int, default=3000)
    ap.add_argument("--samples", type=int, default=1300000)
    ap.add_argument("--timebase", type=int, default=4)   # ~16 ns
    ap.add_argument("--pre", type=int, default=2000)
    ap.add_argument("--probe", action="store_true")
    ap.add_argument("--out", default="/tmp/sqisign_target_finalize.npz")
    a = ap.parse_args()
    if a.probe: a.traces = 2
    if not HEX.exists(): sys.exit("firmware missing (build mode 2)")
    if not HOST.exists(): run([str(APP / "build_host_tool.sh")])

    def pulse_edges(trig, thr=1000):
        hi = trig > thr
        rise = np.where((~hi[:-1]) & hi[1:])[0] + 1
        fall = np.where(hi[:-1] & (~hi[1:]))[0] + 1
        return rise, fall

    scope = cw.scope(); target = cw.target(scope); pico = None
    T = None; records = []; trig_last = None; widths = []; edges = []
    try:
        scope.default_setup()
        print("programming mode-7 firmware ...", flush=True)
        cw.program_target(scope, STM32FProgrammer, str(HEX)); time.sleep(0.2); target.flush()
        setup_vector(target)
        pico = Pico(a.samples, a.timebase, a.pre)
        if not a.probe:
            T = np.lib.format.open_memmap(a.out.replace(".npz", "_traces.npy"), mode="w+",
                                          dtype=np.float32, shape=(a.traces, a.samples))
        t0 = time.time()
        for i in range(a.traces):
            seed = bytes(random.Random(0xC0DE0000 + i).getrandbits(8) for _ in range(16))
            if ss(target, "n", seed, 1, 5000) != bytearray([1]): raise RuntimeError(f"seed {i}")
            pico.arm(); target.flush()
            target.simpleserial_write("c", bytearray())   # 12 s core; trigger fires in finalize
            power, trig, ov = pico.get()
            if target.simpleserial_read("r", 1, timeout=900000) != bytearray([1]):
                raise RuntimeError(f"compute ack {i}")
            rec = read_board_response(target)
            rise, fall = pulse_edges(trig)
            e = np.full(16, -1, np.int64)                 # up to 8 rise + 8 fall (mode 7)
            e[:min(8, len(rise))] = rise[:8]; e[8:8 + min(8, len(fall))] = fall[:8]
            edges.append(e)
            w = int(fall[-1] - rise[0]) if (len(rise) and len(fall)) else 0
            widths.append(w); trig_last = trig
            if T is not None: T[i] = power
            records.append({k: hex(v) for k, v in rec.items()})
            el = time.time() - t0
            print(f"{i+1:4d}/{a.traces} pulses={len(rise)} span={w}({w*pico.interval_ns/1e3:.0f}us) "
                  f"c0={abs(rec['coord0']).bit_length()}b ov={ov&1} {el/(i+1):.1f}s/tr "
                  f"eta={(a.traces-i-1)*el/(i+1)/3600:.1f}h", flush=True)
    finally:
        if pico is not None: pico.close()
        try: target.dis()
        finally: scope.dis()
    if T is not None:
        T.flush()
        np.savez_compressed(a.out, records_json=json.dumps(records),
                            interval_ns=np.float32(pico.interval_ns if pico else 0),
                            trig_last=trig_last if trig_last is not None else np.zeros(1),
                            edges=np.array(edges) if edges else np.zeros((0, 16)),
                            widths=np.array(widths), n=a.traces, samples=a.samples)
        print(f"\nsaved {a.out}; traces -> {a.out.replace('.npz','_traces.npy')}", flush=True)
    else:
        for e in edges:
            r = [int(v) for v in e[:8] if v > 0]
            print(f"PROBE pulses={len(r)} rise_samples={r}", flush=True)
        print(f"PROBE spans(us)={[int(w*pico.interval_ns/1e3) for w in widths]}", flush=True)


if __name__ == "__main__":
    main()
