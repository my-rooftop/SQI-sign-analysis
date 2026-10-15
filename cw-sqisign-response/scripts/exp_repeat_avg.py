#!/usr/bin/env python3
"""
Averaging-gain experiment: K distinct fixed-bit-length coord0 values, each captured
R times. Measures how DATA leakage (LS-limb Hamming weight of coord0) and its
recoverability grow as we average m repeats of the SAME value.

Answers: if the attacker could replay a value m times, how many repeats to recover
the data? (In real SQIsign each signature's coords are unique -> m=1 only.)

content=1, denom=1 so coord0 is serialized early -> small capture window.
"""
import argparse, ctypes, json, pathlib, random, sys, time
import numpy as np
import chipwhisperer as cw
from chipwhisperer.capture.api.programmers import STM32FProgrammer
from picosdk.ps3000a import ps3000a as ps
from picosdk.functions import adc2mV, assert_pico_ok

APP = pathlib.Path(__file__).resolve().parent
FW_DIR = APP.parent / "firmware"
HEX = FW_DIR / "simpleserial-sqisign-response-CW308_STM32F4.hex"
CHUNK = 48
BITLEN = 320  # fixed coord0 size


def ss(t, cmd, data=b"", n=1, tmo=8000):
    t.simpleserial_write(cmd, bytearray(data)); return t.simpleserial_read("r", n, timeout=tmo)


def load_known(t, coord0):
    txt = f"1\n1\n{coord0:x}\n2\n3\n4\n"
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
        print(f"PicoScope: {s.interval_ns:.2f} ns/sample, {samples} samples", flush=True)

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
    ap.add_argument("-K", "--values", type=int, default=64)
    ap.add_argument("-R", "--repeats", type=int, default=48)
    ap.add_argument("--samples", type=int, default=220000)
    ap.add_argument("--timebase", type=int, default=3)
    ap.add_argument("--pre", type=int, default=500)
    ap.add_argument("--out", default="/tmp/sqisign_repeat_avg.npz")
    a = ap.parse_args()
    if not HEX.exists(): sys.exit("firmware missing")

    rng = random.Random(0xA5A5)
    values = [(1 << (BITLEN - 1)) | rng.getrandbits(BITLEN - 1) for _ in range(a.values)]
    Ntot = a.values * a.repeats
    scope = cw.scope(); target = cw.target(scope); pico = None
    T = None; vidx = []; edges = []
    try:
        scope.default_setup()
        print("programming firmware ...", flush=True)
        cw.program_target(scope, STM32FProgrammer, str(HEX)); time.sleep(0.2); target.flush()
        pico = Pico(a.samples, a.timebase, a.pre)
        T = np.lib.format.open_memmap(a.out.replace(".npz", "_traces.npy"), mode="w+",
                                      dtype=np.float32, shape=(Ntot, a.samples))
        order = [(rep, k) for rep in range(a.repeats) for k in range(a.values)]  # interleaved
        t0 = time.time()
        for n, (rep, k) in enumerate(order):
            load_known(target, values[k])
            pico.arm(); target.flush()
            target.simpleserial_write("s", bytearray())
            power, trig, ov = pico.get()
            if target.simpleserial_read("r", 1, timeout=8000) != bytearray([1]):
                raise RuntimeError(f"ack {n}")
            T[n] = power; vidx.append(k)
            r, f = pulse_edges(trig)
            e = np.full(12, -1, np.int64); e[:min(6, len(r))] = r[:6]; e[6:6 + min(6, len(f))] = f[:6]
            edges.append(e)
            if (n + 1) % 500 == 0:
                rate = (n + 1) / (time.time() - t0)
                print(f"{n+1:5d}/{Ntot} rate={rate:.1f}/s eta={(Ntot-n-1)/rate/60:.1f}min", flush=True)
    finally:
        if pico is not None: pico.close()
        try: target.dis()
        finally: scope.dis()
    if T is not None: T.flush()
    np.savez_compressed(a.out, values=np.array(values, dtype=object).astype(str),
                        vidx=np.array(vidx), edges=np.array(edges),
                        interval_ns=np.float32(pico.interval_ns if pico else 0),
                        K=a.values, R=a.repeats, samples=a.samples, bitlen=BITLEN)
    print(f"\nsaved {a.out}", flush=True)


if __name__ == "__main__":
    main()
