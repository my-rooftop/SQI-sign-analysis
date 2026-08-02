#!/usr/bin/env python3
"""
Capture REAL response serialization with per-value trigger pulses (firmware mode 6).
Each trace: new RNG seed -> compute real response (~12 s) -> serialize under 6 trigger
pulses (content, denom, coord0, coord1, coord2, coord3) -> PicoScope captures the whole
window. Ground-truth coord0..3 are read back from the board.

Output npz: traces (N,S), records (hex strings), interval_ns.
Analyse separately (segment the 6 pulses on Ch B, per-coord SNR / recovery).

  ~/.cwvenv/bin/python3 exp_real_coords_capture.py --traces 300 --samples 760000
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
    print("+", " ".join(str(c) for c in cmd), flush=True)
    return subprocess.run(cmd, check=True, **kw)


def ss(t, cmd, data=b"", n=1, timeout_ms=250000):
    t.simpleserial_write(cmd, bytearray(data))
    return t.simpleserial_read("r", n, timeout=timeout_ms)


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
    text = out.decode("ascii", "replace")
    return text, parse_response(text)


def setup_vector(t):
    with VECTOR.open("wb") as f:
        run([str(HOST), "gen"], stdout=f)
    vec = VECTOR.read_bytes()
    if ss(t, "x", b"", 1, 5000) != bytearray([1]):
        raise RuntimeError("reset failed")
    padded = vec + b"\n" * ((-len(vec)) % CHUNK)
    for off in range(0, len(padded), CHUNK):
        if ss(t, "a", padded[off:off + CHUNK], 1, 5000) != bytearray([1]):
            raise RuntimeError(f"write failed at {off}")
    print(f"loaded vector bytes={len(vec)}", flush=True)


class Pico:
    def __init__(self, samples, timebase, pre, cha=5, chb=7):
        self.h = ctypes.c_int16(); self.samples = samples; self.timebase = timebase
        self.pre = pre; self.post = samples - pre; self.cha = cha; self.chb = chb
        st = ps.ps3000aOpenUnit(ctypes.byref(self.h), None)
        try:
            assert_pico_ok(st)
        except Exception:
            if st in (282, 286): assert_pico_ok(ps.ps3000aChangePowerSource(self.h, st))
            else: raise
        assert_pico_ok(ps.ps3000aSetChannel(self.h, 0, 1, 0, cha, 0))  # A power AC
        assert_pico_ok(ps.ps3000aSetChannel(self.h, 1, 1, 1, chb, 0))  # B trig DC
        ps.ps3000aSetChannel(self.h, 2, 0, 1, 0, 0); ps.ps3000aSetChannel(self.h, 3, 0, 1, 0, 0)
        assert_pico_ok(ps.ps3000aSetSimpleTrigger(self.h, 1, 1, 1000, 2, 0, 5000))
        iv = ctypes.c_float(); mx = ctypes.c_int16()
        assert_pico_ok(ps.ps3000aGetTimebase2(self.h, timebase, samples, ctypes.byref(iv), 1, ctypes.byref(mx), 0))
        self.interval_ns = iv.value
        self.maxADC = ctypes.c_int16(); ps.ps3000aMaximumValue(self.h, ctypes.byref(self.maxADC))
        print(f"PicoScope: {self.interval_ns:.2f} ns/sample, {samples} samples "
              f"(~{samples*self.interval_ns/1e6:.2f} ms)", flush=True)

    def arm(self):
        self.bA = (ctypes.c_int16 * self.samples)(); self.bB = (ctypes.c_int16 * self.samples)()
        assert_pico_ok(ps.ps3000aSetDataBuffers(self.h, 0, self.bA, None, self.samples, 0, 0))
        assert_pico_ok(ps.ps3000aSetDataBuffers(self.h, 1, self.bB, None, self.samples, 0, 0))
        assert_pico_ok(ps.ps3000aRunBlock(self.h, self.pre, self.post, self.timebase, 1, None, 0, None, None))

    def get(self):
        ready = ctypes.c_int16(0)
        for _ in range(20_000_000):
            ps.ps3000aIsReady(self.h, ctypes.byref(ready))
            if ready.value: break
        ov = ctypes.c_int16(); got = ctypes.c_int32(self.samples)
        assert_pico_ok(ps.ps3000aGetValues(self.h, 0, ctypes.byref(got), 0, 0, 0, ctypes.byref(ov)))
        p = np.array(adc2mV(self.bA, self.cha, self.maxADC), dtype=np.float32)
        tr = np.array(adc2mV(self.bB, self.chb, self.maxADC), dtype=np.float32)
        return p, tr, ov.value

    def close(self):
        ps.ps3000aStop(self.h); ps.ps3000aCloseUnit(self.h)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-n", "--traces", type=int, default=300)
    ap.add_argument("--samples", type=int, default=760000)
    ap.add_argument("--timebase", type=int, default=2)  # 4 ns
    ap.add_argument("--pre", type=int, default=500)
    ap.add_argument("--out", default="/tmp/sqisign_real_coords.npz")
    args = ap.parse_args()
    if not HEX.exists(): sys.exit("firmware missing (build mode 6)")
    if not HOST.exists(): run([str(APP / "build_host_tool.sh")])

    def pulse_edges(trig, thr=1000):
        hi = trig > thr
        rise = np.where((~hi[:-1]) & hi[1:])[0] + 1
        fall = np.where(hi[:-1] & (~hi[1:]))[0] + 1
        return rise, fall

    scope = cw.scope(); target = cw.target(scope); pico = None
    T = None; records = []; trig_last = None; edges = []
    try:
        scope.default_setup()
        print("programming mode-6 firmware ...", flush=True)
        cw.program_target(scope, STM32FProgrammer, str(HEX))
        time.sleep(0.2); target.flush()
        setup_vector(target)
        pico = Pico(args.samples, args.timebase, args.pre)
        T = np.lib.format.open_memmap(args.out.replace(".npz", "_traces.npy"), mode="w+",
                                      dtype=np.float32, shape=(args.traces, args.samples))
        t0 = time.time()
        for i in range(args.traces):
            seed_rng = random.Random(0xABCD0000 + i)
            seed = bytes(seed_rng.getrandbits(8) for _ in range(16))
            if ss(target, "n", seed, 1, 5000) != bytearray([1]):
                raise RuntimeError(f"seed failed {i}")
            r = ss(target, "c", b"", 1, 900000)          # compute real response (~12s)
            if r != bytearray([1]):
                raise RuntimeError(f"compute failed {i}: {r}")
            text, rec = read_board_response(target)      # ground truth coords
            pico.arm(); target.flush()
            target.simpleserial_write("s", bytearray())  # re-serialize w/ 6 pulses
            power, trig, ov = pico.get()
            if target.simpleserial_read("r", 1, timeout=8000) != bytearray([1]):
                raise RuntimeError(f"serialize ack failed {i}")
            T[i] = power; trig_last = trig
            rise, fall = pulse_edges(trig)
            npulse = len(rise)
            # store 6 rising + 6 falling edges (pad/truncate to 6); -1 if missing
            e = np.full(12, -1, dtype=np.int64)
            e[:min(6, len(rise))] = rise[:6]
            e[6:6 + min(6, len(fall))] = fall[:6]
            edges.append(e)
            rec["raw"] = text
            records.append({k: (hex(v) if isinstance(v, int) else v) for k, v in rec.items()})
            el = time.time() - t0
            print(f"{i+1:4d}/{args.traces} pulses={npulse} "
                  f"c0={abs(rec['coord0']).bit_length()}b c1={abs(rec['coord1']).bit_length()}b "
                  f"c2={abs(rec['coord2']).bit_length()}b c3={abs(rec['coord3']).bit_length()}b "
                  f"ov={ov&1} {el/ (i+1):.1f}s/tr", flush=True)
    finally:
        if pico is not None: pico.close()
        try: target.dis()
        finally: scope.dis()
    if T is not None: T.flush()
    np.savez_compressed(args.out, records_json=json.dumps(records),
                        interval_ns=np.float32(pico.interval_ns if pico else 0),
                        trig_last=trig_last if trig_last is not None else np.zeros(1),
                        edges=np.array(edges) if edges else np.zeros((0, 12)),
                        n=args.traces, samples=args.samples)
    print(f"\nsaved metadata {args.out}; traces -> {args.out.replace('.npz','_traces.npy')}", flush=True)


if __name__ == "__main__":
    main()
