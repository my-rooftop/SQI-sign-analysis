#!/usr/bin/env python3
"""
SQIsign response leakage capture -- PicoScope 3000E (psospa) + CW-Lite + STM32F415
==================================================================================
Same rig as capture_response_leakage_pico.py, but for PicoScope 3000E-series units
(3405E/3415E/...). Those enumerate under the `psospa` driver, NOT `ps3000a`, so the
older script fails on them with PICO_NOT_FOUND (status 3).

  ChipWhisperer-Lite : programs the target, supplies the 7.3728 MHz target clock,
                       drives SimpleSerial.
  PicoScope 3000E    : captures the waveform.
                         Ch A = power (CW308 SMA), AC-coupled
                         Ch B = trigger line, DC-coupled, hardware trigger source

Driver differences vs ps3000a that this file encodes:
  - psospaOpenUnit takes a resolution enum; ranges are given in NANOVOLTS as an
    explicit (min,max) pair rather than a range-code index.
  - psospaGetTimebase returns the sample interval in NANOSECONDS (while
    psospaGetMinimumTimebaseStateless returns SECONDS -- easy to trip over).
  - Channel ranges travel through a c_long in nV, so 2 V is the largest range this
    binding can express. A 3.3 V trigger line therefore clips on Ch B; harmless,
    since we only threshold it for edges, but it does raise B's overflow bit.
  - Measured on a 3415E: 8- and 10-bit open fine, 12-bit is rejected outright
    (status 379). 10-bit is free at 8 ns/sample -- it only costs sample rate at the
    very fastest timebases (2.4 -> 3.2 ns at timebase 2000) -- so it is the default
    here: 0.39 mV/LSB at +/-200 mV versus 3.94 mV at 8-bit / +/-500 mV.

Default build is RESP_TRIGGER_MODE=6 (per-value pulses): one 's' gives six aligned
segments -- lattice_content, denom, coord0, coord1, coord2, coord3.

  ~/.cwvenv/bin/python3 scripts/capture_response_leakage_psospa.py --traces 10
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
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import chipwhisperer as cw
from chipwhisperer.capture.api.programmers import STM32FProgrammer
from picosdk.psospa import psospa as ps
from picosdk.PicoDeviceEnums import picoEnum as enums

APP = pathlib.Path(__file__).resolve().parent
PKG = APP.parent
FW_DIR = PKG / "firmware"
HOST = PKG / "host" / "response_split_host"
HOST_BUILD = PKG / "host" / "build_host_tool.sh"
HEX = FW_DIR / "simpleserial-sqisign-response-CW308_STM32F4.hex"
VECTOR = pathlib.Path("/tmp/sqisign_response_vec.txt")
CHUNK = 48

# ranges the binding can express (nV must fit in c_long -> 2 V ceiling)
RANGES_MV = [10, 20, 50, 100, 200, 500, 1000, 2000]
AUTORANGE_HEADROOM = 1.4    # picked range must exceed peak by this factor
WARMUP_CAPTURES = 5         # AC-coupling settle before any measurement is trusted
PULSE_NAMES = ["lattice_content", "denom", "coord0", "coord1", "coord2", "coord3"]


def run(cmd, **kw):
    print("+", " ".join(str(c) for c in cmd), flush=True)
    return subprocess.run(cmd, check=True, **kw)


def ss(target, cmd, data=b"", n=1, timeout_ms=8000):
    target.simpleserial_write(cmd, bytearray(data))
    return target.simpleserial_read("r", n, timeout=timeout_ms)


def parse_response(text):
    toks = text.split()
    if len(toks) < 7 or toks[0] != "1":
        raise ValueError(f"bad response: {text[:80]!r}")
    v = [int(x, 16) for x in toks[1:7]]
    return dict(lattice_content=v[0], denom=v[1],
                coord0=v[2], coord1=v[3], coord2=v[4], coord3=v[5])


def read_board_response(target):
    r = ss(target, "l", b"", 2, 5000)
    outlen = (r[0] << 8) | r[1]
    out = bytearray()
    for idx in range((outlen + CHUNK - 1) // CHUNK):
        out.extend(ss(target, "g", bytes([(idx >> 8) & 0xFF, idx & 0xFF]), CHUNK, 5000))
    del out[outlen:]
    text = out.decode("ascii", "replace")
    return text, parse_response(text)


# realistic lvl1 coordinate bit-length ranges (RESPONSE_FLOW_AND_MEASUREMENTS.md)
BL = {"content": (516, 528), "coord0": (316, 328), "coord1": (316, 327),
      "coord2": (188, 202), "coord3": (188, 201)}


def gen_response(rng):
    def val(lo, hi):
        b = rng.randint(lo, hi)
        return (1 << (b - 1)) | rng.getrandbits(b - 1)
    return dict(lattice_content=val(*BL["content"]), denom=1,
                coord0=val(*BL["coord0"]), coord1=val(*BL["coord1"]),
                coord2=val(*BL["coord2"]), coord3=val(*BL["coord3"]))


def load_known(target, r):
    """cmd 'k': inject a known response into last_resp, bypassing the ~12 s core.
    Needed on targets whose RAM cannot host the core at all (e.g. STM32F3)."""
    txt = (f"{r['lattice_content']:x}\n{r['denom']:x}\n{r['coord0']:x}\n"
           f"{r['coord1']:x}\n{r['coord2']:x}\n{r['coord3']:x}\n").encode()
    if ss(target, "x", b"", 1, 5000) != bytearray([1]):
        raise RuntimeError("reset failed")
    padded = txt + b"\n" * ((-len(txt)) % CHUNK)
    for off in range(0, len(padded), CHUNK):
        if ss(target, "a", padded[off:off + CHUNK], 1, 5000) != bytearray([1]):
            raise RuntimeError("known-response append failed")
    if ss(target, "k", b"", 1, 5000) != bytearray([1]):
        raise RuntimeError("set-known failed")


def setup_vector(target):
    with VECTOR.open("wb") as f:
        run([str(HOST), "gen"], stdout=f)
    vec = VECTOR.read_bytes()
    if ss(target, "x", b"", 1, 5000) != bytearray([1]):
        raise RuntimeError("reset failed")
    padded = vec + b"\n" * ((-len(vec)) % CHUNK)
    for off in range(0, len(padded), CHUNK):
        if ss(target, "a", padded[off:off + CHUNK], 1, 5000) != bytearray([1]):
            raise RuntimeError(f"vector append failed at {off}")
    print(f"loaded vector bytes={len(vec)} chunks={len(padded) // CHUNK}", flush=True)


class Pico3000E:
    """Ch A = power (AC), Ch B = trigger (DC). Ranges in mV, converted to nV."""

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
        self.interval_ns = interval.value   # psospaGetTimebase reports ns, not s
        print(f"PicoScope 3000E: {self.variant} res={self.res_bits}bit "
              f"{self.interval_ns:.2f} ns/sample, {samples} samples "
              f"(~{samples * self.interval_ns / 1e6:.3f} ms window), pre={pre}", flush=True)

    # -- setup helpers -------------------------------------------------
    def _ok(self, status, what=""):
        if status != 0:
            raise RuntimeError(f"psospa error {status} {what}")
        return status

    def _open(self, bits):
        for want in (bits, 8):
            res = enums.PICO_DEVICE_RESOLUTION[f"PICO_DR_{want}BIT"]
            st = ps.psospaOpenUnit(ctypes.byref(self.h), None, res, None)
            if st == 0 and self.h.value > 0:
                self.res_bits = want
                self.variant = self._info(3)
                return res
            if self.h.value > 0:
                ps.psospaCloseUnit(self.h)
            print(f"  {want}-bit open failed (status {st}), falling back", flush=True)
        raise RuntimeError("could not open any psospa device")

    def _info(self, line):
        buf = ctypes.create_string_buffer(128)
        rl = ctypes.c_int16()
        ps.psospaGetUnitInfo(self.h, buf, 128, ctypes.byref(rl), line)
        return buf.value.decode(errors="replace")

    def set_ranges(self, cha_mv, chb_mv):
        self.cha_mv, self.chb_mv = cha_mv, chb_mv
        full = enums.PICO_BANDWIDTH_LIMITER["PICO_BW_FULL"]
        self._ok(ps.psospaSetChannelOn(self.h, self.A, enums.PICO_COUPLING["PICO_AC"],
                                       -cha_mv * 1000000, cha_mv * 1000000, 0, 0, full), "chA")
        self._ok(ps.psospaSetChannelOn(self.h, self.B, enums.PICO_COUPLING["PICO_DC"],
                                       -chb_mv * 1000000, chb_mv * 1000000, 0, 0, full), "chB")

    def set_trigger(self, threshold_mV=1000, auto_us=5_000_000):
        """auto_us=0 waits forever (needed when the trigger lands ~12 s after arming)."""
        thr = round(threshold_mV * self.maxADC.value / self.chb_mv)
        self._ok(ps.psospaSetSimpleTrigger(self.h, 1, self.B, thr,
                                           enums.PICO_THRESHOLD_DIRECTION["PICO_RISING"],
                                           0, auto_us), "trigger")

    # -- capture -------------------------------------------------------
    def arm(self):
        self.bufA = (ctypes.c_int16 * self.samples)()
        self.bufB = (ctypes.c_int16 * self.samples)()
        dtype = enums.PICO_DATA_TYPE["PICO_INT16_T"]
        raw = enums.PICO_RATIO_MODE["PICO_RATIO_MODE_RAW"]
        action = enums.PICO_ACTION["PICO_CLEAR_ALL"] | enums.PICO_ACTION["PICO_ADD"]
        self._ok(ps.psospaSetDataBuffer(self.h, self.A, ctypes.byref(self.bufA),
                                        self.samples, dtype, 0, raw, action), "bufA")
        self._ok(ps.psospaSetDataBuffer(self.h, self.B, ctypes.byref(self.bufB),
                                        self.samples, dtype, 0, raw,
                                        enums.PICO_ACTION["PICO_ADD"]), "bufB")
        t_ind = ctypes.c_double(0)
        self._ok(ps.psospaRunBlock(self.h, self.pre, self.post, self.timebase,
                                   ctypes.byref(t_ind), 0, None, None), "runBlock")

    def get(self, timeout_s=30.0):
        ready = ctypes.c_int16(0)
        t0 = time.time()
        while ready.value == 0:
            ps.psospaIsReady(self.h, ctypes.byref(ready))
            if time.time() - t0 > timeout_s:
                raise TimeoutError("capture did not complete (no trigger?)")
        n = ctypes.c_uint64(self.samples)
        overflow = ctypes.c_int16(0)
        raw = enums.PICO_RATIO_MODE["PICO_RATIO_MODE_RAW"]
        self._ok(ps.psospaGetValues(self.h, 0, ctypes.byref(n), 1, raw, 0,
                                    ctypes.byref(overflow)), "getValues")
        # vectorised adc -> mV (picosdk's adc2mVV2 is a Python loop; far too slow here)
        a = np.ctypeslib.as_array(self.bufA).astype(np.float32) * (self.cha_mv / self.maxADC.value)
        b = np.ctypeslib.as_array(self.bufB).astype(np.float32) * (self.chb_mv / self.maxADC.value)
        return a.copy(), b.copy(), overflow.value

    def close(self):
        ps.psospaStop(self.h)
        ps.psospaCloseUnit(self.h)


def pulse_edges(trig, thr_mv=1000.0):
    hi = trig > thr_mv
    rise = np.where((~hi[:-1]) & hi[1:])[0] + 1
    fall = np.where(hi[:-1] & (~hi[1:]))[0] + 1
    return rise, fall


def main():
    p = argparse.ArgumentParser(description="SQIsign response capture on PicoScope 3000E")
    p.add_argument("-n", "--traces", type=int, default=10)
    p.add_argument("--samples", type=int, default=400000)
    p.add_argument("--timebase", type=int, default=8000, help="8000 -> ~8 ns/sample")
    p.add_argument("--pre", type=int, default=2000)
    p.add_argument("--cha-mv", type=int, default=500, help="Ch A range, +/- mV")
    p.add_argument("--chb-mv", type=int, default=2000, help="Ch B range, +/- mV")
    p.add_argument("--resolution", type=int, default=10, choices=[8, 10, 12],
                   help="ADC bits. The 3415E supports 8 and 10 only -- it rejects 12 "
                        "at open time (status 379); 10 costs nothing at 8 ns/sample")
    p.add_argument("--no-autorange", action="store_true")
    p.add_argument("--known", action="store_true",
                   help="inject known responses with cmd 'k' instead of running the "
                        "~12 s core; required on targets too small for the core")
    p.add_argument("--skip-flash", action="store_true")
    p.add_argument("--platform", default="CW308_STM32F4",
                   help="target platform; picks the matching .hex in firmware/")
    p.add_argument("--firmware", type=pathlib.Path, default=None,
                   help="explicit .hex path (overrides --platform)")
    p.add_argument("--out", default="/tmp/sqisign_psospa_capture.npz")
    p.add_argument("--plot", default="/tmp/sqisign_psospa_capture.png")
    a = p.parse_args()

    global HEX
    HEX = a.firmware or (FW_DIR / f"simpleserial-sqisign-response-{a.platform}.hex")
    if not HEX.exists():
        sys.exit(f"firmware missing: {HEX}\n"
                 f"build: make -C {FW_DIR} PLATFORM={a.platform} RESP_TRIGGER_MODE=6")
    if not HOST.exists():
        run([str(HOST_BUILD)])

    scope = cw.scope()
    target = cw.target(scope)
    pico = None
    traces, trigs, records, edges = [], [], [], []
    try:
        scope.default_setup()
        if not a.skip_flash:
            print(f"programming {HEX.name} ...", flush=True)
            cw.program_target(scope, STM32FProgrammer, str(HEX))
        # the target must be reset and given time to boot, otherwise the first
        # SimpleSerial command lands mid-banner and comes back as a framing error
        scope.io.nrst = "low"
        time.sleep(0.05)
        scope.io.nrst = "high_z"
        time.sleep(0.5)
        target.flush()
        print("banner:", target.read(timeout=1000), flush=True)

        rng = random.Random(0x5EED)
        if a.known:
            print("known-value mode (cmd 'k'): core is never run", flush=True)
            load_known(target, gen_response(rng))
        else:
            setup_vector(target)
            print("computing one real response (cmd 'c', ~12 s) ...", flush=True)
            t0 = time.time()
            if ss(target, "c", b"", 1, 900000) != bytearray([1]):
                raise RuntimeError("compute failed")
            _txt, rec0 = read_board_response(target)
            print(f"compute OK in {time.time() - t0:.1f}s  "
                  f"coord0={abs(rec0['coord0']).bit_length()}b "
                  f"coord1={abs(rec0['coord1']).bit_length()}b "
                  f"coord2={abs(rec0['coord2']).bit_length()}b "
                  f"coord3={abs(rec0['coord3']).bit_length()}b", flush=True)

        pico = Pico3000E(a.samples, a.timebase, a.pre, a.cha_mv, a.chb_mv, a.resolution)
        pico.set_trigger(1000, auto_us=5_000_000)

        def one_capture():
            pico.arm()
            target.flush()
            target.simpleserial_write("s", bytearray())
            power, trig, ov = pico.get()
            ack = target.simpleserial_read("r", 1, timeout=8000)
            if ack != bytearray([1]):
                raise RuntimeError(f"serialize ack failed: {ack}")
            return power, trig, ov

        # AC coupling settles slowly; two captures was not enough -- the residual
        # transient inflated the measured peak (165 mV vs a true 87 mV) and pushed
        # autorange a whole range step too high.
        for _ in range(WARMUP_CAPTURES):
            one_capture()
        print(f"warm-up: {WARMUP_CAPTURES} captures discarded (AC settle)", flush=True)

        if not a.no_autorange:
            # a percentile, not max(): one stray settling sample must not decide the range
            peaks = []
            for _ in range(2):
                power, _t, _o = one_capture()
                peaks.append(float(np.percentile(np.abs(power), 99.99)))
            peak = max(peaks)
            best = next((r for r in RANGES_MV if peak * AUTORANGE_HEADROOM < r), a.cha_mv)
            if best != a.cha_mv:
                print(f"autorange: peak(p99.99) {peak:.1f} mV -> Ch A +/-{best} mV "
                      f"(was +/-{a.cha_mv})", flush=True)
                pico.set_ranges(best, a.chb_mv)
                pico.set_trigger(1000, auto_us=5_000_000)
                for _ in range(2):
                    one_capture()  # settle after the range change
            else:
                print(f"autorange: keeping Ch A +/-{a.cha_mv} mV "
                      f"(peak(p99.99) {peak:.1f} mV)", flush=True)

        for i in range(a.traces):
            if a.known:
                load_known(target, gen_response(rng))
            power, trig, ov = one_capture()
            text, rec = read_board_response(target)
            rise, fall = pulse_edges(trig)
            traces.append(power)
            trigs.append(trig)
            records.append({k: hex(v) for k, v in rec.items()})
            e = np.full(12, -1, np.int64)
            e[:min(6, len(rise))] = rise[:6]
            e[6:6 + min(6, len(fall))] = fall[:6]
            edges.append(e)
            widths_us = [(int(f) - int(r)) * pico.interval_ns / 1e3
                         for r, f in zip(rise[:6], fall[:6])]
            print(f"{i:03d}: pmV[{power.min():7.1f},{power.max():7.1f}] pulses={len(rise)} "
                  f"widths_us={[f'{w:.0f}' for w in widths_us]} "
                  f"ovA={ov & 1}", flush=True)
    finally:
        if pico is not None:
            pico.close()
        try:
            target.dis()
        finally:
            scope.dis()

    T = np.vstack(traces)
    E = np.array(edges)
    np.savez_compressed(a.out, traces=T, trig_last=trigs[-1], edges=E,
                        records_json=json.dumps(records),
                        interval_ns=np.float32(pico.interval_ns),
                        cha_mv=pico.cha_mv, resolution=pico.res_bits)
    print(f"\nsaved {T.shape} to {a.out}", flush=True)

    # ---- overview plot ----
    trig = trigs[-1]
    rise, fall = pulse_edges(trig)
    fig, ax = plt.subplots(3, 1, figsize=(15, 10))

    ax[0].plot(trig, color="green", lw=0.5)
    for k, (r, f) in enumerate(zip(rise[:6], fall[:6])):
        ax[0].axvspan(r, f, color="orange", alpha=0.20)
        ax[0].text((r + f) / 2, trig.max() * 0.55, PULSE_NAMES[k],
                   ha="center", fontsize=8, rotation=90)
    ax[0].set_title(f"Trigger (Ch B) -- {len(rise)} pulses, last trace")
    ax[0].set_ylabel("mV"); ax[0].grid(alpha=0.3)

    ax[1].plot(traces[-1], lw=0.25, color="C0")
    for r, f in zip(rise[:6], fall[:6]):
        ax[1].axvspan(r, f, color="orange", alpha=0.15)
    ax[1].set_title(f"Power (Ch A) -- last trace, +/-{pico.cha_mv} mV AC, "
                    f"{pico.res_bits}-bit, {pico.interval_ns:.1f} ns/sample")
    ax[1].set_ylabel("mV"); ax[1].grid(alpha=0.3)

    ax[2].plot(T.mean(axis=0), lw=0.4, color="black")
    for r, f in zip(rise[:6], fall[:6]):
        ax[2].axvspan(r, f, color="orange", alpha=0.15)
    ax[2].set_title(f"Mean power over {T.shape[0]} traces")
    ax[2].set_xlabel("sample"); ax[2].set_ylabel("mV"); ax[2].grid(alpha=0.3)

    fig.tight_layout()
    fig.savefig(a.plot, dpi=130)
    plt.close(fig)
    print(f"saved plot to {a.plot}", flush=True)


if __name__ == "__main__":
    main()
