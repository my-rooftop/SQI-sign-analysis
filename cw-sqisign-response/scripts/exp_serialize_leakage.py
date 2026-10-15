#!/usr/bin/env python3
"""
Experiment: does the response serialization leak the coordinate value?
======================================================================
Positive-control leakage test on the PicoScope + CW-Lite + F415 rig.

Firmware command 'q' loads a KNOWN 32-bit value into last_resp.coord[0]
(coord1..3 = 0, denom = 1, content = 1). Command 's' serializes last_resp
(append_ibz -> ASCII hex) under the trigger. So by varying the known value
across traces we can measure, with fully controlled ground truth, whether the
big-int -> hex serialization -- the same code path that emits all four real
coordinates -- produces value-dependent power leakage.

Metrics:
  - SNR keyed on Hamming weight of the 32-bit value (robust, few classes)
  - SNR keyed on each individual byte value (0..255)  -> per-byte recoverability
  - SNR keyed on HW of each byte
  - Pearson correlation with HW (quick sign check)

A clear SNR peak inside the trigger window => the setup captures exploitable
value leakage and coordinate recovery is plausible. Flat SNR => not visible at
this trace count / this region.

Usage:
  ~/.cwvenv/bin/python3 scripts/exp_serialize_leakage.py --traces 3000 --samples 60000
  ~/.cwvenv/bin/python3 scripts/exp_serialize_leakage.py --probe        # 20 traces, size the window
"""

import argparse
import ctypes
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
from picosdk.ps3000a import ps3000a as ps
from picosdk.functions import adc2mV, assert_pico_ok

APP = pathlib.Path(__file__).resolve().parent
FW_DIR = APP.parent / "firmware"
HEX = FW_DIR / "simpleserial-sqisign-response-CW308_STM32F4.hex"


def run(cmd, **kw):
    print("+", " ".join(str(c) for c in cmd), flush=True)
    return subprocess.run(cmd, check=True, **kw)


def ss(target, cmd, data=b"", n=1, timeout_ms=8000):
    target.simpleserial_write(cmd, bytearray(data))
    return target.simpleserial_read("r", n, timeout=timeout_ms)


class Pico:
    def __init__(self, samples, timebase, pre, cha_range=5, chb_range=7):
        self.h = ctypes.c_int16()
        self.samples, self.timebase = samples, timebase
        self.pre, self.post = pre, samples - pre
        self.cha_range, self.chb_range = cha_range, chb_range
        st = ps.ps3000aOpenUnit(ctypes.byref(self.h), None)
        try:
            assert_pico_ok(st)
        except Exception:
            if st in (282, 286):
                assert_pico_ok(ps.ps3000aChangePowerSource(self.h, st))
            else:
                raise
        assert_pico_ok(ps.ps3000aSetChannel(self.h, 0, 1, 0, cha_range, 0))  # A power, AC
        assert_pico_ok(ps.ps3000aSetChannel(self.h, 1, 1, 1, chb_range, 0))  # B trig, DC
        ps.ps3000aSetChannel(self.h, 2, 0, 1, 0, 0)
        ps.ps3000aSetChannel(self.h, 3, 0, 1, 0, 0)
        assert_pico_ok(ps.ps3000aSetSimpleTrigger(self.h, 1, 1, 1000, 2, 0, 5000))
        interval = ctypes.c_float(); mx = ctypes.c_int16()
        assert_pico_ok(ps.ps3000aGetTimebase2(self.h, timebase, samples,
                                              ctypes.byref(interval), 1, ctypes.byref(mx), 0))
        self.interval_ns = interval.value
        self.maxADC = ctypes.c_int16(); ps.ps3000aMaximumValue(self.h, ctypes.byref(self.maxADC))
        print(f"PicoScope: {self.interval_ns:.3f} ns/sample, {samples} samples "
              f"(~{samples*self.interval_ns/1e3:.1f} us)", flush=True)

    def arm(self):
        self.bA = (ctypes.c_int16 * self.samples)()
        self.bB = (ctypes.c_int16 * self.samples)()
        assert_pico_ok(ps.ps3000aSetDataBuffers(self.h, 0, self.bA, None, self.samples, 0, 0))
        assert_pico_ok(ps.ps3000aSetDataBuffers(self.h, 1, self.bB, None, self.samples, 0, 0))
        assert_pico_ok(ps.ps3000aRunBlock(self.h, self.pre, self.post, self.timebase, 1, None, 0, None, None))

    def get(self):
        ready = ctypes.c_int16(0)
        for _ in range(5_000_000):
            ps.ps3000aIsReady(self.h, ctypes.byref(ready))
            if ready.value:
                break
        ov = ctypes.c_int16(); got = ctypes.c_int32(self.samples)
        assert_pico_ok(ps.ps3000aGetValues(self.h, 0, ctypes.byref(got), 0, 0, 0, ctypes.byref(ov)))
        power = np.array(adc2mV(self.bA, self.cha_range, self.maxADC), dtype=np.float32)
        trig = np.array(adc2mV(self.bB, self.chb_range, self.maxADC), dtype=np.float32)
        return power, trig, ov.value

    def close(self):
        ps.ps3000aStop(self.h); ps.ps3000aCloseUnit(self.h)


def snr(traces, labels):
    """partition SNR = var(class means) / mean(class vars), weighted by class count."""
    labels = np.asarray(labels)
    classes = np.unique(labels)
    n, S = traces.shape
    means = np.empty((len(classes), S), dtype=np.float64)
    vars = np.empty((len(classes), S), dtype=np.float64)
    w = np.empty(len(classes), dtype=np.float64)
    for i, c in enumerate(classes):
        m = labels == c
        w[i] = m.sum()
        sub = traces[m]
        means[i] = sub.mean(0)
        vars[i] = sub.var(0)
    w /= w.sum()
    grand = (w[:, None] * means).sum(0)
    signal = (w[:, None] * (means - grand) ** 2).sum(0)   # between-class var
    noise = (w[:, None] * vars).sum(0)                    # within-class var
    return signal / (noise + 1e-12)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("-n", "--traces", type=int, default=3000)
    p.add_argument("--samples", type=int, default=60000)
    p.add_argument("--timebase", type=int, default=1, help="1 -> ~2 ns on this unit")
    p.add_argument("--pre", type=int, default=200)
    p.add_argument("--skip-flash", action="store_true")
    p.add_argument("--probe", action="store_true", help="20 traces, just size the trigger window")
    p.add_argument("--out", default="/tmp/sqisign_serialize_leakage.npz")
    p.add_argument("--plot", default="/tmp/sqisign_serialize_leakage.png")
    args = p.parse_args()
    if args.probe:
        args.traces = 20

    if not HEX.exists():
        sys.exit(f"firmware missing: {HEX} (build mode 1)")

    scope = cw.scope(); target = cw.target(scope); pico = None
    traces = []; values = []
    rng = random.Random(0xC0FFEE)
    try:
        scope.default_setup()
        if not args.skip_flash:
            print("programming mode-1 firmware ...", flush=True)
            cw.program_target(scope, STM32FProgrammer, str(HEX))
        time.sleep(0.2); target.flush()
        if ss(target, "x", b"", 1, 5000) != bytearray([1]):
            raise RuntimeError("reset failed")

        pico = Pico(args.samples, args.timebase, args.pre)

        # 'q' loads a known coord0; 's' serializes under trigger and returns an 'r' ack
        def capture(value):
            payload = bytes([(value >> 24) & 0xFF, (value >> 16) & 0xFF,
                             (value >> 8) & 0xFF, value & 0xFF])
            if ss(target, "q", payload, 1, 5000) != bytearray([1]):
                raise RuntimeError("calibration set failed")
            pico.arm(); target.flush()
            target.simpleserial_write("s", bytearray())
            power, trig, ov = pico.get()
            ack = target.simpleserial_read("r", 1, timeout=8000)
            if ack != bytearray([1]):
                raise RuntimeError(f"serialize failed: {ack}")
            return power, trig, ov

        # warm-up (AC settle)
        for _ in range(3):
            capture(rng.getrandbits(32))

        last_trig = None
        overflow_hits = 0
        t0 = time.time()
        for i in range(args.traces):
            v = rng.getrandbits(32)
            power, trig, ov = capture(v)
            traces.append(power); values.append(v); last_trig = trig
            if ov & 0x1:
                overflow_hits += 1
            if args.probe or (i + 1) % 500 == 0:
                # trigger high-region width
                hi = np.where(trig > 1000)[0]
                span = (hi[-1] - hi[0]) if len(hi) else 0
                print(f"{i+1:5d}/{args.traces} trig_hi_samples={span} "
                      f"({span*pico.interval_ns/1e3:.1f}us) rate={(i+1)/(time.time()-t0):.1f}/s",
                      flush=True)
    finally:
        if pico is not None:
            pico.close()
        try:
            target.dis()
        finally:
            scope.dis()

    T = np.vstack(traces).astype(np.float32)
    V = np.array(values, dtype=np.uint32)
    print(f"\ncaptured {T.shape}, overflow_hits={overflow_hits}", flush=True)

    # trigger window (from last trace)
    hi = np.where(last_trig > 1000)[0]
    win = (int(hi[0]), int(hi[-1])) if len(hi) else (0, T.shape[1])
    print(f"trigger window samples: {win}  (~{(win[1]-win[0])*pico.interval_ns/1e3:.1f} us)", flush=True)

    np.savez_compressed(args.out, traces=T, values=V,
                        interval_ns=np.float32(pico.interval_ns), trig_window=np.array(win))
    print(f"saved raw to {args.out}", flush=True)

    if args.probe:
        print("probe done (no SNR).", flush=True)
        # quick plot of trigger + one trace
        fig, ax = plt.subplots(2, 1, figsize=(13, 5))
        ax[0].plot(last_trig, color="green", lw=0.5); ax[0].set_title("trigger (Ch B)")
        ax[1].plot(traces[0], lw=0.4); ax[1].set_title("power (Ch A) trace 0")
        fig.tight_layout(); fig.savefig(args.plot, dpi=130); plt.close(fig)
        print(f"saved probe plot {args.plot}", flush=True)
        return

    # ---------------- leakage analysis ----------------
    hw = np.array([int(v).bit_count() for v in V], dtype=np.uint8)
    bytesV = np.stack([(V >> (8 * k)) & 0xFF for k in range(4)], axis=1).astype(np.uint8)  # (N,4) little-endian
    bytehw = np.stack([[int(b).bit_count() for b in bytesV[:, k]] for k in range(4)], axis=1).astype(np.uint8)

    snr_hw = snr(T, hw)
    snr_byte = [snr(T, bytesV[:, k]) for k in range(4)]        # 256-class per byte
    snr_bytehw = [snr(T, bytehw[:, k]) for k in range(4)]      # 9-class HW per byte

    # correlation with total HW
    Tc = T - T.mean(0)
    hc = hw.astype(np.float64) - hw.mean()
    corr_hw = (hc @ Tc) / (np.sqrt((hc @ hc) * (Tc * Tc).sum(0)) + 1e-12)

    def peak(x, w=win):
        seg = np.abs(x[w[0]:w[1]])
        j = int(np.argmax(seg)) + w[0]
        return j, float(x[j])

    print("\n=== leakage summary (peaks inside trigger window) ===", flush=True)
    j, val = peak(snr_hw); print(f"SNR[HW32]        peak @ {j:6d}  SNR={val:.4f}", flush=True)
    for k in range(4):
        j, val = peak(snr_byte[k]);   print(f"SNR[byte{k} val] peak @ {j:6d}  SNR={val:.4f}", flush=True)
    for k in range(4):
        j, val = peak(snr_bytehw[k]); print(f"SNR[byte{k} HW]  peak @ {j:6d}  SNR={val:.4f}", flush=True)
    j, val = peak(corr_hw); print(f"|corr(HW32)|     peak @ {j:6d}  corr={val:+.4f}", flush=True)

    # ---------------- plot ----------------
    fig, ax = plt.subplots(4, 1, figsize=(14, 11), sharex=True)
    ax[0].plot(last_trig, color="green", lw=0.5)
    ax[0].axvspan(win[0], win[1], color="orange", alpha=0.12)
    ax[0].set_title("Trigger (Ch B) + serialization window"); ax[0].set_ylabel("mV")
    ax[1].plot(T.mean(0), color="black", lw=0.5)
    ax[1].set_title(f"Mean power trace (N={T.shape[0]})"); ax[1].set_ylabel("mV")
    ax[2].plot(snr_hw, label="SNR HW32", color="C0", lw=0.6)
    for k in range(4):
        ax[2].plot(snr_bytehw[k], lw=0.5, alpha=0.7, label=f"SNR byte{k} HW")
    ax[2].axvspan(win[0], win[1], color="orange", alpha=0.10)
    ax[2].set_title("SNR (Hamming-weight partitions)"); ax[2].set_ylabel("SNR"); ax[2].legend(fontsize=7, ncol=5)
    for k in range(4):
        ax[3].plot(snr_byte[k], lw=0.5, alpha=0.8, label=f"SNR byte{k} value")
    ax[3].axvspan(win[0], win[1], color="orange", alpha=0.10)
    ax[3].set_title("SNR (per-byte 256-value partitions)"); ax[3].set_ylabel("SNR")
    ax[3].set_xlabel("sample"); ax[3].legend(fontsize=7, ncol=4)
    fig.tight_layout(); fig.savefig(args.plot, dpi=130); plt.close(fig)
    print(f"\nsaved analysis plot {args.plot}", flush=True)


if __name__ == "__main__":
    main()
