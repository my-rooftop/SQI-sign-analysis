#!/usr/bin/env python3
"""
SQIsign Response Leakage Capture -- PicoScope + ChipWhisperer-Lite + STM32F4 (F415)
===================================================================================
Hybrid setup, mirroring the msg_recovery_attack_hqc rig:

  ChipWhisperer-Lite : programs the target, supplies the 7.3728 MHz target clock,
                       drives SimpleSerial I/O, and toggles the trigger GPIO.
  PicoScope 3000a    : captures the actual waveform.
                         Ch A = power measurement (from CW308 SMA)
                         Ch B = trigger line (CW-Lite trigger GPIO), used as the
                                PicoScope hardware trigger source.

Measured function (RESP_TRIGGER_MODE, chosen at firmware build time):
  0 = full response core (~12 s, NOT suitable for one PicoScope block)
  1 = response ASCII serialization  <-- default, short (~2.7 ms), sanity target
  2 = sample_finalize (basis*x, denom copy, normalize)
  3 = normalize coord0 division
  4 = normalize coord0 sign multiplication
  5 = coord0 limb-loop artificial probe

For mode 1 we compute ONE response up front (command 'c', ~12 s) and then repeatedly
capture the serialization region (command 's'), which raises the trigger each time.
Use --compute-each (with a mode-0 build) to trigger the whole core instead.

Usage:
  ~/.cwvenv/bin/python3 capture_response_leakage_pico.py --traces 20
  ~/.cwvenv/bin/python3 capture_response_leakage_pico.py --traces 5 --samples 700000
"""

import argparse
import ctypes
import json
import os
import pathlib
import random
import subprocess
import sys
import time

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import chipwhisperer as cw
from chipwhisperer.capture.api.programmers import STM32FProgrammer
from picosdk.ps3000a import ps3000a as ps
from picosdk.functions import adc2mV, assert_pico_ok


APP = pathlib.Path(__file__).resolve().parent
HOST = APP / "response_split_host"
HEX = APP / "simpleserial-sqisign-response-CW308_STM32F4.hex"
VECTOR = pathlib.Path("/tmp/sqisign_response_vec.txt")
CHUNK = 48


# ============================================================
# host tool helpers (input vector + serial protocol)
# ============================================================

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
        run([str(HOST), "gen"], stdout=f)
    vector = VECTOR.read_bytes()
    r = ss_cmd(target, "x", b"", 1, timeout_ms=5000)
    if r != bytearray([1]):
        raise RuntimeError(f"reset failed: {r}")
    padded = vector
    if len(padded) % CHUNK:
        padded += b"\n" * (CHUNK - (len(padded) % CHUNK))
    for off in range(0, len(padded), CHUNK):
        r = ss_cmd(target, "a", padded[off: off + CHUNK], 1, timeout_ms=5000)
        if r != bytearray([1]):
            raise RuntimeError(f"write failed at {off}: {r}")
    print(f"loaded vector bytes={len(vector)} chunks={len(padded) // CHUNK}", flush=True)


# ============================================================
# PicoScope 3000a
# ============================================================

class Pico:
    """Ch A = power (measurement), Ch B = trigger (CW-Lite trigger GPIO)."""

    def __init__(self, samples, timebase, pre_samples, cha_range, chb_range):
        self.chandle = ctypes.c_int16()
        self.status = {}
        self.samples = samples
        self.timebase = timebase
        self.pre = pre_samples
        self.post = samples - pre_samples
        self.cha_range = cha_range
        self.chb_range = chb_range

        self.status["open"] = ps.ps3000aOpenUnit(ctypes.byref(self.chandle), None)
        try:
            assert_pico_ok(self.status["open"])
        except Exception:
            power = self.status["open"]
            if power in (282, 286):  # USB power source
                self.status["power"] = ps.ps3000aChangePowerSource(self.chandle, power)
                assert_pico_ok(self.status["power"])
            else:
                raise

        # Ch A: power, AC coupled (coupling=0) -- matches the HQC rig; DC bias would rail the input
        assert_pico_ok(ps.ps3000aSetChannel(self.chandle, 0, 1, 0, cha_range, 0))
        # Ch B: trigger, DC coupled (coupling=1)
        assert_pico_ok(ps.ps3000aSetChannel(self.chandle, 1, 1, 1, chb_range, 0))
        ps.ps3000aSetChannel(self.chandle, 2, 0, 1, 0, 0)
        ps.ps3000aSetChannel(self.chandle, 3, 0, 1, 0, 0)

        # trigger: Ch B (source=1), rising (dir=2), thresh=1000 ADC counts, autoTrig=5000ms
        assert_pico_ok(ps.ps3000aSetSimpleTrigger(self.chandle, 1, 1, 1000, 2, 0, 5000))

        interval = ctypes.c_float()
        max_samples = ctypes.c_int16()
        assert_pico_ok(ps.ps3000aGetTimebase2(
            self.chandle, timebase, samples, ctypes.byref(interval), 1,
            ctypes.byref(max_samples), 0))
        self.interval_ns = interval.value

        self.maxADC = ctypes.c_int16()
        ps.ps3000aMaximumValue(self.chandle, ctypes.byref(self.maxADC))

        print(f"PicoScope: {self.interval_ns:.3f} ns/sample, {samples} samples "
              f"(~{samples * self.interval_ns / 1e6:.3f} ms window), "
              f"pre={pre_samples}", flush=True)

    def arm(self):
        self.bufA = (ctypes.c_int16 * self.samples)()
        self.bufB = (ctypes.c_int16 * self.samples)()
        assert_pico_ok(ps.ps3000aSetDataBuffers(self.chandle, 0, self.bufA, None, self.samples, 0, 0))
        assert_pico_ok(ps.ps3000aSetDataBuffers(self.chandle, 1, self.bufB, None, self.samples, 0, 0))
        assert_pico_ok(ps.ps3000aRunBlock(
            self.chandle, self.pre, self.post, self.timebase, 1, None, 0, None, None))

    def get(self):
        ready = ctypes.c_int16(0)
        # busy-wait for capture completion (serialize triggers within a few ms)
        for _ in range(2_000_000):
            ps.ps3000aIsReady(self.chandle, ctypes.byref(ready))
            if ready.value != 0:
                break
        overflow = ctypes.c_int16()
        got = ctypes.c_int32(self.samples)
        assert_pico_ok(ps.ps3000aGetValues(
            self.chandle, 0, ctypes.byref(got), 0, 0, 0, ctypes.byref(overflow)))
        power = np.array(adc2mV(self.bufA, self.cha_range, self.maxADC), dtype=np.float32)
        trig = np.array(adc2mV(self.bufB, self.chb_range, self.maxADC), dtype=np.float32)
        return power, trig, overflow.value

    def close(self):
        ps.ps3000aStop(self.chandle)
        ps.ps3000aCloseUnit(self.chandle)


# ============================================================
# main
# ============================================================

def main():
    p = argparse.ArgumentParser(description="SQIsign response leakage capture via PicoScope")
    p.add_argument("-n", "--traces", type=int, default=20)
    p.add_argument("--samples", type=int, default=700000,
                   help="samples per trace (default covers ~2.8 ms serialization at 4 ns)")
    p.add_argument("--timebase", type=int, default=2, help="ps3000a timebase (2 -> ~4 ns with 2 chans on this unit)")
    p.add_argument("--pre-samples", type=int, default=500)
    p.add_argument("--cha-range", type=int, default=5, help="Ch A range code (5=+/-500mV)")
    p.add_argument("--chb-range", type=int, default=7, help="Ch B range code (7=+/-2V)")
    p.add_argument("--compute-each", action="store_true",
                   help="send 'c' (full core) each trace instead of 's' (serialize)")
    p.add_argument("--seed-each", action="store_true")
    p.add_argument("--fixed-seed", action="store_true")
    p.add_argument("--skip-flash", action="store_true")
    p.add_argument("--out", default="/tmp/sqisign_response_pico.npz")
    p.add_argument("--plot", default="/tmp/sqisign_response_pico_overview.png")
    args = p.parse_args()

    if not HEX.exists():
        print(f"ERROR: firmware not found: {HEX}\n"
              f"Build with: make PLATFORM=CW308_STM32F4 RESP_TRIGGER_MODE=1", file=sys.stderr)
        sys.exit(1)
    if not HOST.exists():
        run([str(APP / "build_host_tool.sh")])

    cap_cmd = "c" if args.compute_each else "s"

    # ---- ChipWhisperer-Lite: clock + programming + serial + trigger GPIO ----
    scope = cw.scope()
    target = cw.target(scope)
    pico = None
    traces, records = [], []
    try:
        scope.default_setup()  # sets 7.3728 MHz target clock + trigger IO
        if not args.skip_flash:
            print(f"programming {HEX.name} ...", flush=True)
            cw.program_target(scope, STM32FProgrammer, str(HEX))
        time.sleep(0.2)
        target.flush()
        print("banner:", target.read(timeout=1000), flush=True)

        setup_vector(target)

        if not args.compute_each:
            print("computing one response (command 'c', ~12 s) to populate last_resp ...", flush=True)
            t0 = time.time()
            r = ss_cmd(target, "c", b"", 1, timeout_ms=900000)
            if r != bytearray([1]):
                raise RuntimeError(f"initial compute failed: {r}")
            _t, initial = read_board_response(target)
            print(f"initial compute OK elapsed={time.time() - t0:.1f}s "
                  f"denom_bits={abs(initial['denom']).bit_length()} "
                  f"coord0_bits={abs(initial['coord0']).bit_length()}", flush=True)

        pico = Pico(args.samples, args.timebase, args.pre_samples,
                    args.cha_range, args.chb_range)

        # warm-up: AC coupling caps settle for the first capture or two; discard them
        for _ in range(2):
            pico.arm()
            target.flush()
            target.simpleserial_write(cap_cmd, bytearray())
            pico.get()
            target.simpleserial_read("r", 1, timeout=900000 if args.compute_each else 8000)
        print("warm-up captures done", flush=True)

        for i in range(args.traces):
            if args.seed_each:
                srng = random.Random(0 if args.fixed_seed else i)
                seed = bytes(srng.getrandbits(8) for _ in range(16))
                if ss_cmd(target, "n", seed, 1, timeout_ms=5000) != bytearray([1]):
                    raise RuntimeError(f"seed failed at trace {i}")

            pico.arm()
            target.flush()
            target.simpleserial_write(cap_cmd, bytearray())
            power, trig, overflow = pico.get()
            ack = target.simpleserial_read("r", 1, timeout=900000 if args.compute_each else 8000)
            if ack != bytearray([1]):
                raise RuntimeError(f"{cap_cmd} failed at trace {i}: {ack}")

            text, rec = read_board_response(target)
            rec["raw"] = text
            traces.append(power)
            records.append(rec)

            tspan = float(trig.max() - trig.min())
            over = "  OVERFLOW-A" if (overflow & 0x1) else ""
            print(f"{i:03d}: pmV[{power.min():7.1f},{power.max():7.1f}] "
                  f"trig_span={tspan:6.0f}mV coord0_lsw=0x{limbs32(rec['coord0'], 1)[0]:08x}{over}",
                  flush=True)

    finally:
        if pico is not None:
            pico.close()
        try:
            target.dis()
        finally:
            scope.dis()

    trace_arr = np.vstack(traces)

    # labels (sign / bit length / limb hamming weight), same scheme as the CW-ADC script
    names = ["denom", "coord0", "coord1", "coord2", "coord3"]
    labels = {}
    for name in names:
        vals = [r[name] for r in records]
        labels[f"{name}_sign"] = np.array([1 if int(v) < 0 else 0 for v in vals], dtype=np.uint8)
        labels[f"{name}_bits"] = np.array([abs(int(v)).bit_length() for v in vals], dtype=np.uint16)
        for k in range(4):
            lv = np.array([limbs32(v, k + 1)[k] for v in vals], dtype=np.uint32)
            labels[f"{name}_limb{k}"] = lv
            labels[f"{name}_limb{k}_hw"] = np.array([hw32(v) for v in lv], dtype=np.uint8)

    np.savez_compressed(args.out, traces=trace_arr,
                        records_json=json.dumps(records),
                        interval_ns=np.float32(pico.interval_ns if pico else 0),
                        **labels)
    print(f"\nsaved {trace_arr.shape} to {args.out}", flush=True)

    # ---- overview plot (sanity) ----
    last_trig = trig
    fig, ax = plt.subplots(3, 1, figsize=(14, 9))
    for i in range(min(10, len(traces))):
        ax[0].plot(traces[i], lw=0.3, alpha=0.6)
    ax[0].set_title(f"Power traces (first {min(10, len(traces))} of {len(traces)}) -- Ch A")
    ax[0].set_ylabel("mV"); ax[0].grid(alpha=0.3)

    ax[1].plot(last_trig, color="green", lw=0.4)
    ax[1].axhline(0, color="k", lw=0.3)
    ax[1].set_title("Trigger (Ch B) -- last trace"); ax[1].set_ylabel("mV"); ax[1].grid(alpha=0.3)

    ax[2].plot(trace_arr.mean(axis=0), color="blue", lw=0.6)
    ax[2].set_title(f"Mean power trace ({len(traces)} traces)")
    ax[2].set_xlabel("sample"); ax[2].set_ylabel("mV"); ax[2].grid(alpha=0.3)

    fig.tight_layout()
    fig.savefig(args.plot, dpi=130)
    plt.close(fig)
    print(f"saved overview plot to {args.plot}", flush=True)


if __name__ == "__main__":
    main()
