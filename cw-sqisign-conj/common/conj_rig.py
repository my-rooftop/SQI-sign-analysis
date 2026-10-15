#!/usr/bin/env python3
"""
Shared rig for the SQIsign Algorithm 3.18 line 2 Conjugate experiments.

Hardware, same hybrid setup as the earlier response experiments:

  ChipWhisperer-Lite : programs the STM32F415 on the CW308, supplies the
                       7.3728 MHz target clock, drives SimpleSerial, and provides
                       the trigger GPIO.
  PicoScope 3403D    : captures the waveform.
                         Ch A = power (CW308 SMA)
                         Ch B = trigger GPIO, used as the hardware trigger source

The measured region is chosen at firmware build time via CONJ_TRIG_MODE; see
conj_core.h for the mode table. Mode 2 (default) brackets Algorithm 3.18 line 2.
"""

import ctypes
import pathlib
import subprocess
import sys
import time

import numpy as np

APP = pathlib.Path(__file__).resolve().parent
HOST = APP / "conj_vector_host"
HEX = APP / "simpleserial-sqisign-conj-CW308_STM32F4.hex"
CHUNK = 48

# The board clock actually produced by scope.default_setup() on the CW-Lite:
# clkgen_div=26 off the 192 MHz clkgen gives 192e6/26 = 7.3846 MHz, not the
# 7.3728 MHz the nominal figure suggests. Measured and confirmed against the
# spectrum of a captured segment. set_target_clock() overrides it from the live
# scope object when one is available.
TARGET_CLOCK_HZ = 192_000_000.0 / 26.0


def set_target_clock(scope):
    """Take the real clock from the scope rather than trusting the default."""
    global TARGET_CLOCK_HZ
    try:
        f = float(scope.clock.clkgen_freq)
        if f > 0:
            TARGET_CLOCK_HZ = f
    except Exception:
        pass
    return TARGET_CLOCK_HZ


# ============================================================
# host tool
# ============================================================

def run(cmd, **kwargs):
    print("+", " ".join(str(c) for c in cmd), flush=True)
    return subprocess.run(cmd, check=True, **kwargs)


def ensure_host_tool():
    if not HOST.exists():
        run([str(APP / "build_host_tool.sh")])


def gen_vector(path, seed=None, retries=8):
    """Generate one Alg 3.18 input vector with conj_vector_host.

    Retries with a perturbed seed if the generator fails. Vector generation runs
    upstream quaternion code that can hit rare bad draws, and a capture that
    takes hours should not be lost to one of them.
    """
    ensure_host_tool()
    last = None
    for attempt in range(retries):
        s = seed if attempt == 0 else f"{seed}r{attempt}"
        cmd = [str(HOST), "gen"] + ([str(s)] if s is not None else [])
        try:
            with open(path, "wb") as f:
                subprocess.run(cmd, check=True, stdout=f)
            data = pathlib.Path(path).read_bytes()
            if b"vector" in data:
                if attempt:
                    print(f"  vector seed {seed!r} failed, used {s!r} instead", flush=True)
                return data
            last = "generator produced no vector section"
        except subprocess.CalledProcessError as e:
            last = f"exit {e.returncode}"
    raise RuntimeError(f"vector generation failed after {retries} attempts: {last}")


# ============================================================
# SimpleSerial protocol (see simpleserial-sqisign-conj.c)
# ============================================================

def ss_cmd(target, cmd, data=b"", pay_len=1, timeout_ms=250000):
    target.simpleserial_write(cmd, bytearray(data))
    return target.simpleserial_read("r", pay_len, timeout=timeout_ms)


def reset_target(scope, settle_s=0.4):
    """Pulse nRST on the CW308.

    Needed whenever the board is not reprogrammed (--skip-flash): the previous
    session leaves the firmware mid-protocol and scope.default_setup() re-inits
    the clock and UART under it, so the first commands read garbage.
    """
    try:
        scope.io.nrst = "low"
        time.sleep(0.05)
        scope.io.nrst = "high_z"
    except Exception as e:
        print(f"nRST pulse failed ({e}); relying on flush alone", flush=True)
    time.sleep(settle_s)


def sync(target, attempts=6):
    """Get the serial link into a known state.

    With --skip-flash the board is not reset, so whatever the previous session
    left in the UART is still queued and the first command reads garbage. Flush,
    then send the buffer-clear command until it acks cleanly.
    """
    for i in range(attempts):
        try:
            target.flush()
        except Exception:
            pass
        try:
            if ss_cmd(target, "x", b"", 1, timeout_ms=3000) == bytearray([1]):
                return True
        except Exception:
            pass
        time.sleep(0.1)
    raise RuntimeError(f"could not synchronise with the board after {attempts} attempts; "
                       f"power-cycle the CW308 or drop --skip-flash")


def board_info(target):
    r = ss_cmd(target, "i", b"", 10, timeout_ms=5000)
    f = [(r[2 * i + 1] << 8) | r[2 * i] for i in range(5)]
    return {
        "sizeof_quat_alg_elem": f[0],
        "sizeof_quat_ideal": f[1],
        "ibz_nlimbs": f[2],
        "trig_mode": f[3],
        "out_of_place": bool(f[4] & 0x100),
        "probe_reps": f[4] & 0xFF,
    }


def read_outbuf(target, timeout_ms=10000):
    r = ss_cmd(target, "l", b"", 2, timeout_ms=timeout_ms)
    outlen = (r[0] << 8) | r[1]
    out = bytearray()
    for idx in range((outlen + CHUNK - 1) // CHUNK):
        r = ss_cmd(target, "g", bytes([(idx >> 8) & 0xFF, idx & 0xFF]), CHUNK, timeout_ms=timeout_ms)
        out.extend(r)
    del out[outlen:]
    return out.decode("ascii", errors="replace")


def load_vector(target, vector_bytes):
    if ss_cmd(target, "x", b"", 1, timeout_ms=5000) != bytearray([1]):
        raise RuntimeError("board reset failed")
    padded = vector_bytes
    if len(padded) % CHUNK:
        padded += b"\n" * (CHUNK - (len(padded) % CHUNK))
    for off in range(0, len(padded), CHUNK):
        if ss_cmd(target, "a", padded[off:off + CHUNK], 1, timeout_ms=5000) != bytearray([1]):
            raise RuntimeError(f"vector write failed at offset {off}")


def seed_board(target, seed16):
    if ss_cmd(target, "n", seed16, 1, timeout_ms=5000) != bytearray([1]):
        raise RuntimeError("board seed failed")


def read_profile(target):
    ss_cmd(target, "p", b"", 1, timeout_ms=10000)
    text = read_outbuf(target)
    prof = {}
    for line in text.split():
        if "=" in line:
            k, v = line.split("=", 1)
            try:
                prof[k] = int(v)
            except ValueError:
                pass
    return prof


def parse_result(text):
    """Parse the 'c' command output.

    Layout (see serialize_result() in simpleserial-sqisign-conj.c):
        "1"
        pre.denom, pre.coord0..3                   hex values
        4 limb lines: "<bitlen> <nlimbs> <limb0> ... <limbN-1>"   coord0..3
        post.denom, post.coord0..3                 hex values
        norm, skchall.norm, resp.denom, resp.coord0..3
    """
    lines = [ln.strip() for ln in text.splitlines() if ln.strip()]
    if not lines or lines[0] != "1":
        raise ValueError(f"board reported failure: {text[:120]!r}")
    i = 1

    def take_hex():
        nonlocal i
        v = int(lines[i], 16)
        i += 1
        return v

    pre = {"denom": take_hex(), "coord": [take_hex() for _ in range(4)]}

    limbs = []
    for _ in range(4):
        toks = lines[i].split()
        i += 1
        bitlen, n = int(toks[0]), int(toks[1])
        words = [int(t, 16) for t in toks[2:2 + n]]
        if len(words) != n:
            raise ValueError(f"limb line short: expected {n}, got {len(words)}")
        limbs.append({"bitlen": bitlen, "nlimbs": n, "limbs": words})

    post = {"denom": take_hex(), "coord": [take_hex() for _ in range(4)]}
    norm = take_hex()
    skchall_norm = take_hex()
    resp = {"denom": take_hex(), "coord": [take_hex() for _ in range(4)]}

    return {"pre": pre, "pre_limbs": limbs, "post": post,
            "norm": norm, "skchall_norm": skchall_norm, "resp": resp}


def check_conjugate(res):
    """Line 2 must negate coord[1..3] and leave coord[0] and denom alone."""
    problems = []
    if res["post"]["denom"] != res["pre"]["denom"]:
        problems.append("denom changed")
    if res["post"]["coord"][0] != res["pre"]["coord"][0]:
        problems.append("coord[0] changed")
    for k in (1, 2, 3):
        if res["post"]["coord"][k] != -res["pre"]["coord"][k]:
            problems.append(f"coord[{k}] is not the negation of the input")
    return problems


# ============================================================
# labels
# ============================================================

def hw(x, width=32):
    return int(x & ((1 << width) - 1)).bit_count()


def sign_ext_limbs(words):
    """Number of top limbs that are pure sign extension (0x00000000/0xffffffff)."""
    if not words:
        return 0
    ext = 0xFFFFFFFF if (words[-1] >> 31) else 0x00000000
    k = 0
    for w in reversed(words):
        if w != ext:
            break
        k += 1
    return k


# ============================================================
# PicoScope 3403D
# ============================================================

class Pico:
    """Ch A = power, Ch B = trigger GPIO (hardware trigger source)."""

    def __init__(self, samples, timebase, pre_samples, cha_range, chb_range,
                 auto_trigger_ms=0):
        from picosdk.ps3000a import ps3000a as ps
        from picosdk.functions import adc2mV, assert_pico_ok
        self._ps, self._adc2mV, self._ok = ps, adc2mV, assert_pico_ok

        self.chandle = ctypes.c_int16()
        self.samples = samples
        self.timebase = timebase
        self.pre = pre_samples
        self.post = samples - pre_samples
        self.cha_range = cha_range
        self.chb_range = chb_range

        status = ps.ps3000aOpenUnit(ctypes.byref(self.chandle), None)
        try:
            assert_pico_ok(status)
        except Exception:
            if status in (282, 286):  # running from USB power
                assert_pico_ok(ps.ps3000aChangePowerSource(self.chandle, status))
            else:
                raise

        # Ch A power, AC coupled: the DC bias on the CW308 shunt would rail the input.
        assert_pico_ok(ps.ps3000aSetChannel(self.chandle, 0, 1, 0, cha_range, 0))
        # Ch B trigger, DC coupled.
        assert_pico_ok(ps.ps3000aSetChannel(self.chandle, 1, 1, 1, chb_range, 0))
        ps.ps3000aSetChannel(self.chandle, 2, 0, 1, 0, 0)
        ps.ps3000aSetChannel(self.chandle, 3, 0, 1, 0, 0)

        # Rising edge on Ch B. auto_trigger_ms=0 means wait indefinitely, which is
        # what this target needs: Algorithm 3.18 line 1 runs for seconds before the
        # trigger fires, and any non-zero auto-trigger would return a window of
        # line 1 instead of line 2.
        assert_pico_ok(ps.ps3000aSetSimpleTrigger(
            self.chandle, 1, 1, 1000, 2, 0, auto_trigger_ms))

        interval = ctypes.c_float()
        max_samples = ctypes.c_int16()
        assert_pico_ok(ps.ps3000aGetTimebase2(
            self.chandle, timebase, samples, ctypes.byref(interval), 1,
            ctypes.byref(max_samples), 0))
        self.interval_ns = interval.value
        self.samples_per_cycle = 1e9 / TARGET_CLOCK_HZ / self.interval_ns

        self.maxADC = ctypes.c_int16()
        ps.ps3000aMaximumValue(self.chandle, ctypes.byref(self.maxADC))

        # Allocate the capture buffers once and keep handing the driver the same
        # pointers. Allocating fresh ones per arm() and re-registering them made
        # the capture rate decay trace after trace until the scope stopped
        # triggering at all, which only shows on a fast target where thousands of
        # arms happen in a few minutes.
        self.bufA = (ctypes.c_int16 * self.samples)()
        self.bufB = (ctypes.c_int16 * self.samples)()
        assert_pico_ok(ps.ps3000aSetDataBuffers(self.chandle, 0, self.bufA, None, self.samples, 0, 0))
        assert_pico_ok(ps.ps3000aSetDataBuffers(self.chandle, 1, self.bufB, None, self.samples, 0, 0))

        print(f"PicoScope: {self.interval_ns:.3f} ns/sample "
              f"({self.samples_per_cycle:.1f} samples per "
              f"{TARGET_CLOCK_HZ / 1e6:.4f} MHz target cycle), "
              f"{samples} samples = {samples * self.interval_ns / 1e6:.3f} ms window, "
              f"pre={pre_samples}", flush=True)

    def arm(self, settle_s=0.01):
        """Arm for a block capture.

        The short settle matters for fast targets: RunBlock returns before the
        scope will actually honour a trigger, and a target that fires 20 ms later
        can slip through the gap. It never showed on the Algorithm 3.18 target,
        where the trigger arrives seconds after arming.
        """
        ps = self._ps
        self._ok(ps.ps3000aRunBlock(
            self.chandle, self.pre, self.post, self.timebase, 1, None, 0, None, None))
        if settle_s:
            time.sleep(settle_s)

    def get(self, timeout_s=600.0):
        ps = self._ps
        ready = ctypes.c_int16(0)
        deadline = time.time() + timeout_s
        while time.time() < deadline:
            ps.ps3000aIsReady(self.chandle, ctypes.byref(ready))
            if ready.value != 0:
                break
            time.sleep(0.001)
        if ready.value == 0:
            ps.ps3000aStop(self.chandle)
            raise TimeoutError(f"PicoScope did not trigger within {timeout_s:.0f}s")
        overflow = ctypes.c_int16()
        got = ctypes.c_int32(self.samples)
        self._ok(ps.ps3000aGetValues(
            self.chandle, 0, ctypes.byref(got), 0, 0, 0, ctypes.byref(overflow)))
        power = np.array(self._adc2mV(self.bufA, self.cha_range, self.maxADC), dtype=np.float32)
        trig = np.array(self._adc2mV(self.bufB, self.chb_range, self.maxADC), dtype=np.float32)
        # Stop before the next RunBlock. Driving RunBlock repeatedly without it
        # degrades the capture rate over a few hundred to a few thousand arms
        # until the scope stops honouring the trigger entirely -- only visible on
        # a fast target, where that many arms happen in minutes.
        ps.ps3000aStop(self.chandle)
        return power, trig, overflow.value

    def close(self):
        self._ps.ps3000aStop(self.chandle)
        self._ps.ps3000aCloseUnit(self.chandle)


# ============================================================
# trigger analysis
# ============================================================

def trigger_edges(trig, threshold_mv=None):
    """Return (rising, falling) sample indices of the trigger line."""
    t = np.asarray(trig, dtype=np.float32)
    if threshold_mv is None:
        lo, hi = float(t.min()), float(t.max())
        threshold_mv = lo + 0.5 * (hi - lo)
    high = t > threshold_mv
    d = np.diff(high.astype(np.int8))
    return np.flatnonzero(d == 1) + 1, np.flatnonzero(d == -1) + 1


def trigger_pulses(trig, threshold_mv=None):
    """Pulse list as (start, end, width_samples), pairing each rise with the next fall."""
    rise, fall = trigger_edges(trig, threshold_mv)
    pulses = []
    for r in rise:
        nxt = fall[fall > r]
        if nxt.size:
            pulses.append((int(r), int(nxt[0]), int(nxt[0] - r)))
    return pulses


def expected_pulses(trig_mode, limbs_per_coord):
    """How many trigger pulses each mode should produce, and what each one covers."""
    if trig_mode == 0:
        return 1, ["whole Algorithm 3.18"]
    if trig_mode == 1:
        return 1, ["line 1 EquivalentCoprimeIdeal"]
    if trig_mode == 2:
        return 1, ["line 2 Conjugate (all of it)"]
    if trig_mode == 3:
        return 3, [f"line 2 ibz_neg coord[{k}]" for k in (1, 2, 3)]
    if trig_mode == 4:
        return 5, ["line 2 ibz_copy denom", "line 2 ibz_copy coord[0]"] + \
                  [f"line 2 ibz_neg coord[{k}]" for k in (1, 2, 3)]
    if trig_mode == 5:
        return 1, ["line 5 IdealMultiplication"]
    if trig_mode == 6:
        return 1, ["line 5 RandomBoundedElement"]
    if trig_mode == 7:
        return 3, [f"limb-walk probe coord[{k}]" for k in (1, 2, 3)]
    if trig_mode == 8:
        n = limbs_per_coord or 0
        return n, [f"line 2 ibz_neg coord[1] limb {k}" for k in range(n)]
    return None, []
