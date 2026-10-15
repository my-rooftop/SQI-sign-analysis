# SQIsign Response Coordinate Leakage — Experiment Summary

**Target:** SQIsign lvl1 response-quaternion generation (`sqisign_response_quat_element_core`,
the spec's `RandomEquivalentQuaternion`). Goal: recover the four response coordinates
`coord0..3` of `alpha_rsp = (coord0 + coord1·i + coord2·j + coord3·k)/denom`. The secret
key ideal is consumed in this routine, so its leakage is security-relevant.

**Rig:** PicoScope 3000a captures the power waveform (Ch A, AC-coupled); ChipWhisperer-Lite
programs the target, supplies the 7.3728 MHz clock, drives SimpleSerial, and toggles the
trigger (Ch B). Target = CW308 + STM32F415 (`PLATFORM=CW308_STM32F4`). Python env:
`~/.cwvenv/bin/python3` (chipwhisperer 6.0.0, picosdk, numpy 1.26.4).

**Coordinate sizes (lvl1):** coord0/1 ≈ 316–328 bit (10–11×32-bit limbs), coord2/3 ≈ 188–202 bit
(6–7 limbs), denom 1–2 bit, lattice_content ≈ 520 bit. Total target ≈ 1000 bit / ~130 bytes.
Each signature produces fresh coords (rejection sampling) → a real attack gets ONE trace per value.

---

## Experiment 1 — Measurement chain validation  (`capture_response_leakage_pico.py`)
PicoScope port of the CW-ADC capture. Verified: clean trigger-aligned power traces, 0 overflow,
full serialization window (~2.7 ms) captured. Fixed an AC-coupling clipping bug; added warm-up
captures (first AC-coupled trace settles). **Result: rig works.**

## Experiment 2 — Serialization leakage, known value  (`exp_serialize_leakage.py`, calibrate)
Serialize a known 32-bit value (cmd `q`). First window was too short (missed coord0 region → false
null, confirmed by permutation test). With the full window: **byte3 (MSB) SNR = 8.1**, single-trace
bit-length recovery 99%, MSB-nibble 28%. *Later shown to be a SIZE effect (MSB of a 32-bit value ≈ its magnitude).*

## Experiment 3 — Real coords, 300 traces  (`exp_real_coords_capture.py`, firmware mode 6)
Added firmware **mode 6** = per-value trigger pulses (content, denom, coord0..3) so one capture
yields 6 aligned segments. 300 real responses (12 s each). Per coord, single-trace:

| coord | bit-length SNR | 1-trace bit-length recovery | LS-limb HW SNR (null) |
|---|---|---|---|
| coord0 | 0.27 | 67% (rand 14%) | 0.11 (0.08) |
| coord1 | 0.26 | 43% | 0.16 (0.08) |
| coord2 | 0.31 | 51% | 0.26 (0.07) |
| coord3 | 0.26 | 39% | 0.11 (0.10) |

**Each coord's size leaks and is single-trace recoverable; some LS-limb data leakage above null.**

## Experiment 4 — Large-scale profiling, 5000 traces  (`exp_profile_known.py`, firmware cmd `k`)
Added firmware **cmd `k`** = load a known multi-limb response directly (bypasses the 12 s core),
enabling fast profiling with realistic values. 5000 traces. Trained per-coord templates.

| coord | bit-length (size) SNR / recovery | MS-limb DATA SNR / recovery |
|---|---|---|
| coord0 | 1.29 / **61%** | 0.01 / 4% ❌ |
| coord1 | 1.48 / **56%** | 0.01 / 5% ❌ |
| coord2 | 1.21 / **68%** | 0.01 / 4% ❌ |
| coord3 | 1.29 / **65%** | 0.01 / 5% ❌ |

**Controlled test (fix bit-length → remove size confound):** coord0 LS-limb HW SNR 0.19 (3× null),
all other limbs at noise floor; coord2 all limbs ≈ null. **Leakage is dominated by SIZE
(digit count / timing); per-limb DATA barely leaks.** Size SNR 1.5 vs data SNR 0.19 (~8×).
This explains Exp-2: "byte3 SNR 8.1" was the size effect, not data.

## Experiment 5 — Averaging gain, 3072 traces  (`exp_repeat_avg.py`)
64 fixed-size (320-bit) coord0 values × 48 repeats. Average m repeats of the SAME value:

| m | data SNR (LS-limb HW) | full-value ID acc (64 same-size candidates) |
|---|---|---|
| **1 (real-attack limit)** | 0.086 | **33.4%** (chance 1.6%) |
| 2 | 0.156 | 51.2% |
| 4 | 0.272 | 71.9% |
| 8 | 0.447 | 85.5% |
| 16 | 0.631 | 93.4% |
| 32 | 1.069 | **96.9%** |

**Averaging works strongly IF the value can be replayed** (33% → 97% for 1-of-64 ID). Data SNR ∝ ~m.
**But real SQIsign coords are unique per signature → m = 1 only.** Notably, even m=1 whole-trace
template matching gets 33% among same-size candidates — so single-trace DOES carry weak data info
(more optimistic than Exp-4's single-limb view); the hard part is *blind* recovery from 2^320.

---

## Overall conclusions
1. **Size (bit-length) of each coord leaks strongly** and is single-trace recoverable (56–68%).
   That is only ~log2 bits of info per coord.
2. **Full coordinate values do NOT leak enough for blind single-trace recovery** via power +
   serialization. Per-limb data SNR ≈ noise floor; whole-trace template gives 33% only among a
   known 64-candidate set (not blind).
3. **Averaging raises recovery to ~97%, but requires replaying the same value** — impossible in a
   normal SQIsign attack (fresh coords each signature). Relevant only for profiling or an
   RNG-fixing fault model.

## Where value recovery could still come from
- Better single-trace SNR: **EM probe** (local), CW **synchronous ADC** (no sample jitter), higher gain.
- **Horizontal** exploitation: combine the 6–11 limbs within one trace (spatial repetition, not averaging).
- **Fault / RNG reuse** to force replay, then use the averaging curve above.

## Firmware changes (experimental harness)
- `firmware/simpleserial-sqisign-response.c`: `RESP_TRIGGER_MODE==6` = per-value serialization
  trigger pulses; cmd `k` = load a known response (`[content,denom,coord0..3]` hex) into `last_resp`.
  (Command letter `k` chosen because `v`/`w`/`y` are reserved by SimpleSerial.)

## Scripts
All under `scripts/`: `capture_response_leakage_pico.py`, `exp_serialize_leakage.py`,
`exp_real_coords_capture.py` + `analyze_real_coords.py`/`plot_real_coords.py`,
`exp_profile_known.py` + `analyze_profile.py`, `exp_repeat_avg.py` + `analyze_repeat.py`.
