# SQIsign Algorithm 3.18 side-channel experiments (Round 3)

Upstream `the-sqisign.git` **f417ebd** / the 2026-09-01 spec
(`sqisign-20260901.pdf`), parameter set `p324_3`, STM32F415 (CW308) +
PicoScope 3403D.

Two **different operations** inside Algorithm 3.18 `SampleQuaternionResponse` were
measured. Not two algorithms — two lines of the same algorithm.

| Part | Paper | Target | Result |
|---|---|---|---|
| [`part1_alg318_line2/`](part1_alg318_line2/LINE2_EXPERIMENT.md) | **5.3.1** | line 2 `Conjugate(α_sk,chl)` | 4 sign bits of 614, 0 MSBs |
| [`part2_alg318_line5/`](part2_alg318_line5/RESPONSE_EXPERIMENT.md) | **5.3.2** | line 5 `RandomBoundedElement` (α_rsp) | 4 sign bits of 746, 0 MSBs |

[`extra_alg323_msk/`](extra_alg323_msk/) is separate work against `Msk` at
**Algorithm 3.23 `Sign` line 12**. `Msk` is a fixed secret so CPA does apply
(`Msk[0][1]`, 128/128 bits from 2,000 non-profiled traces), but the full matrix
comes out at 320/512 and reproduces only 2 runs in 3, so it is **unfinished**. It
does not go into section 5.3 of the paper.

The firmware is a **single binary** and `CONJ_TRIG_MODE` selects the target at build
time. The C sources, the rig and the analysers are therefore shared by both parts
and live at the top level and in `common/`.

---

## Everything below is about part 1 (line 2)

For the line 5 (5.3.2) measurement and analysis see
[`part2_alg318_line5/RESPONSE_EXPERIMENT.md`](part2_alg318_line5/RESPONSE_EXPERIMENT.md).

This package targets **Algorithm 3.18 SampleQuaternionResponse, line 2** of the
2026-09-01 spec (`sqisign-20260901.pdf`).

```
Algorithm 3.18 SampleQuaternionResponse                                    (spec p.21)
 1: I'_sk,chl, α_sk,chl <- EquivalentCoprimeIdeal(I_sk,chl, 0, BitLenSkChall, prng_def)
 2: α_sk,chl <- Conjugate(α_sk,chl)                                    <- the target
 3: if GCD(IdealNorm(I_com), IdealNorm(I'_sk,chl)) != 1 then raise
 5: q_rsp, α_rsp <- RandomBoundedElement(IdealMultiplication(I'_sk,chl, I_com), 2^(e_rsp-1), prng_res)
```

## Spec to code

| Spec | Vendored code |
|---|---|
| Alg 3.18 as a whole | `protocol.c` `quat_response_element()` |
| Alg 3.18 line 1 (`EquivalentCoprimeIdeal`) | `ideal.c:846` `quat_ideal_small_equivalent_coprime()` |
| **Alg 3.18 line 2 (`Conjugate`)** | **`protocol.c:37` `quat_alg_conj(sk_chall_quat, sk_chall_quat)`** |
| Alg 4.38 `Conjugate` | `algebra.c:168` `ibz_copy` x2 (denom, coord0) + `ibz_neg` x3 (coord1,2,3) |
| The actual arithmetic | `mp.c:469` `ibz_neg()` — limb-wise two's-complement negation loop |
| Alg 4.74 `IdealMultiplication` | `internal.h:734` `quat_ideal_mul_O0()` |
| Alg 4.83 `RandomBoundedElement` | `internal.h:571` `quat_lattice_sample_from_ball()` |

Upstream hoists the GCD check of Alg 3.18 line 3 out to `sign.c:276`.

## Why this is measurable at all

The 2026 implementation replaced mini-gmp with its own `mp` module, and `ibz_t`
became **fixed-width two's complement**:

```c
struct ibz_t { int bitlen; digit_t limbs[IBZ_NLIMBS]; };   // p324_3, RADIX=32 -> 60

void ibz_neg(ibz_t *neg, const ibz_t *a) {                 // mp.c:469
    int n = NUM_LIMBS(a->bitlen);
    digit_t carry = 1;
    for (int i = 0; i < n; i++) {
        ddigit_t t = (ddigit_t)~a->limbs[i] + carry;
        neg->limbs[i] = t;                                 // every limb is read and written
        carry = t >> NUM_BITS_LIMB;
    }
    neg->bitlen = a->bitlen;
}
```

Under the old mini-gmp, `mpz_neg` was `mpz_set` plus a sign-bit flip, and `mpz_set`
returns early when `r != x` (`mini-gmp.c:1533`), so **an in-place conjugate touched
no limbs at all**. The new implementation runs `LDR → MVN → ADCS → STR` over every
limb even in place, so per-limb HW/HD leakage is there to be measured. The loop count
`n = NUM_LIMBS(bitlen)` depends only on the **declared bound**, not on the value, so
slot positions are identical across traces — and there is no branch.

## Size of the target (measured on the host, `conj_vector_host selftest`)

The table below covers only the three **negated** coordinates; `coord[0]` is copied,
not negated, and is counted separately further down.

| | coord[1] (`i`) | coord[2] (`j`) | coord[3] (`ij`) |
|---|---|---|---|
| Declared width `bitlen` | 401 (fixed) | 401 (fixed) | 401 (fixed) |
| Limbs walked by `ibz_neg` | 13 | 13 | 13 |
| Actual value size | 230–235 bits | 69–73 bits | 69–73 bits |
| Sign-extension limbs | 5 | 10 | 10 |
| **Informative limbs** | **8** | **3** | **3** |

- Total limbs `ibz_neg` walks for these three: **39 = 1248 bits of width**, identical
  on every trace
- Of those, **14 limbs = 448 bits** carry information
- **Actual recovery target for these three ≈ 376 bits** (234 + 71 + 71)

Counting `coord[0]` as well brings the total to the 614 bits used throughout the rest
of this document.

coord[2] and coord[3] are small because the norm form is `a² + b² + p(c² + d²)` with
log₂p = 324, so c and d are about √p times smaller than a and b
(234 − 71 = 163 ≈ 324/2).

**`coord[0]` (`a`, ~234 bits) and denom are `ibz_copy` (`memmove`), not negate.**
Board measurement shows an in-place `ibz_copy` costs **144 cycles (11 per limb)** — it
is not a no-op, so the `a` coordinate also crosses the bus despite not being negated,
and is a leakage candidate. `CONJ_TRIG_MODE=4` puts separate triggers on these two
regions.

### The `bounds` line in vector files

An `ibz_t` carries a declared bit bound which downstream operations inherit, but
`ibz_set_from_str()` (`mp.c:311`) re-derives the bound from the parsed magnitude.
That difference propagates into α_sk,chl and changes how many limbs `ibz_neg` walks at
line 2 (401/13 becomes 319/10 across a hex round-trip). The vector file therefore
carries the bounds of all 14 inputs on a `bounds` line, and both host and firmware
restore them with `ibz_set_bound()` after parsing. Without this the board executes a
narrower region than a real signing device. The firmware reports whether it was
applied via `bounds_applied` in the `p` command.

## Layout

```
vendor/sqisign/p324_3/       upstream embedded sources (the-sqisign.git f417ebd,
                             output of scripts/gen_pqm4_sources.sh ref, plus fips202)
vendor/sqisign/host-helpers/ host-only upstream test helpers (see PROVENANCE)
conj_core.{c,h}              instrumented Algorithm 3.18. upstream
                             quat_response_element() ported line by line, with only
                             triggers and profiling added
simpleserial-sqisign-conj.c  CW308_STM32F4 / STM32F415 simpleserial target
cw_randombytes.c             board RNG shim
host/conj_vector_host.c      vector generation, KAT, ground-truth output (PC)
build_host_tool.sh           host tool build
Makefile                     firmware build

common/                      Python shared by both parts
  conj_rig.py                shared PicoScope 3403D + CW-Lite rig
  check_trigger.py           trigger alignment check
  analyze_conj_leakage.py    TVLA + SNR + correlation
  recover_conj_template.py   profiled templates + per-bit LDA + shuffled-label control
  msb_recover.py             leading-T-bit recovery measurement
  dl_recover.py              deep-learning recovery (MLP)

part1_alg318_line2/          5.3.1 — line 2 Conjugate
  LINE2_EXPERIMENT.md        experiment report
  capture_conj_pico.py       capture — random / fixed / tvla
  capture_line2_pico.py      fast injected-path capture
  figures/                   9 figures

part2_alg318_line5/          5.3.2 — line 5 RandomBoundedElement
  RESPONSE_EXPERIMENT.md     experiment report
  capture_resp_pico.py       capture — random / tvla
  run_resp_analysis.sh       runs all five analysis stages
  make_resp_figures.py       figure generation
  figures/                   3 figures

extra_alg323_msk/            Alg 3.23 line 12 Msk (outside section 5.3, unfinished)
  capture_chall_pico.py
  cpa_chall{,_full,_matrix}.py
```

The firmware links only quaternion + mp + PRNG (ROM 64 KB / 1 MB, static RAM
20.8 KB / 128 KB). No EC, GF or isogeny code is included — the challenge ideal is
built on the host and sent across as data.

## Trigger modes

Selected at build time with `CONJ_TRIG_MODE=n`.

| Mode | Trigger region | Pulses |
|---:|---|---:|
| 0 | whole Algorithm 3.18 | 1 |
| 1 | line 1 `EquivalentCoprimeIdeal` | 1 |
| **2** | **whole of line 2 `Conjugate` (default)** | **1** |
| 3 | line 2, `ibz_neg` on coord[1]/[2]/[3] separately | 3 |
| 4 | line 2, denom copy / coord[0] copy / coord[1..3] neg | 5 |
| 5 | line 5 `IdealMultiplication` | 1 |
| 6 | line 5 `RandomBoundedElement` | 1 |
| 7 | artificial α_sk,chl limb-walk probe (`CONJ_PROBE_REPS`) | 3 |
| 8 | line 2, one pulse **per limb** of coord[1] | 13 |
| 10–14 | Alg 3.23 line 12 `ibz_mat_2x2_eval` (`extra_alg323_msk/`) | varies |
| 15 | line 5, whole formation (152,302 cycles — exceeds the window) | 1 |
| 16 | line 5, `ibz_mat_4x4_eval` only (39,085 cycles — exceeds the window) | 1 |
| 17 | line 5, per-coordinate accumulate (exceeds the window) | 4 |
| **18** | **line 5, per-coordinate limb write-out — the mode 5.3.2 uses** | **4** |

Mode 8 pulses on every iteration of the `ibz_neg` loop, so the capture is segmented
by limb index and alignment uncertainty disappears. Mode 7 is not a real code path —
it exists to establish the leakage ceiling and a profiling baseline.

`CONJ_OUT_OF_PLACE=1` separates the conjugate's source and destination. The default
of 0 is the in-place form from the spec and upstream, and for the reason above the
in-place form leaks perfectly well.

## Running order

### 0. Fix the sizes on the host, and run the KAT
```bash
./build_host_tool.sh
./conj_vector_host selftest 8          # checks invariants, measures limb/bit counts
./conj_vector_host gen myseed > vec.txt
./conj_vector_host run vec.txt         # prints α_sk,chl before/after and the raw limbs
```

`run` verifies: `Norm(α_rsp) = q_rsp · n(I'_sk,chl) · n(I_com)`,
`bitsize(q_rsp) <= e_rsp − 1`, that `Conjugate` is an involution, and the GCD
condition of line 3.

### 1. Build the firmware
```bash
make clean PLATFORM=CW308_STM32F4
make PLATFORM=CW308_STM32F4 CONJ_TRIG_MODE=2
```
Always `make clean` first when changing compile options.

### 2-3. Verify trigger alignment
```bash
~/.cwvenv/bin/python3 common/check_trigger.py --traces 3
```
What it checks:
- that the board's Conjugate is arithmetically correct (only coord[1..3] flip sign;
  coord[0] and denom unchanged)
- that the trigger pulse count matches what the mode should produce
- that the trigger-high width, converted to target cycles, matches the **DWT cycle
  count** the firmware reports
- that cycles per limb are consistent with the cost of one `ibz_neg` iteration
- that the trigger voltage swing is sufficient

Nothing after this step means anything until it passes. Alg 3.18 line 1 runs first
for several seconds, so the PicoScope waits indefinitely on the hardware trigger
**with auto-trigger disabled** (`auto_trigger_ms=0`).

### 4. Capture
```bash
# for SNR and profiling: a fresh input and a fresh seed on every trace
~/.cwvenv/bin/python3 part1_alg318_line2/capture_conj_pico.py --mode random --traces 200 --skip-flash

# for TVLA: fixed and random interleaved trace by trace
~/.cwvenv/bin/python3 part1_alg318_line2/capture_conj_pico.py --mode tvla --traces 200 --skip-flash
```
`--crop` stores only the triggered region. Each trace costs one full run of Alg 3.18,
most of which is line 1, so measure the per-run time with `check_trigger.py` first and
budget accordingly.

Interleaving fixed and random prevents temperature and supply drift from being read as
leakage.

### 5. Quantify the leakage
```bash
~/.cwvenv/bin/python3 common/analyze_conj_leakage.py \
    --tvla /tmp/sqisign_conj_tvla.npz \
    --snr  /tmp/sqisign_conj_random.npz
```

- **TVLA** — fixed vs random Welch t. `|t| > 4.5` declares leakage. The same output
  also carries a **fixed vs fixed** half-split; if that one crosses 4.5 the rig is
  drifting and the fixed-vs-random result cannot be trusted.
- **SNR** — per-limb SNR with the limb's HW as the class. HW is binomial(32, ½) so the
  extreme classes are nearly empty and the SNR argmax is unstable; by default HW is
  grouped into 8 bins (`--hw-bins`) and classes with fewer traces than `--min-class`
  (default 5) are dropped.
- **correlation** — per-sample Pearson r against `HW(limb)` and `HW(~limb)`. Its argmax
  is more stable than SNR's, so the point of interest is taken at the `|corr|` peak and
  SNR is read at that sample. Two independent statistics spiking together is the
  evidence wanted. The noise floor is reported as `4/√N`.
- Sign-extension limbs are excluded as `label constant, no SNR`.

Rule of thumb: if the `|corr|` peak clears the `4/√N` floor and SNR is roughly 0.1 or
above, CPA against a fixed value would converge in hundreds of traces; 0.01–0.1 means
thousands; below that, improve the measurement itself (probe, range, averaging) before
attempting recovery.

## Board measurements (STM32F415 @ 7.3846 MHz)

| Region | Cycles | Share |
|---|---:|---:|
| Whole Algorithm 3.18 | 73–92 M (10–12.5 s) | 100% |
| line 1 `EquivalentCoprimeIdeal` | 43–50 M | ~53% |
| line 5 `RandomBoundedElement` | 41–42 M | ~46% |
| line 5 `IdealMultiplication` | 0.19 M | 0.2% |
| `Intersect` (before line 1) | 0.25 M | 0.3% |
| **line 2 `Conjugate`** | **1,121 (152 µs)** | **0.0013%** |

Inside line 2: `coord[0]` copy 144 cycles, `coord[1]/[2]/[3]` negate 247 cycles each
→ **19.0 cycles per limb** x 13 limbs.

The mode 2 trigger starts exactly at the pre-sample position and its width has
**sd = 0.0 across traces** (34,739 samples = 1024 cycles at 4 ns). The ratio against
the 1121 cycles DWT reports is 0.91; the difference is bitlen/limb profiling code that
sits outside the trigger but inside the instrumented region. No alignment work is
needed.

## Verification status

- Host Algorithm 3.18: passes every upstream invariant (`selftest: OK`)
- ARM firmware: builds without warnings (ROM 64 KB, static RAM 20.8 KB)
- **Board**: `check_trigger.py` passes — conjugate arithmetic verified, pulse count,
  width against DWT, cycles per limb, and bitlen 401 / 13 limbs reproduced
  (`bounds_applied=1`)
- **Signal**: zero overflow at Ch A ±100 mV; the spectrum shows the 7.384 MHz target
  clock and its 14.8 / 29.5 / 36.9 / 44.3 MHz harmonics, and folding at the limb
  period (644 samples) produces a reproducible loop signature in the average
- TVLA/SNR/correlation numeric paths: validated on synthetic data with known HW
  leakage injected. The correlation locator recovered the time position of 7 of 8
  injected limbs exactly, and all 25 sign-extension limbs were correctly excluded

## Measured leakage results

### TVLA — 60 fixed / 60 random, CONJ_TRIG_MODE=2

```
fixed vs random : max|t| = 9.92 at sample 13757 (cycle ~406)
                  2647 of 34735 samples exceed 4.5 (7.62%)
fixed vs fixed  : max|t| = 4.23, 0 samples over threshold  OK
verdict: LEAKAGE DETECTED
```

The control (fixed vs fixed) does not cross the threshold, so this is data-dependent
leakage rather than rig drift. 120 traces were enough to detect it.

The t-trace is not uniform but **clusters into 4 bands**. Band centres sit at roughly
175 / 428 / 635 / 893 cycles, about 250 cycles apart, matching the measured 247 cycles
per coordinate. The four bands appear to correspond to the `coord[0]` copy and the
`coord[1]`/`[2]`/`[3]` negations, with the quiet stretches between them being
profiling calls and loop setup.

The `coord[0]` band also reaches |t| ≈ 9.9, consistent with the disassembly analysis
that the `a` coordinate leaks even though it is only copied by `memmove` rather than
negated.

The band-to-coordinate mapping is confirmed by measurement in the SNR results below.

### Per-limb SNR / correlation — 200 random traces, CONJ_TRIG_MODE=2

Excluding the 30 sign-extension limbs, **all 22 informative limbs clear the 0.283
noise floor.**

| Coordinate | Limbs | \|corr(HW)\| | SNR@POI | Position (cycle) |
|---|---|---|---|---|
| coord[0] (`a`, `ibz_copy`) | 0–6 | 0.46 – 0.80 | 0.45 – 3.69 | 121 → 143 |
| coord[1] (`b`, `ibz_neg`) | 0–6 | 0.57 – 0.72 | 0.90 – 2.65 | 262 → 358 |
| coord[2] (`c`, `ibz_neg`) | 0–1 | 0.40, 0.49 | 0.86, 0.92 | 524, 540 |
| coord[3] (`d`, `ibz_neg`) | 0–1 | 0.44, 0.50 | 0.64, 0.71 | 784, 802 |

**The limbs march in time order.** coord[1]'s limbs 0→6 land in sequence at samples
8862 → 12106, about 540 samples ≈ 16 cycles apart (consistent with the measured
19 cycles/limb). Per-limb separation in time is confirmed by data.

That also fixes the per-coordinate regions: coord[0] copy = cycles 121–161,
coord[1] neg = 262–438, coord[2] = 524–652, coord[3] = 784–834. The 262 / 260 cycle
gaps between coordinate starts match the measured 247 cycles plus profiling overhead.

#### Two cautions when reading these numbers

**The top limb of each coordinate (\|corr\| 0.94–0.97) is sign leakage.** coord[0]
limb 7, coord[1] limb 7, coord[2] limb 2 and coord[3] limb 2 all read 0.94–0.97, but
their peak samples (5458 / 14830 / 22065 / 28227) coincide with the peaks of the
sign-extension limbs. The values fall short of the limb boundary, so the top limb is
mostly sign bits. **The genuine per-limb HW leakage is the 0.40–0.80 of the 18 lower
limbs.** `classify_limb()` in `analyze_conj_leakage.py` sorts limbs into
constant / sign / informative and drops sign from the ranking — without that, the
sign-extension limbs fill the top of the table at an identical 0.97.

**The "13 traces per limb" estimate is optimistic.** It applies the Mangard formula at
median r = 0.618, which measures *separation at the correct hypothesis*. A real attack
enumerating 2^32 hypotheses needs far more, and splitting into bytes shrinks the
distinguisher margin.

Even so, |corr| 0.40–0.80 and SNR 0.6–3.7 are well above the usual "exploitable"
threshold of SNR ≈ 0.1.

## Recovery experiment (common/recover_conj_template.py)

### Why a template attack and not CPA

α_sk,chl is redrawn for every signature (Alg 3.18 line 1). There is no fixed secret to
average across traces, so CPA does not apply. Collecting many traces serves
**profiling**; the real attack has to recover from the **single trace** belonging to
the signature under attack.

This experiment profiles on the same board with ground truth known — the easiest
condition an attacker could face — so the results should be read as an upper bound.

### The leakage function is `HW(32-bit word)` and nothing else

| Leakage model | 200 traces | 20,000 traces |
|---|---|---|
| `HW(word)` | **0.40 – 0.80**, certain | certain |
| Individual bits (binary LDA) | 0 of 576 above | **481/576 above** (mean 54.5%, baseline 50.7%) |
| Per-byte HW (partial correlation) | 3 of 72 marginal | — |
| HD between consecutive limbs | 2 of 17 marginal | — |

**Bit-level leakage does exist.** The 0/576 at 200 traces was insufficient power (only
40 attack traces). Remeasured at 20,000 traces, 481/576 significantly beat the
baseline, and the network sees the same thing. So the correct statement is not
"impossible in principle" but **"present and meaningless"**.

Meaningless in size, that is. 54.5% accuracy, even perfectly calibrated, is
`1 − H(0.545) = 0.0058 bits` per bit, or at most 3.3 bits over all 576. The measured
network figure is **0.044 bits** per 32-bit limb, about 0.8 bits over 18 limbs.
Accuracy alone easily reads as "recovered", so the verdict has to be made on
information.

The mechanism appears to be the carry chain in `ibz_neg` — the carry leaving a limb
depends on whether that limb is zero, which is not a function symmetric in the bits.
A small bit-level dependency riding on top of HW is therefore expected.

Per-byte HW showing raw correlation at half of `HW(word)` is a mathematical
consequence of `HW(word) = Σ HW(byte)`, not independent byte leakage. Removing
`HW(word)` by partial correlation makes it disappear.

### Information-theoretic ceiling

```
maximum extractable per limb = H(HW of a 32-bit word) = 3.55 bits
18 clean limbs x 3.55                =  64 bits   <- the ceiling
measured (templates, 20,000 traces)  =  10 bits
measured (network, 20,000 traces)    =  29 bits
secret to recover                    = 610 bits
search space remaining               ~ 546 bits
```

**Single-trace recovery of α_sk,chl is not possible.** The leakage itself is strong
(TVLA max|t| = 9.92) but cannot be used for recovery.

The reason is structural: **`Conjugate` does not mix the secret with a known value.**
There is no divide-and-conquer structure like AES SubBytes, where a known plaintext
plus a guessed key byte predicts an intermediate — it is plain bitwise complement. So
all a trace gives up is `HW(limb)`, and 3.55 bits per limb is the ceiling. More traces
do not raise it — they only improve the accuracy of the HW estimate.

### Two analysis traps

**Sign-extension limbs**: 30 of the 52 limbs are either `0x00000000` or `0xffffffff`,
carrying 1 sign bit per coordinate and 4 bits in total. They share an identical label
and show |corr| 0.97 at the same sample; without filtering, this becomes the
overestimate "5 limbs leak at 0.98". `classify_limb()` sorts them into
constant / sign / informative.

**Max-over-samples threshold**: `4/√N` is the threshold for **one pre-chosen sample**.
Compared against the **maximum** over 34,734 samples it declares leakage even in pure
noise. The correct criterion is `√(2 ln S)/√N`, which at N=200, S=34734 is **0.323**,
not 0.283. The analyser prints both and judges on the latter.

## How many bits line 2 gives up (measured at 20,000 traces)

### Fast path

Running the whole of Alg 3.18 on the board costs 13.5 s/trace, which rules out a
statistically meaningful campaign. Only line 2 is measured, so the **host runs the
real Alg 3.18** to produce α (`conj_vector_host alphas`, 2000 of them in 7.4 s) and
the board executes line 2 alone (`u`/`z`/`q` commands). The executed code and the
trigger placement are identical and only the surrounding instructions differ; capture
start verifies this with `alpha loaded exactly, Conjugate correct`.
**238 ms/trace — 57 times faster.**

### Results

Profiled templates, 16,000 profiling / 4,000 attack (held out), single-trace attack.

| Target | Accuracy | Baseline | z | Verdict |
|---|---|---|---|---|
| **Sign (1 bit per coordinate)** | **100.0%** | 50% | +62 | **fully recovered** |
| HW of the MSB limb | 26.6–28.6% | 11–14% | +23 to +33 | recovered |
| HW of the LSB limb (control) | 19.7–21.9% | 14% | +11 to +15 | recovered |
| magnitude bit-length | 27.4% | 27.5% | −0.0 | at baseline |
| **Leading T bits (T = 1…8)** | same as baseline | | | **not recovered** |

Total HW information is **10.35 bits across 18 limbs** (ceiling 64 bits).

```
recovered:      4 sign bits (completely) + 28.75 bits of HW information (network)  ->  about 33 bits
not recovered:  0 high-order bits,  0.04 bits/limb on individual bits (meaningless)
whole secret:   614 bits
```

HW information shrinks the search space by a factor of 2^28.75; it does not name
bits. The only bits named and correct are the 4 signs.

### Deep-learning attack (common/dl_recover.py)

Templates and LDA assume a leakage model (HW). A neural network assumes none, so if
anything beyond HW is present it should find it. MLP, 14,000 train / 2,000 validation
/ 4,000 test, the epoch chosen on the validation split (choosing a model on the test
set introduces optimistic bias), with a shuffled-label control.

| Target | Templates | **Network** |
|---|---:|---:|
| HW information (sum over 18 limbs) | 10.35 bits | **28.75 bits** |
| A single bit (one 32-bit limb) | — | 0.044 bits |
| **Leading T bits (T = 1…8)** | 0 | **0** |

**HW extraction improves 2.8x** (16% → 45% of the 64-bit ceiling), pulling 1.5–2.2
bits each from 15 of the 18 limbs. Combining many samples non-linearly without a model
assumption beats Gaussian templates.

**MSBs are still 0 with a network.** All four coordinates, every T = 1…8, give z ≈ 0
or below and I = 0.000.

### Why the MSBs do not come out

HW is recovered but **does not convert into MSB bits**. Two things block it at once.

1. **The magnitude bit-length is not recoverable** (z ≈ 0). Without knowing which bit
   is the MSB, no high-order bit can be named. `bitlen` is fixed at 401, so the limb
   count is independent of the value, and the magnitude only changes the *values* of
   the high limbs, leaving no structural trace.
2. **HW does not say which bits inside a limb are set.** It is symmetric in the bits,
   so knowing the high limb's HW does not determine the arrangement within it.

This is consistent with the earlier finding that the leakage function is
`HW(32-bit word)` and nothing else.

### Can this feed an MSB-based reconstruction algorithm?

**Not on its own.** With 0 high-order bits there is nothing to hand an algorithm that
takes MSBs as input. The **4 sign bits do come out with 100% certainty**, so if the
signs can be used as a constraint, that is usable.

### Running it

```bash
./conj_vector_host alphas 20000 run1 > /tmp/alphas.txt
make clean PLATFORM=CW308_STM32F4 && make PLATFORM=CW308_STM32F4 CONJ_TRIG_MODE=2
~/.cwvenv/bin/python3 part1_alg318_line2/capture_line2_pico.py --alphas /tmp/alphas.txt -n 20000
~/.cwvenv/bin/python3 common/msb_recover.py   --npz /tmp/sqisign_line2.npz
~/.cwvenv/bin/python3 common/recover_conj_template.py --npz /tmp/sqisign_line2.npz --no-bits
```

## A second target: Sign Algorithm 3.23 line 12 (recovery succeeds)

After Alg 3.18 line 2 proved unrecoverable, this target was chosen by looking for a
point where a known value meets the secret. Alg 3.18 line 5
(`IdealMultiplication`) was considered too, but `I'_sk,chl` and `I_com` are **both
secret** (recovering `E_com` from a signature still does not yield the ideals), so it
hit the same wall.

```c
// sign.c :: compute_challenge_ideal_signature
ibz_copy_digits(&vec.v[1], sig->chall_coeff, NWORDS_ORDER);  // in the signature = public
ibz_mat_2x2_eval(&vec, &sk->mat_BAcan_to_BA0_two, &vec);     // secret matrix x public
  └─ dim2.c:50  ibz_mul(&prod, &mat->m[0][1], &vec->v[1]);   // <- the target
       └─ mp.c:455  prod_ij = (ddigit_t)a.limbs[i] * b.limbs[j];
                    res.limbs[0] = low32(a0 * b0)
```

| | Alg 3.18 line 2 | **Sign line 12** |
|---|---|---|
| Secret | redrawn per signature | **fixed** (part of the private key) |
| Known input | none | `chall_coeff` |
| Leakage function | `HW(limb)` only | `HW(low32(a₀·b₀))` |
| TVLA max\|t\| | 9.92 | **22.02** (control 4.03) |
| Cost per trace | 13.5 s | **65 ms** |
| Result | 64-bit ceiling on 610 bits | **a 32-bit limb recovered** |

### Measurements

- Board result matches the host (`c1/c2 match`); the target multiply is 1560 cycles =
  211 µs
- Sweeping all 16 limb combinations with the correct model, only `(i=0, j=0)` shows
  `|corr| = 0.38` at sample 9904; the other 15 sit at the noise floor. Removing
  `HW(b0)` by partial correlation leaves it at 0.3786, so it is leakage from the
  product itself and not a transfer of the public value's leakage
- For reference, `HW(b0)` alone leaks more strongly at 0.43, but carries no secret

### The attack: two-round candidate-list CPA

Plain CPA fails. It took three stages.

1. **One-byte hypotheses (256) — fails.** Predicting only 8 of 32 bits loses to ghost
   peaks. Going from 100 to 500 traces **does not improve the rank**, which stays at
   3rd–5th. The limit is the model, not the statistics.
2. **16-bit hypotheses (65536) — never reaches first place.** Even at 10,000 traces the
   rank is **fixed at 4**. Projecting out the common `HW(b0)` component changes
   nothing.
3. **Candidate list plus a full 32-bit distinguisher — succeeds.** The usable fact is
   that the correct answer, while not first, is **always within the top 4**. Keeping
   the top K from round one, completing each with its high 16 bits and scoring with
   the full 32-bit model, floats the correct answer to first place.

| Traces | Round-1 rank of the truth | Round-2 outcome |
|---:|---:|---|
| 250 | 2 / 65536 | fails |
| **500** | 4 / 65536 | **recovered** |
| 1000 / 2000 / 4000 | 4 / 65536 | recovered |

K = 16 and K = 256 both succeed, so 16 candidates are enough.

### What remains

What was recovered is **limb 0 (the low 32 bits)** of `Msk[0][1]`. The remaining limbs
follow from `res[1] = low32(a₀b₁ + a₁b₀ + carry)` by repeating the same procedure with
`a₀` known, but that is not implemented yet. `Msk[1][1]` takes the same approach under
trigger mode 13.

> Note: subsequent work extended this to the full matrix and reached 320/512 bits, but
> reproduced on only 2 of 3 disjoint trace ranges, and the HW model cannot separate the
> ×2 and ±1 near-degeneracies. Treat this section as the reliable single-limb result
> and the full-matrix extension as unfinished.

### Usage

```bash
make clean PLATFORM=CW308_STM32F4
make PLATFORM=CW308_STM32F4 CONJ_TRIG_MODE=12
~/.cwvenv/bin/python3 extra_alg323_msk/capture_chall_pico.py --mode tvla   --traces 1000 --skip-flash
~/.cwvenv/bin/python3 extra_alg323_msk/capture_chall_pico.py --mode random --traces 2000 --skip-flash \
    --timebase 3 --samples 30000 --out /tmp/sqisign_chall_cpa.npz
~/.cwvenv/bin/python3 extra_alg323_msk/cpa_chall.py --npz /tmp/sqisign_chall_cpa.npz --find-poi --curve 500,1000,2000
```

## Scope settings worth watching

- **Ch A range**: the CW308 shunt output is around ±30 mV, so the ±500 mV default uses
  only 5% of the ADC and the waveform looks like pure noise. `--cha-range 3` (±100 mV)
  is the default here.
- **Ch B range**: the trigger GPIO is 3.3 V, so a ±2 V range overflows.
  `--chb-range 8` (±5 V) is the default. Only the edge is used so it is not fatal, but
  the overflow flag stays set.
- **Target clock**: CW-Lite `default_setup()` produces 192 MHz / 26 = **7.3846 MHz**,
  not the commonly quoted 7.3728 MHz (confirmed against the captured spectrum peak).
  The scripts read the measured value from the scope with
  `rig.set_target_clock(scope)`.
- **No auto-trigger**: line 1 runs first for several seconds, so `auto_trigger_ms=0`
  waits indefinitely on the hardware trigger. Keeping a rig's usual 5-second
  auto-trigger captures the line 1 region instead.
- **`--skip-flash` needs an nRST**: without reprogramming, the previous session's
  firmware is left mid-protocol and `scope.default_setup()` reinitialises the UART
  underneath it, so the first command breaks. `rig.reset_target(scope)` pulses nRST.

## Throughput

One trace = one run of Alg 3.18 = **about 12–14.4 s** (including vector generation,
transfer and reading back the result). TVLA at 120 traces ≈ 29 min, SNR at 200 traces
≈ 48 min.

## Hardware wiring

- ChipWhisperer-Lite: programs the STM32F415, supplies the 7.3846 MHz target clock,
  SimpleSerial, trigger GPIO
- PicoScope 3403D Ch A: CW308 SMA power measurement (AC coupled — a DC bias saturates
  the input)
- PicoScope 3403D Ch B: trigger GPIO (DC coupled, hardware trigger source)

The default timebase 2 gives ~4 ns/sample on this instrument with 2 channels, about 34
samples per clock at 7.3846 MHz, which is enough to separate limbs.

## Vendoring provenance

See `vendor/sqisign/PROVENANCE.txt`. It is the `src/pqm4/sqisign_p324_3/ref` tree
generated by `scripts/gen_pqm4_sources.sh ref` from upstream `the-sqisign.git` commit
`f417ebd` (2026-09-21), plus `fips202.{c,h}` and `randombytes.h`. The level names
changed: `lvl1 → p324_3`, `lvl3 → p500_27`, `lvl5 → p664_17`.

## Known deviations

- Host vector generation builds inputs at the quaternion level rather than running a
  real keygen/commit. The challenge ideal uses the same call as `sign.c`
  (`id2iso_kernel_dlogs_to_ideal_even`), but `vec` is sampled directly instead of being
  derived from the private key matrix — the distribution has the same shape and no
  keygen is needed.
- A merge-conflict marker is left in upstream around `quaternion.h:436`. It sits inside
  a doc comment and does not affect compilation.
