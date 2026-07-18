# SQIsign Response Core Flow and Measurement Notes

This note summarizes the part currently ported to the CW308 STM32F405 target:
the response quaternion generation corresponding to the spec step
`alpha_rsp <- RandomEquivalentQuaternion(...)`.

## Ported Scope

The board code does not run the full signing protocol. It runs the response core:

```text
inputs:
  secret_ideal     = Isk
  challenge_ideal  = Ichl
  commit_ideal     = Icom

output:
  lattice_content
  resp.denom
  resp.coord[0]
  resp.coord[1]
  resp.coord[2]
  resp.coord[3]
```

The spec-level response quaternion is:

```text
alpha_rsp = (coord0 + coord1*i + coord2*j + coord3*k) / denom
```

`lattice_content` is not part of `alpha_rsp`; it is a bound/verification helper.

## Input Format

The PC sends three ideals as ASCII hex integers:

```text
secret_ideal
challenge_ideal
commit_ideal
```

Each ideal has 18 `ibz_t` values:

```text
norm
lattice.denom
lattice.basis[0][0]
...
lattice.basis[3][3]
```

Total input size is 54 hex integer tokens.

## Output Structure

The output type is `quat_alg_elem_t`, defined in:

```text
src/quaternion/ref/generic/include/quaternion.h
```

Definition:

```c
typedef struct quat_alg_elem
{
    ibz_t denom;
    ibz_vec_4_t coord;
} quat_alg_elem_t;
```

`ibz_vec_4_t` is:

```c
typedef ibz_t ibz_vec_4_t[4];
```

`ibz_t` is:

```c
typedef mpz_t ibz_t;
```

On the STM32F4 build:

```text
GMP_LIMB_BITS = 32
sizeof(ibz_t header) = 12 bytes
sizeof(quat_alg_elem_t header) = 60 bytes
```

The actual big integer limbs are separately allocated through `_mp_d`.

## Response Flow

### 1. Build Challenge-Secret Ideal

Code:

```c
quat_lideal_inter(&challenge_secret,
                  challenge_ideal,
                  secret_ideal,
                  &QUATALG_PINFTY);
```

Meaning:

```text
challenge_secret corresponds to the Isk * Ichl side.
```

This is where `secret_ideal` is heavily used.

### 2. Conjugate Commit Lattice

Code:

```c
quat_lattice_conjugate_without_hnf(&commit_conj,
                                   &commit_ideal->lattice);
```

Meaning:

```text
commit_conj is the conjugate/dual-side lattice for Icom.
```

### 3. Build Response Lattice

Code:

```c
quat_lattice_intersect(&hom_challenge_to_commit,
                       &challenge_secret.lattice,
                       &commit_conj);
```

Meaning:

```text
L = Icom intersection (Isk * Ichl)
```

The response quaternion is sampled from this lattice.

### 4. Compute Lattice Content

Code:

```c
ibz_mul(lattice_content,
        &challenge_secret.norm,
        &commit_ideal->norm);
```

Meaning:

```text
lattice_content = norm(challenge_secret) * norm(commit_ideal)
```

### 5. Set Sampling Bound

Code:

```c
bound = (2^SQIsign_response_length - 1) * lattice_content;
```

For the current board build:

```text
SQISIGN_LEVEL = lvl1
SQIsign_response_length = 126
bound = (2^126 - 1) * lattice_content
```

Other levels:

```text
lvl1: 126
lvl3: 192
lvl5: 253
```

### 6. Rejection Sampling

Code:

```c
quat_lattice_sample_from_ball(resp_quat,
                              &hom_challenge_to_commit,
                              &QUATALG_PINFTY,
                              &bound);
```

Internal flow:

```text
1. Compute Gram matrix G of the lattice.
2. Correct radius by lattice denominator.
3. Compute a bounding parallelogram box.
4. Sample random candidate x[0..3] in the box.
5. Apply U transform to x.
6. Evaluate quadratic form tmp = qf_G(x).
7. Reject if tmp == 0 or tmp > radius.
8. Repeat until accepted.
```

Rejection count is not fixed, so absolute timing from command start to finalization is not stable.

### 7. Convert Accepted Sample to Quaternion Coordinates

After an accepted `x`:

```c
ibz_mat_4x4_eval(&res->coord, &lattice->basis, &x);
ibz_copy(&res->denom, &lattice->denom);
```

Meaning:

```text
coord = lattice.basis * accepted_x
denom = lattice.denom
```

This is where the response coordinate candidate is first materialized.

### 8. Normalize

Normalization:

```text
gcd = gcd(content(coord[0..3]), denom)

denom  /= gcd
coord0 /= gcd
coord1 /= gcd
coord2 /= gcd
coord3 /= gcd

if denom < 0:
  denom *= -1
  coord[i] *= -1
```

This gives the normalized representation of:

```text
alpha_rsp = coord / denom
```

## Observed Value Ranges

PC sampling of 1000 protocol executions, lvl1:

```text
lattice_content : 516..528 bits
denom           : 1..2 bits
coord0          : 315..328 bits
coord1          : 315..327 bits
coord2          : 188..202 bits
coord3          : 188..201 bits
```

KAT first case:

```text
lattice_content : 521 bits
denom           : 1 bit
coord0          : 322 bits
coord1          : 323 bits
coord2          : 197 bits
coord3          : 197 bits
```

Approximate storage on STM32F4, using 32-bit limbs:

```text
denom  : 1 limb
coord0 : 10..11 limbs
coord1 : 10..11 limbs
coord2 : 6..7 limbs
coord3 : 6..7 limbs
```

## Cycle Measurements

Target clock:

```text
DF_CPU = 7,372,800 Hz
ADC clock = 29,454,545 Hz
1 target cycle ~= 4 ADC samples
```

Measured KAT profile:

```text
handler total:       88,648,839 cycles ~= 12.025 s
parse input:             91,450 cycles ~= 0.012 s
core compute:        88,526,640 cycles ~= 12.008 s
serialize output:        20,159 cycles ~= 0.003 s
```

Core breakdown:

```text
quat_lideal_inter:          19,835,239 cycles ~= 2.690 s
commit conjugate:                3,591 cycles ~= 0.0005 s
lattice_intersect:          30,338,231 cycles ~= 4.115 s
content mul:                     2,082 cycles ~= 0.0003 s

sample_total:               38,319,601 cycles ~= 5.198 s
  sample_setup:             37,854,139 cycles ~= 5.135 s
  rejection loop:              317,472 cycles ~= 0.043 s
  sample_finalize:             142,701 cycles ~= 0.019 s
    normalize_total:           107,696 cycles ~= 0.0146 s
      normalize_gcd:            96,241 cycles ~= 0.0131 s
      coord0_div:                1,841 cycles ~= 0.00025 s
      coord1..3_div:             5,346 cycles ~= 0.00073 s
      sign mul:                  3,493 cycles ~= 0.00047 s
```

The example run had:

```text
rejection_iters = 2
```

Absolute offsets from command start depend on rejection count. Trigger-relative offsets are more meaningful.

## Trigger-Relative Regions

`RESP_TRIGGER_MODE=2` starts at `sample_finalize`:

```text
mode2 sample 0 .. 140,020:
  basis*x + denom copy

mode2 sample 140,020 .. 570,804:
  normalize_total

mode2 sample 524,984 .. 532,348:
  coord0_div

mode2 sample 532,348 .. 553,732:
  coord1..3_div

mode2 sample 553,732 .. 567,704:
  sign multiplication
```

`RESP_TRIGGER_MODE=3` starts at `coord0_div`:

```text
mode3 sample 0 .. 7,364:
  coord0_div
```

## Preliminary TVLA Results

Welch fixed-vs-random t-test, 20 fixed traces and 20 random traces per region.

Fixed group:

```text
same vector, same RNG seed each trace
```

Random group:

```text
same vector, different RNG seed per trace
```

Results:

```text
mode2 offset 0:
  region: sample_finalize start, first 20k ADC samples
  max |t| = 1.624
  |t| > 4.5 hits = 0

mode2 offset 520k:
  region: coord0_div vicinity inside mode2 trigger
  max |t| = 1.333
  |t| > 4.5 hits = 0

mode3 coord0_div:
  region: direct coord0_div trigger
  max |t| = 2.095
  |t| > 4.5 hits = 0
```

Interpretation:

```text
No first-order fixed-vs-random TVLA signal was observed in these three tested output-finalization regions at 20+20 traces.
```

This does not prove absence of leakage. It means the currently tested regions and sample size do not show a strong first-order mean difference.

## Useful Artifacts

Cycle/profile-enabled firmware commands:

```text
p  profile text into output buffer
l  output length
g  output chunk
```

TVLA trace files:

```text
/tmp/tvla_mode2_offset0_fixed20.npz
/tmp/tvla_mode2_offset0_random20.npz
/tmp/tvla_mode2_offset520k_fixed20.npz
/tmp/tvla_mode2_offset520k_random20.npz
/tmp/tvla_mode3_coord0div_fixed20.npz
/tmp/tvla_mode3_coord0div_random20.npz
```

Single-position sanity traces:

```text
/tmp/sqisign_mode2_offset0.npz
/tmp/sqisign_mode2_offset520k.npz
/tmp/sqisign_mode3_coord0div_offset0.npz
```
