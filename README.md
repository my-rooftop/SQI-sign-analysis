# Power analysis on SQIsign quaternion arithmetic

Side-channel measurement and recovery experiments against the quaternion arithmetic
of **SQIsign Round 3** — upstream [`the-sqisign.git`](https://github.com/SQISign/the-sqisign)
commit `f417ebd` (2026-09-21) against the 2026-09-01 specification
(`sqisign-20260901.pdf`), parameter set `p324_3`.

Everything here is measured on real hardware: STM32F415 (Cortex-M4) on a CW308
carrier, clocked at 192/26 = 7.3846 MHz, captured with a PicoScope 3403D at
4 ns/sample.

---

## Results

Two different operations inside **Algorithm 3.18 `SampleQuaternionResponse`** were
measured — not two algorithms, two lines of the same one.

| | line 2 `Conjugate` | line 5 `RandomBoundedElement` |
|---|---|---|
| Secret to infer | 614 bits (α_sk,chl) | 746 bits (α_rsp) |
| Operation | in-place negation | multiply-accumulate formation |
| Leakage detected (TVLA) | yes, max\|t\| 9.92 | yes, max\|t\| 25.62 |
| …at matched trace counts | 9.92 | 10.49 ± 0.84 — **indistinguishable** |
| Hamming-weight information | +9.9 bits | +8.9 bits |
| **Signs** | **4/4 bits, single trace** | **4/4 bits, 99.75–99.95%** |
| Magnitude bit-length | not recoverable | not recoverable |
| **High-order bits (MSBs)** | **0 bits** | **0 bits** |
| Report | [part 1](cw-sqisign-conj/part1_alg318_line2/LINE2_EXPERIMENT.md) | [part 2](cw-sqisign-conj/part2_alg318_line5/RESPONSE_EXPERIMENT.md) |

**Both operations leak, and neither gives up the high-order bits that an
MSB-based reconstruction needs.** What an attacker secures is the 4 sign bits per
quaternion; the Hamming-weight information only shrinks the search space (614 → 604
and 746 → 737 bits respectively).

The two targets differ in size, in operation structure, and in what fraction of their
limbs is sign extension — and reach the same conclusion. That points at a limit of the
**Hamming-weight leakage model itself** rather than at either implementation.

The sharpest single piece of evidence is the deep-learning result in part 2: the same
network on the same traces learns the Hamming weight well (13.7% → 48%, z = +64) and
the full byte not at all (0.71% → 0.63%, 0 of 3 limbs significant). Naming a byte
requires separating values *inside* one Hamming-weight class, so a model that learns
the weight and not the byte has found the weight and nothing else.

---

## Layout

```
cw-sqisign-conj/              the measurement package (Round 3, p324_3)
├── Makefile, *.c, *.h        firmware -- one binary; CONJ_TRIG_MODE picks the target
├── host/, vendor/            host tool and vendored upstream sources
├── common/                   rig and analysers shared by both parts
├── part1_alg318_line2/       section 5.3.1 -- capture, report, figures
├── part2_alg318_line5/       section 5.3.2 -- capture, driver, report, figures
└── extra_alg323_msk/         Algorithm 3.23 line 12 (see below) -- unfinished

sqisign-20260901.pdf          the specification these experiments target
```

Start with [`cw-sqisign-conj/README.md`](cw-sqisign-conj/README.md) for the
instrumentation, the trigger modes, and how to reproduce a campaign.

### A third target, not part of the above

[`extra_alg323_msk/`](cw-sqisign-conj/extra_alg323_msk/) attacks `Msk` at
**Algorithm 3.23 `Sign` line 12**, where a fixed secret meets a known input — so
unlike the two above, CPA applies there. A single 32-bit limb of `Msk[0][1]` recovers
reliably from 500 traces. The full-matrix extension reaches 320/512 bits but
reproduces on only 2 of 3 disjoint trace ranges, and the Hamming-weight model cannot
separate the ×2 and ±1 near-degeneracies. **It is unfinished and should not be cited
as a result.**

---

## Notes on method

A few practices that the measurements here depend on, recorded because getting them
wrong produced believable but false results during this work:

- **Compare TVLA at matched trace counts.** t grows as √N, so 25.62 on 1,000 traces
  and 9.92 on 120 are the same leakage, not a 2.5× difference.
- **Judge recovery on information, not accuracy.** 54.5% per-bit accuracy beats its
  baseline significantly and is worth 0.006 bits.
- **Run a shuffled-label control beside any information estimate.** Choosing points of
  interest by correlation and fitting templates at those same points extracts bits
  from pure noise — 5 bits of them, on a small enough set.
- **Watch for the sign in disguise.** A boundary limb has thousands of distinct values
  yet a Hamming weight decided almost entirely by the sign, and will otherwise top the
  ranking as the strongest "value" leak.
- **Use the max-over-samples threshold.** `4/√N` is the bound for one pre-chosen
  sample; against a maximum over tens of thousands, use `√(2 ln S)/√N`.

Both experiment reports carry a section on the traps encountered and the rig settings
that actually changed results.

---

## Reproducing

Requires an STM32F415 on CW308, a ChipWhisperer-Lite, a PicoScope 3403D, and the
`arm-none-eabi` toolchain.

```bash
cd cw-sqisign-conj
./build_host_tool.sh
make PLATFORM=CW308_STM32F4 CONJ_TRIG_MODE=2      # or 18 for line 5
~/.cwvenv/bin/python3 common/check_trigger.py --traces 3
```

Full per-part instructions are in each part's report.
