# SQIsign Algorithm 3.18 line 2 — `Conjugate(α_sk,chl)` 부채널 실험

2026-09-01 spec(`sqisign-20260901.pdf`)의 **Algorithm 3.18 SampleQuaternionResponse, line 2**를
타겟으로 하는 측정 패키지다.

```
Algorithm 3.18 SampleQuaternionResponse                                    (spec p.21)
 1: I'_sk,chl, α_sk,chl ← EquivalentCoprimeIdeal(I_sk,chl, 0, BitLenSkChall, prng_def)
 2: α_sk,chl ← Conjugate(α_sk,chl)                                    ← 측정 대상
 3: if GCD(IdealNorm(I_com), IdealNorm(I'_sk,chl)) ≠ 1 then raise
 5: q_rsp, α_rsp ← RandomBoundedElement(IdealMultiplication(I'_sk,chl, I_com), 2^(e_rsp−1), prng_res)
```

구 spec 기반 실험은 `../cw-sqisign-response/`에 그대로 남아 있다. 이 패키지는 신규 구현을
별도로 vendoring하므로 둘은 독립적이다.

## spec → 코드 대응

| spec | vendored 코드 |
|---|---|
| Alg 3.18 전체 | `protocol.c` `quat_response_element()` |
| Alg 3.18 line 1 (`EquivalentCoprimeIdeal`) | `ideal.c:846` `quat_ideal_small_equivalent_coprime()` |
| **Alg 3.18 line 2 (`Conjugate`)** | **`protocol.c:37` `quat_alg_conj(sk_chall_quat, sk_chall_quat)`** |
| Alg 4.38 `Conjugate` | `algebra.c:168` `ibz_copy`×2 (denom, coord0) + `ibz_neg`×3 (coord1,2,3) |
| 실제 산술 | `mp.c:469` `ibz_neg()` — limb 단위 2의 보수 부정 루프 |
| Alg 4.74 `IdealMultiplication` | `internal.h:734` `quat_ideal_mul_O0()` |
| Alg 4.83 `RandomBoundedElement` | `internal.h:571` `quat_lattice_sample_from_ball()` |

Alg 3.18 line 3의 GCD 검사는 upstream이 `sign.c:276`으로 끌어내 두었다.

## 왜 측정 가능한가

2026 구현은 mini-gmp를 자체 `mp` 모듈로 교체했고 `ibz_t`가 **고정폭 2의 보수**가 되었다:

```c
struct ibz_t { int bitlen; digit_t limbs[IBZ_NLIMBS]; };   // p324_3, RADIX=32 → 60

void ibz_neg(ibz_t *neg, const ibz_t *a) {                 // mp.c:469
    int n = NUM_LIMBS(a->bitlen);
    digit_t carry = 1;
    for (int i = 0; i < n; i++) {
        ddigit_t t = (ddigit_t)~a->limbs[i] + carry;
        neg->limbs[i] = t;                                 // 모든 limb를 읽고 쓴다
        carry = t >> NUM_BITS_LIMB;
    }
    neg->bitlen = a->bitlen;
}
```

구 mini-gmp에서는 `mpz_neg`가 `mpz_set` + 부호비트 반전이고 `mpz_set`에 `r != x` 조기반환이
있어서(`mini-gmp.c:1533`) **in-place conjugate가 limb을 전혀 건드리지 않았다**. 신규 구현은
in-place여도 전체 limb에 대해 `LDR → MVN → ADCS → STR`가 돌기 때문에 limb별 HW/HD 누설이
그대로 나온다. 또한 루프 횟수 `n = NUM_LIMBS(bitlen)`이 값이 아니라 **선언된 경계값**에만
의존하므로 trace마다 슬롯 위치가 동일하다 — 분기도 없다.

## 복구 대상 크기 (호스트 실측, `conj_vector_host selftest`)

| | coord[1] (`i`) | coord[2] (`j`) | coord[3] (`ij`) |
|---|---|---|---|
| 선언 폭 `bitlen` | 401 (고정) | 401 (고정) | 401 (고정) |
| `ibz_neg` 순회 limb | 13 | 13 | 13 |
| 실제 값 크기 | 230–235 bit | 69–73 bit | 69–73 bit |
| sign-extension limb | 5 | 10 | 10 |
| **정보 있는 limb** | **8** | **3** | **3** |

- `ibz_neg`가 도는 총 limb: **39개 = 1248 bit 폭**, 전 trace 동일
- 그중 **정보를 담은 limb은 14개 = 448 bit**
- **실제 복구 대상 ≈ 376 bit** (234 + 71 + 71)

coord[2]/[3]가 작은 이유: norm form이 `a² + b² + p(c² + d²)`이고 log₂p = 324이므로
c, d는 a, b보다 약 √p배 작다 (234 − 71 = 163 ≈ 324/2).

**coord[0] (`a`, ~234 bit)와 denom은 negate가 아니라 `ibz_copy`(= `memmove`)** 다.
in-place `memmove`에서 실제로 데이터가 움직이는지는 libc 구현에 달려 있어 trace로 확인해야
한다 — `CONJ_TRIG_MODE=4`가 이 두 구간에 별도 트리거를 넣는다.

## 구성

```
vendor/sqisign/p324_3/       upstream 임베디드 소스 (the-sqisign.git f417ebd,
                             scripts/gen_pqm4_sources.sh ref 산출물 + fips202)
vendor/sqisign/host-helpers/ 호스트 전용 upstream 테스트 헬퍼 (PROVENANCE 참고)
conj_core.{c,h}              계측된 Algorithm 3.18. upstream quat_response_element()를
                             line 단위로 그대로 옮기고 트리거/프로파일만 추가
simpleserial-sqisign-conj.c  CW308_STM32F4 / STM32F415 simpleserial 타겟
cw_randombytes.c             보드 RNG shim
host/conj_vector_host.c      벡터 생성 + KAT + ground truth 출력 (PC)
build_host_tool.sh           호스트 도구 빌드
Makefile                     펌웨어 빌드
conj_rig.py                  PicoScope 3403D + CW-Lite 공용 리그
check_trigger.py             (3) 트리거 정렬 검증
capture_conj_pico.py         (2) 캡처 — random / fixed / tvla
analyze_conj_leakage.py      (4) TVLA + SNR + correlation
```

펌웨어는 quaternion + mp + PRNG만 링크한다 (ROM 64 KB / 1 MB, 정적 RAM 20.8 KB / 128 KB).
EC·GF·isogeny 코드는 들어가지 않는다 — challenge ideal은 호스트에서 만들어 데이터로 보낸다.

## 트리거 모드

빌드 시 `CONJ_TRIG_MODE=n`으로 고른다.

| Mode | 트리거 구간 | 펄스 수 |
|---:|---|---:|
| 0 | Algorithm 3.18 전체 | 1 |
| 1 | line 1 `EquivalentCoprimeIdeal` | 1 |
| **2** | **line 2 `Conjugate` 전체 (기본값)** | **1** |
| 3 | line 2, `ibz_neg` coord[1]/[2]/[3] 각각 | 3 |
| 4 | line 2, denom copy / coord[0] copy / coord[1..3] neg | 5 |
| 5 | line 5 `IdealMultiplication` | 1 |
| 6 | line 5 `RandomBoundedElement` | 1 |
| 7 | α_sk,chl limb 순회 인공 probe (`CONJ_PROBE_REPS`) | 3 |
| 8 | line 2, coord[1]의 **limb마다** 1펄스 | 13 |

Mode 8은 `ibz_neg` 루프 반복마다 펄스를 주므로 캡처가 limb 인덱스로 분절되어 정렬
불확실성이 사라진다. Mode 7은 실제 코드 경로가 아니라 누설 상한/프로파일링 기준선용이다.

`CONJ_OUT_OF_PLACE=1`은 conjugate 대상과 출력을 분리한다. 기본값 0이 spec·upstream 그대로의
in-place이며, 위에서 설명한 이유로 in-place도 충분히 누설한다.

## 실행 순서

### 0. 호스트에서 크기 확정 + KAT
```bash
./build_host_tool.sh
./conj_vector_host selftest 8          # 불변식 검증 + limb/bit 수 실측
./conj_vector_host gen myseed > vec.txt
./conj_vector_host run vec.txt         # α_sk,chl 전/후 값과 raw limb 출력
```

`run`은 다음을 검증한다: `Norm(α_rsp) = q_rsp · n(I'_sk,chl) · n(I_com)`,
`bitsize(q_rsp) ≤ e_rsp − 1`, `Conjugate`가 involution인지, line 3의 GCD 조건.

### 1. 펌웨어 빌드
```bash
make clean PLATFORM=CW308_STM32F4
make PLATFORM=CW308_STM32F4 CONJ_TRIG_MODE=2
```
컴파일 옵션을 바꿀 때는 반드시 `make clean`을 먼저 한다.

### 2·3. 트리거 정렬 검증
```bash
~/.cwvenv/bin/python3 check_trigger.py --traces 3
```
확인 항목:
- 보드의 Conjugate가 산술적으로 맞는지 (coord[1..3]만 부호 반전, coord[0]·denom 불변)
- 트리거 펄스 수가 모드별 기대값과 일치하는지
- 트리거 high 폭을 target cycle로 환산한 값이 펌웨어가 보고한 **DWT cycle 카운트**와 맞는지
- limb당 cycle 수가 `ibz_neg` 한 반복 비용과 일관되는지
- 트리거 전압 스윙이 충분한지

이 단계가 통과하지 않으면 뒤 단계는 의미가 없다. Alg 3.18 line 1이 먼저 수 초간 돌기
때문에 PicoScope는 **auto-trigger 없이**(`auto_trigger_ms=0`) 하드웨어 트리거를 무한 대기한다.

### 4. 캡처
```bash
# SNR/프로파일링용: 매 trace마다 새 입력 + 새 시드
~/.cwvenv/bin/python3 capture_conj_pico.py --mode random --traces 200 --skip-flash

# TVLA용: fixed / random을 trace 단위로 교차
~/.cwvenv/bin/python3 capture_conj_pico.py --mode tvla --traces 200 --skip-flash
```
`--crop`을 주면 트리거 구간만 저장한다. 각 trace는 Alg 3.18 1회 실행 비용이 들고 그 대부분은
line 1이므로, 먼저 `check_trigger.py`로 1회 소요 시간을 재고 시간을 계산할 것.

fixed/random 교차(interleave)는 온도·전원 드리프트가 누설로 오인되는 것을 막는다.

### 5. 누설 정량화
```bash
~/.cwvenv/bin/python3 analyze_conj_leakage.py \
    --tvla /tmp/sqisign_conj_tvla.npz \
    --snr  /tmp/sqisign_conj_random.npz
```

- **TVLA** — fixed vs random Welch t. `|t| > 4.5`면 누설 판정. 같은 출력에
  **fixed vs fixed** 반쪽 검정도 함께 나오는데, 이쪽이 4.5를 넘으면 리그가 드리프트하는
  것이므로 fixed-vs-random 결과를 신뢰할 수 없다.
- **SNR** — limb의 HW를 클래스로 하는 per-limb SNR. HW는 binomial(32, ½)이라 극단
  클래스가 거의 비어 SNR argmax가 불안정하므로, 기본적으로 HW를 8개 bin으로 묶고
  (`--hw-bins`) trace 수가 `--min-class`(기본 5) 미만인 클래스를 버린다.
- **correlation** — `HW(limb)`, `HW(~limb)`에 대한 per-sample Pearson r. argmax가 SNR보다
  안정적이므로 관심 지점(POI)은 `|corr|` 피크로 잡고 SNR은 그 지점에서 읽는다.
  두 통계가 같은 샘플에서 함께 튀는 것이 원하는 증거다. 노이즈 하한은 `4/√N`으로 보고한다.
- sign-extension limb은 `label constant, no SNR`로 배제된다.

판정 기준: `|corr|` 피크가 `4/√N` 하한을 넘고 SNR이 대략 0.1 이상이면 CPA가 수백 trace로
수렴할 수준, 0.01–0.1이면 수천 trace, 그 이하면 복구 시도 전에 측정 자체(프로브·레인지·
평균)를 개선해야 한다.

## 검증 상태

- 호스트 Algorithm 3.18: upstream 불변식 전부 통과 (`selftest: OK`)
- ARM 펌웨어: 경고 없이 빌드 (ROM 64 KB, 정적 RAM 20.8 KB)
- TVLA/SNR/correlation 수치 경로: 알려진 HW 누설을 심은 합성 데이터로 검증.
  correlation locator가 주입 8개 limb 중 7개의 시간 위치를 정확히 복원했고,
  sign-extension limb 25개는 모두 올바르게 배제됨
- **실측 trace는 아직 없음** — 2~5단계는 하드웨어(PicoScope 3403D + CW-Lite + CW308/STM32F415)
  연결 후 실행해야 한다

## 하드웨어 배선

`../cw-sqisign-response/`의 기존 리그와 동일하다.

- ChipWhisperer-Lite: STM32F415 프로그래밍, 7.3728 MHz 타겟 클럭, SimpleSerial, 트리거 GPIO
- PicoScope 3403D Ch A: CW308 SMA 전력 측정 (AC 결합 — DC 바이어스가 입력을 포화시킨다)
- PicoScope 3403D Ch B: 트리거 GPIO (DC 결합, 하드웨어 트리거 소스)

기본 timebase 2는 이 장비에서 2채널 기준 ~4 ns/sample이고, 7.3728 MHz 클럭당 약 34 샘플이라
limb 단위 분리에 충분하다.

## vendoring 출처

`vendor/sqisign/PROVENANCE.txt` 참고. upstream `the-sqisign.git` 커밋 `f417ebd`
(2026-09-21)에서 `scripts/gen_pqm4_sources.sh ref`로 생성한 `src/pqm4/sqisign_p324_3/ref`
트리에 `fips202.{c,h}`와 `randombytes.h`를 더한 것이다. 레벨 이름이 바뀌었다:
`lvl1 → p324_3`, `lvl3 → p500_27`, `lvl5 → p664_17`.

## 알려진 편차

- 호스트 벡터 생성은 실제 keygen/commit 대신 quaternion 수준에서 입력을 만든다. challenge
  ideal은 `sign.c`와 같은 호출(`id2iso_kernel_dlogs_to_ideal_even`)을 쓰지만, `vec`는 비밀키
  행렬에서 유도하는 대신 직접 샘플링한다 — 분포 형태는 같고 keygen이 필요 없다.
- upstream `quaternion.h:436` 부근에 병합 충돌 마커가 남아 있다. doc 주석 안이라 컴파일에는
  영향이 없다.
