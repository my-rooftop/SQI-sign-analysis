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
보드 실측 결과 in-place `ibz_copy`가 **144 cycle (limb당 11 cycle)** 을 소비한다 — no-op이
아니므로 negate되지 않는 `a` 좌표도 버스에 실리고, 누설 후보다. `CONJ_TRIG_MODE=4`가 이
두 구간에 별도 트리거를 넣는다.

### 벡터 파일의 `bounds` 줄

`ibz_t`는 선언된 비트 경계값을 갖고 하위 연산이 이를 상속하는데, `ibz_set_from_str()`
(`mp.c:311`)은 파싱한 크기에서 경계값을 다시 잡는다. 이 차이가 α_sk,chl까지 전파되어
line 2에서 `ibz_neg`가 도는 limb 수를 바꾼다 (hex 왕복 시 401/13 → 319/10). 그래서 벡터
파일이 14개 입력값의 경계값을 `bounds` 줄로 함께 싣고, 호스트와 펌웨어가 파싱 후
`ibz_set_bound()`로 복원한다. 이게 없으면 보드는 실제 서명 장치보다 좁은 구간을 실행한다.
펌웨어는 `p` 명령의 `bounds_applied`로 적용 여부를 보고한다.

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

## 보드 실측 (STM32F415 @ 7.3846 MHz)

| 구간 | cycles | 비중 |
|---|---:|---:|
| Algorithm 3.18 전체 | 73–92 M (10–12.5 s) | 100% |
| line 1 `EquivalentCoprimeIdeal` | 43–50 M | ~53% |
| line 5 `RandomBoundedElement` | 41–42 M | ~46% |
| line 5 `IdealMultiplication` | 0.19 M | 0.2% |
| `Intersect` (line 1 이전) | 0.25 M | 0.3% |
| **line 2 `Conjugate`** | **1,121 (152 µs)** | **0.0013%** |

line 2 내부: `coord[0]` copy 144 cycle, `coord[1]/[2]/[3]` negate 각 247 cycle
→ **limb당 19.0 cycle** × 13 limb.

트리거(mode 2)는 pre-sample 위치에 정확히 붙고 폭이 **trace 간 sd = 0.0** 이다
(34,739 샘플 = 1024 cycle @ 4 ns). DWT가 보고한 1121 cycle 대비 비율 0.91인데, 차이는
트리거 밖에 있지만 계측 구간 안에 있는 bitlen/limb 프로파일링 코드다. 정렬 작업이 필요 없다.

## 검증 상태

- 호스트 Algorithm 3.18: upstream 불변식 전부 통과 (`selftest: OK`)
- ARM 펌웨어: 경고 없이 빌드 (ROM 64 KB, 정적 RAM 20.8 KB)
- **보드**: `check_trigger.py` 통과 — conjugate 산술 검증, 펄스 수, DWT 대비 폭,
  limb당 cycle, bitlen 401/13 limb 재현 (`bounds_applied=1`)
- **신호**: Ch A ±100 mV에서 overflow 0, 스펙트럼에 7.384 MHz 타겟 클럭과
  14.8 / 29.5 / 36.9 / 44.3 MHz 고조파. limb 주기(644 샘플)로 접은 평균에 재현되는
  루프 시그니처가 나온다
- TVLA/SNR/correlation 수치 경로: 알려진 HW 누설을 심은 합성 데이터로 검증.
  correlation locator가 주입 8개 limb 중 7개의 시간 위치를 정확히 복원했고,
  sign-extension limb 25개는 모두 올바르게 배제됨

## 누설 측정 결과 (실측)

### TVLA — fixed 60 / random 60, CONJ_TRIG_MODE=2

```
fixed vs random : max|t| = 9.92 at sample 13757 (cycle ~406)
                  2647 of 34735 samples exceed 4.5 (7.62%)
fixed vs fixed  : max|t| = 4.23, 0 samples over threshold  OK
verdict: LEAKAGE DETECTED
```

대조군(fixed vs fixed)이 임계값을 넘지 않으므로 리그 드리프트가 아니라 데이터 의존
누설이다. 120 trace만으로 검출됐다.

t-trace는 균일하지 않고 **4개 밴드로 뭉친다**. 밴드 중심은 대략 175 / 428 / 635 / 893
cycle이고 간격이 ~250 cycle로, 측정된 좌표당 비용 247 cycle과 일치한다. 네 밴드가 각각
`coord[0]` copy와 `coord[1]`/`[2]`/`[3]` negate에 대응하는 것으로 보이며, 밴드 사이
조용한 구간은 프로파일링 호출과 루프 셋업이다.

`coord[0]` 밴드에도 |t| ≈ 9.9가 뜨는데, 이는 negate하지 않고 `memmove`로 복사만 하는
`a` 좌표도 누설한다는 디스어셈블리 분석과 일치한다.

밴드 ↔ 좌표 대응은 아래 SNR 결과가 실측으로 확정한다.

### per-limb SNR / correlation — random 200 trace, CONJ_TRIG_MODE=2

sign-extension limb 30개를 제외한 **정보 있는 limb 22개 전부가 노이즈 하한 0.283을 넘는다.**

| 좌표 | limb | \|corr(HW)\| | SNR@poi | 위치 (cycle) |
|---|---|---|---|---|
| coord[0] (`a`, `ibz_copy`) | 0–6 | 0.46 – 0.80 | 0.45 – 3.69 | 121 → 143 |
| coord[1] (`b`, `ibz_neg`) | 0–6 | 0.57 – 0.72 | 0.90 – 2.65 | 262 → 358 |
| coord[2] (`c`, `ibz_neg`) | 0–1 | 0.40, 0.49 | 0.86, 0.92 | 524, 540 |
| coord[3] (`d`, `ibz_neg`) | 0–1 | 0.44, 0.50 | 0.64, 0.71 | 784, 802 |

**limb이 시간순으로 행진한다.** coord[1]의 limb 0→6이 샘플 8862 → 12106에 순서대로,
간격 ~540 샘플 ≈ 16 cycle로 찍힌다 (실측 19 cycle/limb과 일치). limb별 시간 분리가
데이터로 확인된다.

좌표별 구간도 이로써 확정된다: coord[0] copy = cycle 121–161, coord[1] neg = 262–438,
coord[2] = 524–652, coord[3] = 784–834. 좌표 시작 간격 262 / 260 cycle이 측정된
247 cycle + 프로파일링 오버헤드와 맞는다.

#### 수치를 읽을 때 주의할 두 가지

**각 좌표의 최상위 limb(\|corr\| 0.94–0.97)은 부호 누설이다.** coord[0] limb 7,
coord[1] limb 7, coord[2] limb 2, coord[3] limb 2가 0.94–0.97로 나오지만, 이들의 피크
샘플(5458 / 14830 / 22065 / 28227)이 sign-extension limb들의 피크와 같다. 값 크기가
limb 경계에 못 미쳐 상위 limb이 대부분 부호 비트이기 때문이다. **진짜 per-limb HW
누설은 하위 limb 18개의 \|corr\| 0.40–0.80** 이다. `analyze_conj_leakage.py`의
`classify_limb()`이 limb을 constant / sign / informative로 분류하고 sign은 순위에서
제외한다 — 이걸 안 하면 sign-ext limb들이 동일한 0.97로 상위를 채운다.

**"limb당 13 trace" 추정치는 낙관적이다.** Mangard 공식을 median r = 0.618에 적용한
값인데, 이는 *정답 가설에서의 분리*를 재는 것이다. 실제로 2^32 가설을 열거하는 공격은
훨씬 더 필요하고, 바이트 단위로 쪼개도 구별자 여유는 줄어든다.

그래도 \|corr\| 0.40–0.80, SNR 0.6–3.7은 일반적인 "활용 가능" 기준(SNR ≈ 0.1)을 한참
웃돈다.

## 복구 실험 (recover_conj_template.py)

### 왜 CPA가 아니라 템플릿 공격인가

α_sk,chl은 서명마다 새로 뽑힌다 (Alg 3.18 line 1). 고정된 비밀을 여러 trace로 평균 내는
구조가 아니므로 CPA가 성립하지 않는다. trace를 많이 모으는 것은 **프로파일링** 용이고,
실제 공격은 해당 서명의 **trace 한 장**에서 복구해야 한다.

이 실험은 같은 보드에서 ground truth를 알고 프로파일링한다 — 공격자가 마주할 가장 쉬운
조건이므로 결과는 상한으로 읽어야 한다.

### 누설 함수는 `HW(32-bit word)` 하나뿐이다

| 누설 모델 | 200 trace | 20,000 trace |
|---|---|---|
| `HW(word)` | **0.40 – 0.80**, 확실 | 확실 |
| 개별 비트 (binary LDA) | 576개 중 0개 초과 | **481/576 초과** (평균 54.5%, 기준선 50.7%) |
| 바이트별 HW (편상관) | 72개 중 3개 경계선 | — |
| 연속 limb 간 HD | 17개 중 2개 경계선 | — |

**비트 수준 누설은 존재한다.** 200 trace에서 0/576이었던 것은 검정력 부족이었다
(공격 trace가 40개뿐). 20,000 trace로 다시 재면 481/576이 기준선을 유의하게 넘고,
신경망도 같은 것을 본다. 즉 "원리적으로 불가능"이 아니라 **"존재하지만 무의미"** 가 맞다.

크기가 그렇다. 정확도 54.5%는 완벽히 보정돼도 비트당 `1 − H(0.545) = 0.0058 bit`이고
576비트 전체로도 최대 3.3 bit다. 신경망 실측은 32비트 limb 하나당 **0.044 bit**,
18 limb 합계 약 0.8 bit다. 정확도만 보면 "복구됐다"고 오판하기 쉬우므로 정보량으로
판정해야 한다.

기전은 `ibz_neg`의 캐리 체인으로 보인다 — limb 밖으로 나가는 캐리는 그 limb이 0인지에
의존하는데, 이는 비트에 대해 대칭인 함수가 아니다. 따라서 HW 위에 작은 비트 수준
의존성이 얹히는 것은 예상 범위다.

바이트별 HW의 원상관이 0.5배로 나오는 것은 `HW(word) = Σ HW(byte)`의 수학적 귀결이지
독립적인 바이트 누설이 아니다. 편상관으로 `HW(word)`를 제거하면 사라진다.

### 정보이론적 상한

```
limb당 추출 가능 최대 = H(HW of 32-bit word) = 3.55 bit
18개 깨끗한 limb × 3.55            =  64 bit   ← 천장
실측 (템플릿, 20,000 trace)         =  10 bit
실측 (신경망, 20,000 trace)         =  29 bit
복구 대상                          = 610 bit
남는 탐색 공간                      ≈ 546 bit
```

**단일 trace로 α_sk,chl 복구는 불가능하다.** 누설 자체는 강하지만(TVLA max|t| = 9.92)
복구에는 쓸 수 없다.

이유는 구조적이다: **`Conjugate`는 비밀을 알려진 값과 섞지 않는다.** AES SubBytes처럼
"알려진 평문 + 추측한 키 바이트 → 중간값 예측"하는 분할정복 구조가 없고 단순 비트 보수다.
그래서 trace가 흘리는 것은 `HW(limb)` 하나이고, limb당 3.55 bit가 천장이다. trace를
늘려도 이 천장은 올라가지 않는다 — HW 추정 정확도만 올라간다.

### 분석 함정 두 가지

**sign-extension limb**: 52개 limb 중 30개는 `0x00000000`/`0xffffffff` 둘 중 하나라
좌표당 부호 1 bit, 총 4 bit만 담는다. 이들은 동일한 레이블을 공유해 같은 샘플에서
|corr| 0.97로 나타나며, 걸러내지 않으면 "limb 5개가 0.98로 샌다"는 과대평가가 된다.
`classify_limb()`이 constant / sign / informative로 분류한다.

**max-over-samples 임계값**: `4/√N`은 **미리 정한 표본 하나**에 대한 임계값이다. 34,734개
샘플에 대한 **최댓값**과 비교하면 순수 노이즈에서도 누설 판정이 난다. 올바른 기준은
`√(2 ln S)/√N`이고, N=200·S=34734에서 0.283이 아니라 **0.323**이다. 분석기는 두 값을
모두 출력하고 후자로 판정한다.

## line 2에서 몇 비트까지 복구되는가 (20,000 trace 실측)

### 고속 경로

Alg 3.18 전체를 보드에서 돌리면 13.5 s/trace라 통계적으로 의미 있는 캠페인이 불가능하다.
측정 구간은 line 2뿐이므로, **호스트에서 실제 Alg 3.18을 돌려** α를 뽑고
(`conj_vector_host alphas`, 2000개에 7.4 s) 보드는 line 2만 실행한다
(`u`/`z`/`q` 명령). 실행 코드와 트리거 위치는 동일하고 주변 명령만 다르며,
캡처 시작 시 `alpha loaded exactly, Conjugate correct`로 검증한다.
**238 ms/trace — 57배 빠르다.**

### 결과

프로파일드 템플릿, 프로파일링 16,000 / 공격 4,000 (held-out), 단일 trace 공격.

| 대상 | 정확도 | 기준선 | z | 판정 |
|---|---|---|---|---|
| **부호 (좌표당 1 bit)** | **100.0%** | 50% | +62 | **완전 복구** |
| MSB limb의 HW | 26.6–28.6% | 11–14% | +23~33 | 복구 |
| LSB limb의 HW (대조) | 19.7–21.9% | 14% | +11~15 | 복구 |
| magnitude bit-length | 27.4% | 27.5% | −0.0 | 기준선 |
| **상위 T비트 (T = 1…8)** | 기준선과 동일 | | | **복구 안 됨** |

HW 정보량 합계는 **18 limb에서 10.35 bit** (천장 64 bit).

```
복구됨:  부호 4 bit (완전) + HW 정보 28.75 bit (신경망)  ->  약 33 bit
복구 안 됨:  상위 비트 0 bit,  개별 비트 0.04 bit/limb (무의미)
전체 비밀:  614 bit
```

HW 정보는 탐색공간을 2^28.75 배 줄이는 것이지 지목해서 맞춘 비트가 아니다.
지목해서 맞춘 비트는 부호 4개뿐이다.

### 딥러닝 공격 (dl_recover_line2.py)

템플릿과 LDA는 누설 모델(HW)을 가정한다. 신경망은 가정하지 않으므로 HW 외의 것이
있다면 찾아낸다. MLP, 학습 14,000 / 검증 2,000 / 테스트 4,000, 에폭은 검증 split으로
선택(테스트셋으로 모델을 고르면 낙관 편향이 들어간다), 레이블 섞기 대조군 포함.

| 대상 | 템플릿 | **신경망** |
|---|---:|---:|
| HW 정보량 (18 limb 합) | 10.35 bit | **28.75 bit** |
| 개별 비트 (32비트 limb 1개) | — | 0.044 bit |
| **상위 T비트 (T=1…8)** | 0 | **0** |

**HW 추출은 2.8배 좋아진다** (천장 64 bit의 16% → 45%). 18 limb 중 15개에서 1.5–2.2 bit씩
뽑는다. 모델 가정 없이 여러 샘플을 비선형 결합하는 쪽이 가우시안 템플릿보다 낫다.

**MSB는 신경망으로도 0이다.** 네 좌표 × T=1…8 전부 z ≈ 0 이하, I = 0.000.

### 왜 MSB가 안 나오는가

HW는 복구되는데 **MSB 비트로 변환되지 않는다.** 두 가지가 동시에 막는다.

1. **magnitude bit-length가 복구되지 않는다** (z ≈ 0). MSB가 몇 번째 비트인지 모르면
   상위 비트를 지정할 수 없다. `bitlen`은 401로 고정이라 limb 개수가 값에 무관하고,
   magnitude는 상위 limb의 *값*만 바꾸므로 구조적 흔적을 남기지 않는다.
2. **HW는 limb 안에서 어느 비트가 켜졌는지 말해주지 않는다.** HW는 비트에 대해
   대칭이므로, 상위 limb의 HW를 알아도 그 안의 비트 배치는 결정되지 않는다.

앞서 확인한 "누설 함수가 `HW(32-bit word)` 하나뿐"이라는 결과와 일관된다.

### 상위 비트 기반 복원 알고리즘에 쓸 수 있는가

**이 연산만으로는 어렵다.** 상위 비트가 0 bit이므로 MSB를 입력으로 받는 알고리즘에
줄 것이 없다. 다만 **부호 4 bit는 100% 확실하게** 나오므로, 부호를 제약으로 쓸 수 있는
형태라면 활용 가능하다.

### 실행

```bash
./conj_vector_host alphas 20000 run1 > /tmp/alphas.txt
make clean PLATFORM=CW308_STM32F4 && make PLATFORM=CW308_STM32F4 CONJ_TRIG_MODE=2
~/.cwvenv/bin/python3 capture_line2_pico.py --alphas /tmp/alphas.txt -n 20000
~/.cwvenv/bin/python3 msb_recover_line2.py   --npz /tmp/sqisign_line2.npz
~/.cwvenv/bin/python3 recover_conj_template.py --npz /tmp/sqisign_line2.npz --no-bits
```

## 두 번째 타겟: Sign Algorithm 3.23 line 12 (복구 성공)

Alg 3.18 line 2가 복구 불가로 판명된 뒤, 알려진 값이 비밀과 만나는 지점을 찾아 옮긴
타겟이다. Alg 3.18 line 5(`IdealMultiplication`)도 검토했으나 `I'_sk,chl`과 `I_com`이
**둘 다 비밀**이라(서명에서 `E_com`을 복원해도 이데알은 못 얻는다) 같은 벽이었다.

```c
// sign.c :: compute_challenge_ideal_signature
ibz_copy_digits(&vec.v[1], sig->chall_coeff, NWORDS_ORDER);  // 서명에 포함 = 공개
ibz_mat_2x2_eval(&vec, &sk->mat_BAcan_to_BA0_two, &vec);     // 비밀 행렬 × 공개
  └─ dim2.c:50  ibz_mul(&prod, &mat->m[0][1], &vec->v[1]);   // ★ 타겟
       └─ mp.c:455  prod_ij = (ddigit_t)a.limbs[i] * b.limbs[j];
                    res.limbs[0] = low32(a0 * b0)
```

| | Alg 3.18 line 2 | **Sign line 12** |
|---|---|---|
| 비밀 | 서명마다 새로 뽑힘 | **고정** (비밀키 일부) |
| 알려진 입력 | 없음 | `chall_coeff` |
| 누설 함수 | `HW(limb)` 뿐 | `HW(low32(a₀·b₀))` |
| TVLA max\|t\| | 9.92 | **22.02** (대조군 4.03) |
| trace 비용 | 13.5 s | **65 ms** |
| 결과 | 610 bit 중 64 bit 천장 | **32 bit limb 복구** |

### 측정

- 보드 결과가 호스트와 일치(`c1/c2 match`), 타겟 곱셈 1560 cycle = 211 µs
- 정답 모델로 16개 limb 조합을 훑으면 `(i=0, j=0)`만 sample 9904에서 `|corr| = 0.38`,
  나머지 15개는 노이즈 하한. `HW(b0)`를 편상관으로 제거해도 0.3786으로 유지되므로
  공개값 누설의 전이가 아니라 곱 자체의 누설이다
- 참고로 `HW(b0)` 단독은 0.43으로 더 강하게 새지만 비밀을 담지 않는다

### 공격: 2라운드 후보 목록 CPA

단순 CPA는 실패한다. 세 단계를 거쳤다.

1. **1바이트 가설(256개) — 실패.** 32비트 중 8비트만 예측해 ghost에 진다.
   trace를 100→500으로 늘려도 순위가 3~5위에서 **개선되지 않는다**. 통계가 아니라 모델
   한계다.
2. **16비트 가설(65536개) — 1위는 못 잡는다.** 10,000 trace까지 가도 **rank 4 고정**.
   `HW(b0)` 공통 성분을 투영 제거해도 그대로다.
3. **후보 목록 + 전체 32비트 판별 — 성공.** 정답이 1위는 아니어도 **항상 상위 4위 안**에
   있다는 점을 쓴다. 1라운드에서 상위 K개를 남기고, 각각을 상위 16비트로 완성해 전체
   32비트 모델로 점수를 매기면 정답이 1위로 떠오른다.

| trace 수 | 1라운드 정답 순위 | 2라운드 결과 |
|---:|---:|---|
| 250 | 2 / 65536 | 실패 |
| **500** | 4 / 65536 | **복구 성공** |
| 1000 / 2000 / 4000 | 4 / 65536 | 복구 성공 |

K = 16과 K = 256 모두 성공하므로 후보는 16개로 충분하다.

### 남은 범위

복구한 것은 `Msk[0][1]`의 **limb 0 (하위 32비트)** 이다. 나머지 limb은
`res[1] = low32(a₀b₁ + a₁b₀ + carry)`에서 `a₀`를 안 상태로 같은 절차를 반복하면 되지만
아직 구현하지 않았다. `Msk[1][1]`은 트리거 모드 13으로 같은 방식이 적용된다.

### 사용법

```bash
make clean PLATFORM=CW308_STM32F4
make PLATFORM=CW308_STM32F4 CONJ_TRIG_MODE=12
~/.cwvenv/bin/python3 capture_chall_pico.py --mode tvla   --traces 1000 --skip-flash
~/.cwvenv/bin/python3 capture_chall_pico.py --mode random --traces 2000 --skip-flash \
    --timebase 3 --samples 30000 --out /tmp/sqisign_chall_cpa.npz
~/.cwvenv/bin/python3 cpa_chall.py --npz /tmp/sqisign_chall_cpa.npz --find-poi --curve 500,1000,2000
```

## 스코프 설정에서 주의할 점

- **Ch A 레인지**: CW308 션트 출력이 ±30 mV 수준이라 기본 ±500 mV로 잡으면 ADC를 5%만
  쓰고 파형이 노이즈로만 보인다. `--cha-range 3` (±100 mV)이 기본값이다.
- **Ch B 레인지**: 트리거 GPIO가 3.3 V라 ±2 V 레인지는 overflow한다. `--chb-range 8`
  (±5 V)이 기본값이다. 엣지만 쓰므로 치명적이진 않지만 overflow 플래그가 계속 뜬다.
- **타겟 클럭**: CW-Lite `default_setup()`은 192 MHz / 26 = **7.3846 MHz**를 낸다.
  흔히 인용되는 7.3728 MHz가 아니다 (캡처 스펙트럼 피크와 일치). 스크립트는
  `rig.set_target_clock(scope)`로 스코프에서 실측값을 읽는다.
- **auto-trigger 없음**: line 1이 먼저 수 초간 돌기 때문에 `auto_trigger_ms=0`으로
  하드웨어 트리거를 무한 대기한다. 기존 리그의 5초 auto-trigger를 그대로 쓰면 line 1
  구간이 잡힌다.
- **`--skip-flash` 사용 시 nRST 필요**: 재프로그래밍하지 않으면 이전 세션의 펌웨어가
  프로토콜 중간 상태로 남아 있고 `scope.default_setup()`이 그 밑에서 UART를 다시
  초기화하므로 첫 명령이 깨진다. `rig.reset_target(scope)`가 nRST를 펄스한다.

## 처리량

trace 1개 = Alg 3.18 1회 = **약 12–14.4초** (벡터 생성·전송·결과 읽기 포함).
TVLA 120 trace ≈ 29분, SNR 200 trace ≈ 48분.

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
