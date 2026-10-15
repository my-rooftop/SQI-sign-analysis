# CW SQIsign response 실험 패키지

이 폴더는 SQIsign level 1 서명 과정 중 response quaternion 생성 단계인
`RandomEquivalentQuaternion(...)` 코어를 CW308 STM32F405에서 실행하고,
결과 검증·사이클 프로파일링·전력 파형 수집을 하기 위한 독립 실행형 실험
패키지다.

SQIsign 전체 서명 구현을 보드에서 실행하는 프로그램은 아니다. 보드가 받는
입력은 이미 만들어진 secret/challenge/commit ideal 세 개이며, 출력은 다음
response quaternion과 검증 보조값이다.

```text
alpha_rsp = (coord0 + coord1*i + coord2*j + coord3*k) / denom
보조값    = lattice_content
```

상세 알고리즘, 측정 위치, 기존 측정값은
[`RESPONSE_FLOW_AND_MEASUREMENTS.md`](RESPONSE_FLOW_AND_MEASUREMENTS.md)에 정리되어
있다.

> 주의: `cw_randombytes.c`는 재현 가능한 실험을 위한 단순 결정론적 RNG다.
> 암호학적으로 안전하지 않으므로 실제 키나 운영 환경에 사용하면 안 된다.

## 폴더 구성

```text
firmware/   STM32F4 타깃 펌웨어 (arm-none-eabi, 32-bit limb)
host/       PC 참조 구현 및 검증 도구 (gcc, 64-bit limb)
scripts/    파형 수집·실험·분석 Python 스크립트
vendor/     포함된 SQIsign / ChipWhisperer 소스와 라이선스
```

`cw_randombytes.c`는 `firmware/`에 있지만 호스트 도구도 같은 파일을 컴파일한다.
PC와 보드가 같은 결정론적 실험 RNG를 쓰게 하려는 의도적인 공유다.

## 폴더 독립성

`cw-sqisign-response/` 폴더만 다른 위치로 복사해도 상위 SQIsign 저장소 없이
호스트 도구와 STM32F4 펌웨어를 다시 빌드할 수 있다. 필요한 SQIsign 소스와
CW308 STM32F4용 ChipWhisperer 펌웨어 소스는 `vendor/` 아래에 포함되어 있다.
폴더 이름을 바꿔도 Python 스크립트는 동작한다. 다만 `firmware/`, `host/`,
`scripts/`, `vendor/`의 상대 위치는 그대로 유지해야 한다.

단, 컴파일러, Python 패키지, USB 드라이버와 실제 ChipWhisperer 장비는 운영체제에
설치하거나 연결해야 한다.

## 의존성

확인된 개발 환경은 Ubuntu 24.04 계열, Python 3.12.3, GCC 13.3.0,
GNU Make 4.3, Arm GNU Toolchain 13.2.1이다.

필수 소프트웨어는 다음과 같다.

- 호스트 도구 빌드: `bash`, `gcc`, 표준 C 라이브러리
- 펌웨어 재빌드: `make`, `arm-none-eabi-gcc`, `arm-none-eabi-binutils`
- 보드 제어: Python 3, `chipwhisperer==6.0.0`
- 파형 저장·분석: `numpy==1.26.4`
- 장비: ChipWhisperer capture 장비, CW308 UFO baseboard, STM32F405 target board
- 권한: ChipWhisperer USB 장치에 접근할 수 있는 udev 규칙 또는 동등한 권한

Python 환경은 다음처럼 준비한다.

```bash
cd cw-sqisign-response
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
```

이 저장소에는 소스만 들어 있고 프리빌트 바이너리(`host/response_split_host`,
`.hex` 펌웨어 등)는 포함하지 않는다. 호스트 도구는 `host/build_host_tool.sh`,
펌웨어는 `make -C firmware PLATFORM=CW308_STM32F4`로 빌드한다. 아래 스크립트들은
실행 시 필요한 산출물을 자동으로 빌드한다.

## 가장 빠른 정상 동작 확인

### 1. 호스트 도구 빌드와 확인

```bash
cd cw-sqisign-response
./host/build_host_tool.sh
./host/response_split_host gen > /tmp/sqisign_response_vec.txt
./host/response_split_host compute-vector /tmp/sqisign_response_vec.txt \
  > /tmp/sqisign_response_host.txt
```

`gen`은 고정 seed로 SQIsign key generation/commit/challenge를 수행해 보드 입력
ideal 세 개를 만든다. `compute-vector`는 같은 입력에 대해 PC에서 response 코어를
실행한다. response 생성에는 rejection sampling이 들어가므로 PC와 보드의 출력값을
단순히 바이트 단위로 비교하지 않는다.

### 2. 보드 flash, 실행, 검증

장비를 연결한 뒤 다음을 실행한다.

```bash
python scripts/pc_board_response_kat.py
```

이 스크립트는 다음 작업을 한 번에 수행한다.

1. 호스트에서 입력 벡터 생성
2. STM32F4 펌웨어 빌드 후 `.hex` flash
3. SimpleSerial로 입력 전송
4. 보드 response 계산 및 프로파일 읽기
5. 결과가 올바른 lattice에 속하고 norm bound를 만족하는지 PC에서 검증

마지막 출력이 아래처럼 나오면 정상이다.

```text
content=OK
contains=OK
bound=OK
OK
```

이미 같은 펌웨어가 올라가 있으면 flash를 생략할 수 있다.

```bash
python scripts/pc_board_response_kat.py --skip-flash
```

`scripts/program_and_test.py`는 response 계산 대신 flash와 SimpleSerial `i` 명령
응답만 확인하는 짧은 연결 점검용 스크립트다.

```bash
python scripts/program_and_test.py
```

## 펌웨어 재빌드

필요한 ChipWhisperer 빌드 파일은 `vendor/chipwhisperer-firmware/`에 포함되어
있으므로 외부 `FIRMWAREPATH` 설정은 필요 없다.

```bash
make -C firmware clean PLATFORM=CW308_STM32F4
make -C firmware PLATFORM=CW308_STM32F4
```

결과 펌웨어는 다음 파일이다.

```text
firmware/simpleserial-sqisign-response-CW308_STM32F4.hex
```

다른 ChipWhisperer 펌웨어 트리를 사용해야 할 때만 경로를 덮어쓴다.

```bash
make -C firmware PLATFORM=CW308_STM32F4 \
  FIRMWAREPATH=/path/to/chipwhisperer/firmware/mcu
```

컴파일 옵션을 바꾸면 기존 object가 그대로 재사용되지 않도록 반드시 먼저
`make -C firmware clean PLATFORM=CW308_STM32F4`를 실행한다.

## 파형 캡처

기본 명령은 이미 계산된 response를 다시 직렬화하는 구간을 캡처한다.

```bash
python scripts/capture_response_leakage.py \
  --traces 20 \
  --samples 80000 \
  --out response_serialize.npz
```

스크립트는 기본적으로 trigger mode 1 펌웨어를 빌드하고 flash한다. 전체 response
계산이나 내부 세부 구간을 측정하려면 원하는 mode로 clean build한 뒤
`--skip-build --compute-each`를 함께 사용한다.

```bash
make -C firmware clean PLATFORM=CW308_STM32F4
make -C firmware PLATFORM=CW308_STM32F4 RESP_TRIGGER_MODE=2
python scripts/capture_response_leakage.py \
  --skip-build --compute-each \
  --traces 20 --samples 80000 \
  --out response_finalize.npz
```

Trigger mode는 다음과 같다.

| Mode | trigger 구간 |
|---:|---|
| 0 | response core 전체 |
| 1 | 결과 ASCII 직렬화 |
| 2 | sample finalize: `basis*x`, denom 복사, normalize |
| 3 | normalize 중 `coord0` division |
| 4 | normalize 중 `coord0` sign multiplication |
| 5 | `coord0` limb 순회용 인공 probe |

Mode 5 반복 횟수는 빌드 시 `RESP_PROBE_REPS`로 지정한다.

```bash
make -C firmware clean PLATFORM=CW308_STM32F4
make -C firmware PLATFORM=CW308_STM32F4 RESP_TRIGGER_MODE=5 RESP_PROBE_REPS=100
```

주요 캡처 옵션은 다음과 같다.

- `--compute-each`: trace마다 response를 새로 계산
- `--seed-each`: trace마다 보드 RNG seed 설정
- `--fixed-seed`: `--seed-each`와 함께 모든 trace에 같은 seed 사용
- `--offset`, `--samples`: ADC 캡처 시작 위치와 길이
- `--skip-build`, `--skip-flash`: 기존 빌드 또는 flash 재사용
- `--calibrate`: 실제 response 대신 알려진 32-bit 값의 직렬화 파형 확인

출력 `.npz`에는 `traces`, JSON 응답 기록, 각 좌표의 sign/bit length/32-bit limb와
Hamming weight label이 저장된다.

## 호스트 도구 명령

```text
host/response_split_host gen > vector.txt
host/response_split_host compute-vector vector.txt > response.txt
host/response_split_host verify vector.txt board_response.txt
host/response_split_host bench-vector vector.txt N
host/response_split_host sample N > response_stats.txt
```

- `gen`: 보드용 54개 hex 정수 입력 생성
- `compute-vector`: PC에서 주어진 ideal 세 개의 response 계산
- `verify`: 보드 출력의 content, lattice membership, norm bound 검증
- `bench-vector`: 같은 입력을 `N`회 계산해 PC cycle 통계 출력
- `sample`: 프로토콜을 `N`회 수행해 response 값 범위 통계 생성

## 주요 파일

```text
firmware/
  Makefile                         STM32F4 펌웨어 빌드
  simpleserial-sqisign-response.c  SimpleSerial 명령 및 입출력 버퍼
  response_core.c/.h               response 코어, trigger, cycle profile
  cw_randombytes.c                 실험용 결정론적 RNG (호스트 도구와 공유)

host/
  response_split_host.c            벡터 생성·PC 계산·검증 도구
  build_host_tool.sh               호스트 도구 빌드

scripts/
  program_and_test.py              flash/통신 smoke test
  pc_board_response_kat.py         보드 end-to-end KAT
  capture_response_leakage.py      ChipWhisperer ADC 파형 캡처
  capture_response_leakage_pico.py PicoScope 파형 캡처
  exp_serialize_leakage.py         직렬화 누출 (알려진 32-bit 값)
  exp_real_coords_capture.py       실제 response 좌표 캡처 (mode 6)
  exp_profile_known.py             알려진 값 대량 프로파일링 (cmd k)
  exp_repeat_avg.py                반복 평균화 이득 측정
  exp_target_finalize.py           finalize 구간 캡처 (mode 2/7)
  analyze_real_coords.py           실제 좌표 SNR·복원 분석
  analyze_profile.py               프로파일링 템플릿 학습·평가
  analyze_repeat.py                평균화 이득 분석
  plot_real_coords.py              실제 좌표 결과 플롯

vendor/
  sqisign/                         포함된 SQIsign 소스와 라이선스
  chipwhisperer-firmware/          CW308 STM32F4 빌드 최소 트리와 라이선스
```

`firmware/objdir-CW308_STM32F4/`, `.dep/`, `.elf`, `.map`, `.lss`, `.sym`,
`.bin`, `.eep`, `.hex`와 `host/response_split_host` 실행 파일은 빌드 산출물이므로
저장소에 포함하지 않으며 `.gitignore`로 제외한다. 빌드하면 위 절차대로 다시
생성된다.

## SimpleSerial 프로토콜

일반 사용자는 Python 스크립트를 쓰면 되며, 직접 호스트를 작성할 때의 명령은
다음과 같다. 모든 응답 command는 `r`이다.

| 명령 | 입력 | 기능 |
|---|---:|---|
| `i` | 0 B | 주요 구조체 크기/함수 주소 하위 값 반환 |
| `x` | 0 B | 입력·출력 상태 초기화 |
| `a` | 48 B | ASCII hex 입력 chunk 추가 |
| `c` | 0 B | 입력 parse 후 response 계산 |
| `n` | 16 B | 실험용 RNG seed 설정 |
| `s` | 0 B | 마지막 response 재직렬화 |
| `l` | 0 B | 출력 길이를 big-endian 2 B로 반환 |
| `g` | 2 B | big-endian chunk index의 출력 48 B 읽기 |
| `p` | 0 B | cycle profile을 출력 버퍼에 생성 |
| `q` | 4 B | 직렬화 calibration 값 설정 |

입력과 출력의 정확한 정수 순서 및 response 수학적 흐름은
`RESPONSE_FLOW_AND_MEASUREMENTS.md`를 참고한다.
