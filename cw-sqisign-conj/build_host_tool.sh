#!/bin/bash
# Build the host-side vector generator / KAT checker for the Alg 3.18 line 2
# Conjugate experiment. Uses the same vendored sources as the firmware, so host
# and board agree bit for bit.
set -euo pipefail

APP_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SQ="${APP_DIR}/vendor/sqisign/p324_3"
OUT="${APP_DIR}/conj_vector_host"

# Same flags upstream's scripts/gen_pqm4_sources.sh applies to this source set,
# minus TARGET_ARM, plus CONJ_HOST_BUILD so the trigger calls compile out.
CPPFLAGS=(
  -DRADIX_32
  -DRADIX=32
  -DSQISIGN_BUILD_TYPE_REF
  -DSQISIGN_GF_IMPL_REF
  -DSQISIGN_VARIANT=p324_3
  -DNDEBUG
  -DDISABLE_NAMESPACING
  -DSQISIGN_SINGLE_THREADED
  -DENABLE_SIGN
  -DCONJ_HOST_BUILD
)

# Quaternion + mp + prng reach Algorithm 3.18 on their own; the EC/GF files below
# are needed only so the host can call id2iso_kernel_dlogs_to_ideal_even() to build
# the challenge ideal exactly as sign.c does. The firmware receives that ideal as
# data and links none of them.
SRC=(
  "${APP_DIR}/host/conj_vector_host.c"
  "${APP_DIR}/conj_core.c"
  "${SQ}/algebra.c"
  "${SQ}/dim2.c"
  "${SQ}/dim4.c"
  "${SQ}/integers.c"
  "${SQ}/lattice.c"
  "${SQ}/finit.c"
  "${SQ}/ideal.c"
  "${SQ}/protocol.c"
  "${SQ}/qlapoty.c"
  "${SQ}/stdorder.c"
  "${SQ}/lll_applications.c"
  "${SQ}/lll_lg2.c"
  "${SQ}/lll_dim4.c"
  "${SQ}/lehmer_xgcd.c"
  "${SQ}/gaussian_xgcd.c"
  "${SQ}/mp.c"
  "${SQ}/modqx.c"
  "${SQ}/quaternion_data.c"
  "${SQ}/torsion_constants.c"
  "${SQ}/prng.c"
  "${SQ}/fips202.c"
  "${SQ}/mem.c"
  "${SQ}/id2iso.c"
  "${SQ}/endomorphism_action.c"
  "${SQ}/ec.c"
  "${SQ}/ec_params.c"
  "${SQ}/basis.c"
  "${SQ}/e0_basis.c"
  "${SQ}/biextension.c"
  "${SQ}/fp2.c"
  "${SQ}/fp.c"
  "${SQ}/fp_p324_3.c"
  "${APP_DIR}/vendor/sqisign/host-helpers/helpers.c"
  "${APP_DIR}/vendor/sqisign/host-helpers/random_input_generation.c"
)

CFLAGS=(-O2 -std=c11 -Wall -Wno-unused-function)
INC=(-I"${APP_DIR}" -I"${SQ}" -I"${APP_DIR}/vendor/sqisign/host-helpers")

echo "building ${OUT}"
cc "${CFLAGS[@]}" "${CPPFLAGS[@]}" "${INC[@]}" "${SRC[@]}" -o "${OUT}" -lm
echo "ok: ${OUT}"
