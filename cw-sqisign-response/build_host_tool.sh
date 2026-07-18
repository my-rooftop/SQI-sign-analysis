#!/usr/bin/env bash
set -euo pipefail

APP_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SQ_ROOT="$APP_DIR/vendor/sqisign"

cd "$APP_DIR"

gcc -O2 -g \
  -DMINI_GMP -DGMP_LIMB_BITS=64 -DRADIX_64 \
  -DSQISIGN_BUILD_TYPE_REF -DSQISIGN_GF_IMPL_REF -DSQISIGN_VARIANT=lvl1 \
  -DTARGET_AMD64 -DTARGET_OS_OTHER -DNDEBUG -DDISABLE_NAMESPACING -DENABLE_SIGN \
  -I"$SQ_ROOT/include" \
  -I"$SQ_ROOT/src/common/generic/include" \
  -I"$SQ_ROOT/src/mini-gmp" \
  -I"$SQ_ROOT/src/mp/ref/generic/include" \
  -I"$SQ_ROOT/src/precomp/ref/lvl1/include" \
  -I"$SQ_ROOT/src/quaternion/ref/generic/include" \
  -I"$SQ_ROOT/src/quaternion/ref/generic/internal_quaternion_headers" \
  -I"$SQ_ROOT/src/gf/ref/include" \
  -I"$SQ_ROOT/src/gf/ref/lvl1" \
  -I"$SQ_ROOT/src/gf/ref/lvlx" \
  -I"$SQ_ROOT/src/ec/ref/lvlx" \
  -I"$SQ_ROOT/src/ec/ref/include" \
  -I"$SQ_ROOT/src/hd/ref/include" \
  -I"$SQ_ROOT/src/hd/ref/lvlx" \
  -I"$SQ_ROOT/src/id2iso/ref/include" \
  -I"$SQ_ROOT/src/id2iso/ref/lvlx" \
  -I"$SQ_ROOT/src/signature/ref/include" \
  -I"$SQ_ROOT/src/signature/ref/lvlx" \
  -I"$SQ_ROOT/src/verification/ref/include" \
  -I"$SQ_ROOT/src/verification/ref/lvlx" \
  response_split_host.c \
  cw_randombytes.c \
  "$SQ_ROOT/src/common/generic/tools.c" \
  "$SQ_ROOT/src/common/generic/fips202.c" \
  "$SQ_ROOT/src/common/generic/mem.c" \
  "$SQ_ROOT/src/mini-gmp/mini-gmp.c" \
  "$SQ_ROOT/src/mini-gmp/mini-gmp-extra.c" \
  "$SQ_ROOT/src/mp/ref/generic/mp.c" \
  "$SQ_ROOT/src/gf/ref/lvlx/fp.c" \
  "$SQ_ROOT/src/gf/ref/lvlx/fp2.c" \
  "$SQ_ROOT/src/gf/ref/lvl1/fp_p5248_64.c" \
  "$SQ_ROOT/src/precomp/ref/lvl1/ec_params.c" \
  "$SQ_ROOT/src/precomp/ref/lvl1/e0_basis.c" \
  "$SQ_ROOT/src/precomp/ref/lvl1/hd_splitting_transforms.c" \
  "$SQ_ROOT/src/precomp/ref/lvl1/torsion_constants.c" \
  "$SQ_ROOT/src/precomp/ref/lvl1/quaternion_data.c" \
  "$SQ_ROOT/src/precomp/ref/lvl1/endomorphism_action.c" \
  "$SQ_ROOT/src/ec/ref/lvlx/ec.c" \
  "$SQ_ROOT/src/ec/ref/lvlx/ec_jac.c" \
  "$SQ_ROOT/src/ec/ref/lvlx/xisog.c" \
  "$SQ_ROOT/src/ec/ref/lvlx/xeval.c" \
  "$SQ_ROOT/src/ec/ref/lvlx/isog_chains.c" \
  "$SQ_ROOT/src/ec/ref/lvlx/basis.c" \
  "$SQ_ROOT/src/ec/ref/lvlx/biextension.c" \
  "$SQ_ROOT/src/hd/ref/lvlx/hd.c" \
  "$SQ_ROOT/src/hd/ref/lvlx/theta_structure.c" \
  "$SQ_ROOT/src/hd/ref/lvlx/theta_isogenies.c" \
  "$SQ_ROOT/src/verification/ref/lvlx/common.c" \
  "$SQ_ROOT/src/verification/ref/lvlx/encode_verification.c" \
  "$SQ_ROOT/src/verification/ref/lvlx/verify.c" \
  "$SQ_ROOT/src/id2iso/ref/lvlx/id2iso.c" \
  "$SQ_ROOT/src/id2iso/ref/lvlx/dim2id2iso.c" \
  "$SQ_ROOT/src/signature/ref/lvlx/encode_signature.c" \
  "$SQ_ROOT/src/signature/ref/lvlx/keygen.c" \
  "$SQ_ROOT/src/quaternion/ref/generic/intbig.c" \
  "$SQ_ROOT/src/quaternion/ref/generic/algebra.c" \
  "$SQ_ROOT/src/quaternion/ref/generic/ideal.c" \
  "$SQ_ROOT/src/quaternion/ref/generic/dim4.c" \
  "$SQ_ROOT/src/quaternion/ref/generic/dim2.c" \
  "$SQ_ROOT/src/quaternion/ref/generic/integers.c" \
  "$SQ_ROOT/src/quaternion/ref/generic/lattice.c" \
  "$SQ_ROOT/src/quaternion/ref/generic/lat_ball.c" \
  "$SQ_ROOT/src/quaternion/ref/generic/finit.c" \
  "$SQ_ROOT/src/quaternion/ref/generic/printer.c" \
  "$SQ_ROOT/src/quaternion/ref/generic/lll/rationals.c" \
  "$SQ_ROOT/src/quaternion/ref/generic/lll/l2.c" \
  "$SQ_ROOT/src/quaternion/ref/generic/lll/lll_verification.c" \
  "$SQ_ROOT/src/quaternion/ref/generic/lll/lll_applications.c" \
  "$SQ_ROOT/src/quaternion/ref/generic/normeq.c" \
  "$SQ_ROOT/src/quaternion/ref/generic/hnf/ibz_division.c" \
  "$SQ_ROOT/src/quaternion/ref/generic/hnf/hnf_internal.c" \
  "$SQ_ROOT/src/quaternion/ref/generic/hnf/hnf.c" \
  "$SQ_ROOT/src/quaternion/ref/generic/test/random_input_generation.c" \
  -lm \
  -o response_split_host
