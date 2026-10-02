TARGET = simpleserial-sqisign-conj
APP_DIR := $(dir $(abspath $(lastword $(MAKEFILE_LIST))))

# Vendored upstream embedded source set, produced by the-sqisign.git
# scripts/gen_pqm4_sources.sh -- see vendor/sqisign/PROVENANCE.txt.
SQISIGN_VARIANT_DIR ?= p324_3
SQ := $(APP_DIR)/vendor/sqisign/$(SQISIGN_VARIANT_DIR)

SRC += simpleserial-sqisign-conj.c
SRC += conj_core.c
SRC += chall_core.c
SRC += cw_randombytes.c

VPATH += :$(SQ)

# Reaching Algorithm 3.18 needs the quaternion module, the 2026 mp module (which
# replaced mini-gmp), the SHAKE-based PRNG and the level's quaternion constants.
# No EC, GF, isogeny or signature-packaging code is linked: the challenge ideal
# is built on the host and sent over as data.
SRC += algebra.c
SRC += dim2.c
SRC += dim4.c
SRC += integers.c
SRC += lattice.c
SRC += finit.c
SRC += ideal.c
SRC += protocol.c
SRC += qlapoty.c
SRC += stdorder.c
SRC += lll_applications.c
SRC += lll_lg2.c
SRC += lll_dim4.c
SRC += lehmer_xgcd.c
SRC += gaussian_xgcd.c
SRC += mp.c
SRC += modqx.c
SRC += quaternion_data.c
SRC += torsion_constants.c
SRC += prng.c
SRC += fips202.c
SRC += mem.c

EXTRAINCDIRS += .
EXTRAINCDIRS += $(SQ)

CRYPTO_TARGET = NONE
OPT = s
SS_VER ?= SS_VER_1_1

# 2  = Algorithm 3.18 line 2, Conjugate(alpha_sk,chl)        -- see conj_core.h
# 12 = Algorithm 3.23 line 12, ibz_mul(Msk[0][1], achl)      -- see chall_core.h
CONJ_TRIG_MODE ?= 2
CONJ_PROBE_REPS ?= 1
CONJ_OUT_OF_PLACE ?= 0

# Same definitions upstream's gen_pqm4_sources.sh applies to this source set.
# fips202.c is compiled here rather than supplied by pqm4, so the
# shake256_inc_ctx_reset alias that config.mk sets is deliberately omitted.
CFLAGS_LAST += -DRADIX_32
CFLAGS_LAST += -DRADIX=32
CFLAGS_LAST += -DSQISIGN_BUILD_TYPE_REF
CFLAGS_LAST += -DSQISIGN_GF_IMPL_REF
CFLAGS_LAST += -DSQISIGN_VARIANT=$(SQISIGN_VARIANT_DIR)
CFLAGS_LAST += -DTARGET_ARM
CFLAGS_LAST += -DTARGET_OS_OTHER
CFLAGS_LAST += -DNDEBUG
CFLAGS_LAST += -DDISABLE_NAMESPACING
CFLAGS_LAST += -DSQISIGN_SINGLE_THREADED
CFLAGS_LAST += -DENABLE_SIGN
CFLAGS_LAST += -DCONJ_TRIG_MODE=$(CONJ_TRIG_MODE)
CFLAGS_LAST += -DCONJ_PROBE_REPS=$(CONJ_PROBE_REPS)
CFLAGS_LAST += -DCONJ_OUT_OF_PLACE=$(CONJ_OUT_OF_PLACE)

FIRMWAREPATH ?= $(APP_DIR)/vendor/chipwhisperer-firmware

include $(FIRMWAREPATH)/Makefile.inc
