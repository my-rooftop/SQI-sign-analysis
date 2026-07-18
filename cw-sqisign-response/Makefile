TARGET = simpleserial-sqisign-response
APP_DIR := $(dir $(abspath $(lastword $(MAKEFILE_LIST))))
SQ_ROOT ?= $(APP_DIR)/vendor/sqisign/src

SQISIGN_LEVEL ?= lvl1
LVL_PRECOMP = $(SQ_ROOT)/precomp/ref/$(SQISIGN_LEVEL)

SRC += simpleserial-sqisign-response.c
SRC += response_core.c
SRC += cw_randombytes.c

VPATH += :$(SQ_ROOT)/common/generic
VPATH += :$(SQ_ROOT)/mini-gmp
VPATH += :$(SQ_ROOT)/quaternion/ref/generic
VPATH += :$(SQ_ROOT)/quaternion/ref/generic/lll
VPATH += :$(SQ_ROOT)/quaternion/ref/generic/hnf
VPATH += :$(LVL_PRECOMP)

SRC += tools.c
SRC += mini-gmp.c
SRC += mini-gmp-extra.c
SRC += intbig.c
SRC += dim2.c
SRC += dim4.c
SRC += lattice.c
SRC += ideal.c
SRC += lat_ball.c
SRC += algebra.c
SRC += finit.c
SRC += hnf.c
SRC += hnf_internal.c
SRC += ibz_division.c
SRC += l2.c
SRC += lll_verification.c
SRC += quaternion_data.c

EXTRAINCDIRS += .
EXTRAINCDIRS += $(APP_DIR)/vendor/sqisign/include
EXTRAINCDIRS += $(SQ_ROOT)/common/generic/include
EXTRAINCDIRS += $(SQ_ROOT)/mini-gmp
EXTRAINCDIRS += $(SQ_ROOT)/precomp/ref/$(SQISIGN_LEVEL)/include
EXTRAINCDIRS += $(SQ_ROOT)/quaternion/ref/generic/include
EXTRAINCDIRS += $(SQ_ROOT)/quaternion/ref/generic/internal_quaternion_headers

CRYPTO_TARGET = NONE
OPT = s
SS_VER ?= SS_VER_1_1
RESP_TRIGGER_MODE ?= 0
RESP_PROBE_REPS ?= 1

CFLAGS_LAST += -DMINI_GMP
CFLAGS_LAST += -DGMP_LIMB_BITS=32
CFLAGS_LAST += -DRADIX_32
CFLAGS_LAST += -DSQISIGN_BUILD_TYPE_REF
CFLAGS_LAST += -DSQISIGN_GF_IMPL_REF
CFLAGS_LAST += -DSQISIGN_VARIANT=$(SQISIGN_LEVEL)
CFLAGS_LAST += -DTARGET_ARM
CFLAGS_LAST += -DTARGET_OS_OTHER
CFLAGS_LAST += -DNDEBUG
CFLAGS_LAST += -DDISABLE_NAMESPACING
CFLAGS_LAST += -DRESP_TRIGGER_MODE=$(RESP_TRIGGER_MODE)
CFLAGS_LAST += -DRESP_PROBE_REPS=$(RESP_PROBE_REPS)

FIRMWAREPATH ?= $(APP_DIR)/vendor/chipwhisperer-firmware

include $(FIRMWAREPATH)/Makefile.inc
