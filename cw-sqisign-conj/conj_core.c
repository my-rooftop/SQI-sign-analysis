#include "conj_core.h"

#include <internal.h>
#include <mp_internal.h>
#include <quaternion_constants.h>
#include <quaternion_data.h>

#ifdef CONJ_HOST_BUILD
#define trigger_high() ((void)0)
#define trigger_low() ((void)0)
#else
#include "hal.h"
#endif

#define CONJ_NUM_LIMBS(bits) ((int)(((bits) + (NUM_BITS_LIMB - 1)) / (NUM_BITS_LIMB)))

volatile conj_profile_t conj_profile;
static volatile uint32_t conj_probe_sink;

#if CONJ_OUT_OF_PLACE
static quat_alg_elem_t conj_scratch;
#endif

/* ------------------------------------------------------------------ cycles */

#define DWT_CTRL (*(volatile uint32_t *)0xE0001000u)
#define DWT_CYCCNT (*(volatile uint32_t *)0xE0001004u)
#define DEMCR (*(volatile uint32_t *)0xE000EDFCu)

void
conj_profile_reset(void)
{
    volatile uint8_t *p = (volatile uint8_t *)&conj_profile;
    for (unsigned i = 0; i < sizeof(conj_profile); i++) {
        p[i] = 0;
    }
}

void
conj_cycles_init(void)
{
#ifndef CONJ_HOST_BUILD
    DEMCR |= (1u << 24);
    DWT_CYCCNT = 0;
    DWT_CTRL |= 1u;
#endif
}

uint32_t
conj_cycles_read(void)
{
#ifdef CONJ_HOST_BUILD
    return 0;
#else
    return DWT_CYCCNT;
#endif
}

static uint32_t
cycles_since(uint32_t start)
{
    return conj_cycles_read() - start;
}

/* --------------------------------------------------------------- probes */

/* Artificial limb walk over one coordinate. Mirrors the memory access pattern
 * ibz_neg performs, repeated CONJ_PROBE_REPS times to raise SNR. Used by mode 7
 * to establish the leakage ceiling for the real single-pass Conjugate. */
static void
probe_ibz_limbs(const ibz_t *x)
{
    int n = CONJ_NUM_LIMBS(x->bitlen);
    uint32_t acc = (uint32_t)n;

    trigger_high();
    for (int rep = 0; rep < CONJ_PROBE_REPS; rep++) {
        for (int i = 0; i < n; i++) {
            acc ^= (uint32_t)x->limbs[i];
            acc = (acc << 3) | (acc >> 29);
        }
    }
    trigger_low();

    conj_probe_sink ^= acc;
}

/* Byte-for-byte equivalent of upstream ibz_neg (mp.c:469) with one trigger pulse
 * per limb, so a single capture can be aligned against the limb index. Only used
 * by mode 8; every other mode calls the vendored ibz_neg. */
static void
ibz_neg_per_limb_triggered(ibz_t *neg, const ibz_t *a)
{
    int n = CONJ_NUM_LIMBS(a->bitlen);
    digit_t carry = 1;
    for (int i = 0; i < n; i++) {
        trigger_high();
        ddigit_t t = (ddigit_t)~a->limbs[i] + carry;
        neg->limbs[i] = (digit_t)t;
        carry = (digit_t)(t >> NUM_BITS_LIMB);
        trigger_low();
    }
    neg->bitlen = a->bitlen;
}

/* ------------------------------------------- Algorithm 3.18, line 2 target */

/* Conjugate (Algorithm 4.38): (a + ib + jc + kd)/e -> (a - ib - jc - kd)/e.
 *
 * Expanded from quat_alg_conj (upstream algebra.c:168) so the three ibz_neg
 * calls can be individually triggered and profiled. The arithmetic performed is
 * identical to the vendored function.
 */
static void
conjugate_triggered(quat_alg_elem_t *conj, const quat_alg_elem_t *x)
{
    uint32_t t;

    for (int i = 0; i < 4; i++) {
        conj_profile.conj_bitlen[i] = (uint32_t)x->coord.v[i].bitlen;
        conj_profile.conj_limbs[i] = (uint32_t)CONJ_NUM_LIMBS(x->coord.v[i].bitlen);
    }

#if CONJ_TRIG_MODE == 2
    trigger_high();
#endif

#if CONJ_TRIG_MODE == 4
    trigger_high();
#endif
    ibz_copy(&conj->denom, &x->denom);
#if CONJ_TRIG_MODE == 4
    trigger_low();
    trigger_high();
#endif
    t = conj_cycles_read();
    ibz_copy(&conj->coord.v[0], &x->coord.v[0]);
    conj_profile.conj_coord_cycles[0] = cycles_since(t);
#if CONJ_TRIG_MODE == 4
    trigger_low();
#endif

    for (int i = 1; i < 4; i++) {
#if CONJ_TRIG_MODE == 3 || CONJ_TRIG_MODE == 4
        trigger_high();
#endif
        t = conj_cycles_read();
#if CONJ_TRIG_MODE == 8
        if (i == 1) {
            ibz_neg_per_limb_triggered(&conj->coord.v[i], &x->coord.v[i]);
        } else {
            ibz_neg(&conj->coord.v[i], &x->coord.v[i]);
        }
#else
        ibz_neg(&conj->coord.v[i], &x->coord.v[i]);
#endif
        conj_profile.conj_coord_cycles[i] = cycles_since(t);
#if CONJ_TRIG_MODE == 3 || CONJ_TRIG_MODE == 4
        trigger_low();
#endif
    }

#if CONJ_TRIG_MODE == 2
    trigger_low();
#endif
}

void
conj_run_conjugate(quat_alg_elem_t *x, quat_alg_elem_t *pre)
{
    if (pre != NULL) {
        quat_alg_elem_copy(pre, x);
    }
    uint32_t t = conj_cycles_read();
#if CONJ_OUT_OF_PLACE
    conjugate_triggered(&conj_scratch, x);
    quat_alg_elem_copy(x, &conj_scratch);
#else
    conjugate_triggered(x, x);
#endif
    conj_profile.conjugate_cycles = cycles_since(t);
}

/* ---------------------------------------------- Algorithm 3.18, full flow */

void
conj_response_element(quat_ideal_t *sk_chall_ideal,
                      quat_alg_elem_t *resp_quat,
                      ibz_t *norm,
                      quat_alg_elem_t *sk_chall_quat,
                      quat_alg_elem_t *sk_chall_quat_pre,
                      const quat_ideal_t *sk_ideal,
                      const quat_ideal_t *ideal_chall_two,
                      const quat_alg_elem_t *chall_split,
                      const quat_ideal_t *ideal_commit,
                      prng_domain_ctx_t *prng_domain,
                      int e)
{
    ibz_t lattice_content, bound;
    quat_lattice_t prod;
    ibz_init(&lattice_content);
    ibz_init(&bound);
    quat_lattice_init(&prod);

    uint32_t t_total = conj_cycles_read();
    uint32_t t;

#if CONJ_TRIG_MODE == 0
    trigger_high();
#endif

    /* Intersect (Algorithm 3.17): sk_chall_ideal = ideal_chall_two ^ sk_ideal */
    t = conj_cycles_read();
    quat_ideal_intersect_O0(sk_chall_ideal, ideal_chall_two, sk_ideal, chall_split);
    conj_profile.intersect_cycles = cycles_since(t);

    /* line 1: I'_sk,chl, alpha_sk,chl <- EquivalentCoprimeIdeal(I_sk,chl, 0, BitLenSkChall, prng_def) */
#if CONJ_TRIG_MODE == 1
    trigger_high();
#endif
    t = conj_cycles_read();
    int resp = quat_ideal_small_equivalent_coprime(
        sk_chall_quat, sk_chall_ideal, sk_chall_ideal, &ibz_const_zero, &PRNG_default_domain);
    conj_profile.equiv_coprime_cycles = cycles_since(t);
    conj_profile.equiv_ok = (uint32_t)resp;
#if CONJ_TRIG_MODE == 1
    trigger_low();
#endif
    assert(resp);

    /* line 2: alpha_sk,chl <- Conjugate(alpha_sk,chl)   <-- measurement target */
    if (sk_chall_quat != NULL) {
        /* Ground truth for side-channel labels: the value entering line 2. */
        if (sk_chall_quat_pre != NULL) {
            quat_alg_elem_copy(sk_chall_quat_pre, sk_chall_quat);
        }

#if CONJ_TRIG_MODE == 7
        probe_ibz_limbs(&sk_chall_quat->coord.v[1]);
        probe_ibz_limbs(&sk_chall_quat->coord.v[2]);
        probe_ibz_limbs(&sk_chall_quat->coord.v[3]);
#endif

        conj_run_conjugate(sk_chall_quat, NULL);
    }

    /* line 5a: IdealMultiplication(I'_sk,chl, I_com) */
#if CONJ_TRIG_MODE == 5
    trigger_high();
#endif
    t = conj_cycles_read();
    quat_ideal_mul_O0(&prod, sk_chall_ideal, ideal_commit);
    conj_profile.ideal_mul_cycles = cycles_since(t);
#if CONJ_TRIG_MODE == 5
    trigger_low();
#endif

    ibz_div_2exp(&lattice_content, &prod.basis.m[0][0], 1);

    /* line 5b: RandomBoundedElement(..., 2^(e_rsp - 1), prng_res) */
    {
        ibz_mul_2exp(&bound, &ibz_const_one, e);
        ibz_sub(&bound, &bound, &ibz_const_one);
        ibz_set_bound(&bound, e + 1);

#if CONJ_TRIG_MODE == 6
        trigger_high();
#endif
        t = conj_cycles_read();
        int ok = quat_lattice_sample_from_ball(resp_quat, norm, &prod, &bound, prng_domain);
        conj_profile.sample_cycles = cycles_since(t);
        conj_profile.sample_ok = (uint32_t)ok;
#if CONJ_TRIG_MODE == 6
        trigger_low();
#endif
        assert(ok);
        ibz_set_bound(norm, e + 1);
    }

#if CONJ_TRIG_MODE == 0
    trigger_low();
#endif

    conj_profile.total_cycles = cycles_since(t_total);
}
