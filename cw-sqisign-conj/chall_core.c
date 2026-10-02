#include "chall_core.h"
#include "conj_core.h" /* CONJ_TRIG_MODE, cycle counter */

#include <internal.h>
#include <mp_internal.h>

#ifdef CONJ_HOST_BUILD
#define trigger_high() ((void)0)
#define trigger_low() ((void)0)
#else
#include "hal.h"
#endif

#define CHALL_NUM_LIMBS(bits) ((int)(((bits) + (NUM_BITS_LIMB - 1)) / (NUM_BITS_LIMB)))

volatile chall_profile_t chall_profile;

void
chall_profile_reset(void)
{
    volatile uint8_t *p = (volatile uint8_t *)&chall_profile;
    for (unsigned i = 0; i < sizeof(chall_profile); i++) {
        p[i] = 0;
    }
}

/* Schoolbook multiply with one trigger pulse per outer-loop iteration, i.e. per
 * limb of the secret operand. Arithmetic identical to upstream ibz_mul
 * (mp.c:431); used only by mode 14, which wants the capture segmented by secret
 * limb so the byte-serial attack knows which window belongs to which limb.
 *
 * Both operands here are non-negative (Msk entries are reduced mod 2^128 and
 * normalised positive by keygen, achl is a hash output), so the sign handling
 * upstream does with ibz_cneg is a no-op and is left out.
 */
static void
ibz_mul_per_limb_triggered(ibz_t *prod, const ibz_t *a, const ibz_t *b)
{
    ibz_t res = { 0 };
    int n = a->bitlen + b->bitlen - 1;
    if (n > IBZ_NLIMBS * NUM_BITS_LIMB) {
        n = IBZ_NLIMBS * NUM_BITS_LIMB;
    }
    res.bitlen = n;

    int n1 = CHALL_NUM_LIMBS(a->bitlen);
    int n2 = CHALL_NUM_LIMBS(b->bitlen);
    for (int i = 0; i < n1; i++) {
        trigger_high();
        digit_t carry = 0;
        int n3 = n2 < (IBZ_NLIMBS - i) ? n2 : (IBZ_NLIMBS - i);
        for (int j = 0; j < n3; j++) {
            ddigit_t prod_ij = (ddigit_t)a->limbs[i] * b->limbs[j];
            ddigit_t acc = (ddigit_t)res.limbs[i + j] + (digit_t)prod_ij + carry;
            res.limbs[i + j] = (digit_t)acc;
            carry = (digit_t)(prod_ij >> NUM_BITS_LIMB) + (digit_t)(acc >> NUM_BITS_LIMB);
        }
        if (i + n3 < IBZ_NLIMBS) {
            res.limbs[i + n3] = carry;
        }
        trigger_low();
    }
    *prod = res;
}

/* One multiply of the four, with the trigger placed per CONJ_TRIG_MODE.
 * `slot` is 0..3 for [0][0], [0][1], [1][0], [1][1]. */
static void
mul_slot(ibz_t *prod, const ibz_t *m, const ibz_t *v, int slot)
{
    uint32_t t = conj_cycles_read();

#if CONJ_TRIG_MODE == 11
    trigger_high();
#endif
#if CONJ_TRIG_MODE == 12
    if (slot == 1) {
        trigger_high();
    }
#endif
#if CONJ_TRIG_MODE == 13
    if (slot == 3) {
        trigger_high();
    }
#endif

#if CONJ_TRIG_MODE == 14
    if (slot == 1) {
        ibz_mul_per_limb_triggered(prod, m, v);
    } else {
        ibz_mul(prod, m, v);
    }
#else
    ibz_mul(prod, m, v);
#endif

#if CONJ_TRIG_MODE == 12
    if (slot == 1) {
        trigger_low();
    }
#endif
#if CONJ_TRIG_MODE == 13
    if (slot == 3) {
        trigger_low();
    }
#endif
#if CONJ_TRIG_MODE == 11
    trigger_low();
#endif

    chall_profile.mul_cycles[slot] = conj_cycles_read() - t;
}

void
chall_mat_2x2_eval_triggered(ibz_vec_2_t *res, const ibz_mat_2x2_t *mat, const ibz_vec_2_t *vec)
{
    ibz_t prod;
    ibz_vec_2_t matvec;
    ibz_init(&prod);
    ibz_vec_2_init(&matvec);

    chall_profile.secret_bitlen = (uint32_t)mat->m[0][1].bitlen;
    chall_profile.known_bitlen = (uint32_t)vec->v[1].bitlen;
    chall_profile.secret_limbs = (uint32_t)CHALL_NUM_LIMBS(mat->m[0][1].bitlen);
    chall_profile.known_limbs = (uint32_t)CHALL_NUM_LIMBS(vec->v[1].bitlen);

    uint32_t t_all = conj_cycles_read();
#if CONJ_TRIG_MODE == 10
    trigger_high();
#endif

    /* Same order as upstream dim2.c:48-55. */
    mul_slot(&prod, &mat->m[0][0], &vec->v[0], 0);
    ibz_copy(&matvec.v[0], &prod);
    mul_slot(&prod, &mat->m[0][1], &vec->v[1], 1);
    ibz_add(&matvec.v[0], &matvec.v[0], &prod);
    mul_slot(&prod, &mat->m[1][0], &vec->v[0], 2);
    ibz_copy(&matvec.v[1], &prod);
    mul_slot(&prod, &mat->m[1][1], &vec->v[1], 3);
    ibz_add(&matvec.v[1], &matvec.v[1], &prod);

#if CONJ_TRIG_MODE == 10
    trigger_low();
#endif
    chall_profile.eval_cycles = conj_cycles_read() - t_all;

    ibz_copy(&res->v[0], &matvec.v[0]);
    ibz_copy(&res->v[1], &matvec.v[1]);
}
