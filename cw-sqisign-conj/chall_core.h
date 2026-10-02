/* Instrumented SQIsign Algorithm 3.23 (Sign) line 12:
 *
 *     (c1, c2) <- Msk * [1, achl]^T
 *
 * Spec:  sqisign-20260901.pdf, Algorithm 3.23 line 12 (p.23)
 * Upstream: the-sqisign.git f417ebd
 *             sign.c::compute_challenge_ideal_signature  -> ibz_mat_2x2_eval()
 *             dim2.c:42  ibz_mat_2x2_eval
 *             mp.c:431   ibz_mul, schoolbook over 32-bit limbs
 *
 * Why this target rather than Algorithm 3.18 line 2
 * -------------------------------------------------
 * alpha_sk,chl is redrawn for every signature, so no amount of averaging helps
 * and the leakage ceiling is whatever a single trace carries. Msk is different:
 * it is part of the secret key and is the same in every signature, while
 * achl = chall_coeff comes from HASH(pk, E_com, msg) and ships inside the
 * signature. A fixed secret multiplied by a known, varying input is the
 * structure CPA needs, and the schoolbook multiply gives a byte-serial
 * divide-and-conquer: (a*b) mod 2^8 depends only on the low bytes of a and b.
 *
 * Trigger modes (CONJ_TRIG_MODE, continuing the table in conj_core.h):
 *   10  the whole ibz_mat_2x2_eval
 *   11  one pulse per ibz_mul             -> 4 aligned segments
 *   12  ibz_mul(Msk[0][1], achl) only     -- primary target
 *   13  ibz_mul(Msk[1][1], achl) only
 *   14  one pulse per outer-loop iteration of the ibz_mul at [0][1], i.e. per
 *       secret limb -> 4 segments, one per limb of the secret entry
 *
 * Msk[0][0] and Msk[1][0] multiply the constant 1, so they carry no
 * known-input structure; modes 12 and 13 target the two entries that meet achl.
 */
#ifndef CW_SQISIGN_CHALL_CORE_H
#define CW_SQISIGN_CHALL_CORE_H

#include <stdint.h>

#include <quaternion.h>

typedef struct
{
    uint32_t eval_cycles;    /* the whole ibz_mat_2x2_eval */
    uint32_t mul_cycles[4];  /* each ibz_mul: [0][0], [0][1], [1][0], [1][1] */
    uint32_t secret_limbs;   /* NUM_LIMBS of the targeted secret entry */
    uint32_t known_limbs;    /* NUM_LIMBS of achl */
    uint32_t secret_bitlen;
    uint32_t known_bitlen;
} chall_profile_t;

extern volatile chall_profile_t chall_profile;

void chall_profile_reset(void);

/* res = mat * vec, instrumented. Mirrors upstream ibz_mat_2x2_eval (dim2.c:42)
 * operation for operation; only triggers and cycle counts are added. */
void chall_mat_2x2_eval_triggered(ibz_vec_2_t *res, const ibz_mat_2x2_t *mat, const ibz_vec_2_t *vec);

#endif /* CW_SQISIGN_CHALL_CORE_H */
