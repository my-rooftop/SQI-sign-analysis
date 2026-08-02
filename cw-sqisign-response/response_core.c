#include "response_core.h"

#include "hal.h"
#include "internal.h"

#include <encoded_sizes.h>
#include <quaternion_data.h>
#include <rng.h>

#ifndef RESP_TRIGGER_MODE
#define RESP_TRIGGER_MODE 0
#endif

#ifndef RESP_PROBE_REPS
#define RESP_PROBE_REPS 1
#endif

volatile response_profile_t response_profile;
static volatile uint32_t response_probe_sink;

#define DWT_CTRL (*(volatile uint32_t *)0xE0001000u)
#define DWT_CYCCNT (*(volatile uint32_t *)0xE0001004u)
#define DEMCR (*(volatile uint32_t *)0xE000EDFCu)

void response_profile_reset(void)
{
    volatile uint8_t *p = (volatile uint8_t *)&response_profile;
    for (unsigned i = 0; i < sizeof(response_profile); i++) {
        p[i] = 0;
    }
}

void response_cycles_init(void)
{
    DEMCR |= (1u << 24);
    DWT_CYCCNT = 0;
    DWT_CTRL |= 1u;
}

uint32_t response_cycles_read(void)
{
    return DWT_CYCCNT;
}

static uint32_t cycles_since(uint32_t start)
{
    return response_cycles_read() - start;
}

static void probe_ibz_limb_loop(const ibz_t *x)
{
    const __mpz_struct *mp = (const __mpz_struct *)(*x);
    int n = mp->_mp_size;
    if (n < 0) {
        n = -n;
    }
    uint32_t acc = (uint32_t)n;
#if RESP_TRIGGER_MODE == 5
    trigger_high();
#endif
    for (int rep = 0; rep < RESP_PROBE_REPS; rep++) {
        for (int i = 0; i < n; i++) {
            acc ^= (uint32_t)mp->_mp_d[i];
            acc = (acc << 3) | (acc >> 29);
        }
    }
#if RESP_TRIGGER_MODE == 5
    trigger_low();
#endif
    response_probe_sink ^= acc;
}

static void normalize_response_triggered(quat_alg_elem_t *x)
{
    uint32_t t_total = response_cycles_read();
    uint32_t t = t_total;
    ibz_t gcd, sign, r;
    ibz_init(&gcd);
    ibz_init(&sign);
    ibz_init(&r);

    t = response_cycles_read();
    ibz_vec_4_content(&gcd, &x->coord);
    ibz_gcd(&gcd, &gcd, &x->denom);
    ibz_div(&x->denom, &r, &x->denom, &gcd);
    response_profile.normalize_gcd_cycles = cycles_since(t);

#if RESP_TRIGGER_MODE == 3 || RESP_TRIGGER_MODE == 7
    trigger_high();
#endif
    t = response_cycles_read();
    ibz_div(&x->coord[0], &r, &x->coord[0], &gcd);
    response_profile.normalize_coord0_div_cycles = cycles_since(t);
#if RESP_TRIGGER_MODE == 3 || RESP_TRIGGER_MODE == 7
    trigger_low();
#endif
    t = response_cycles_read();
    for (int i = 1; i < 4; i++) {
#if RESP_TRIGGER_MODE == 7
        trigger_high();
#endif
        ibz_div(&x->coord[i], &r, &x->coord[i], &gcd);
#if RESP_TRIGGER_MODE == 7
        trigger_low();
#endif
    }
    response_profile.normalize_coord_rest_div_cycles = cycles_since(t);

    probe_ibz_limb_loop(&x->coord[0]);

    ibz_set(&sign, 2 * (0 > ibz_cmp(&ibz_const_zero, &x->denom)) - 1);
#if RESP_TRIGGER_MODE == 4
    trigger_high();
#endif
    t = response_cycles_read();
    ibz_mul(&x->coord[0], &sign, &x->coord[0]);
#if RESP_TRIGGER_MODE == 4
    trigger_low();
#endif
    for (int i = 1; i < 4; i++) {
        ibz_mul(&x->coord[i], &sign, &x->coord[i]);
    }
    ibz_mul(&x->denom, &sign, &x->denom);
    response_profile.normalize_sign_cycles = cycles_since(t);

    ibz_finalize(&gcd);
    ibz_finalize(&sign);
    ibz_finalize(&r);
    response_profile.normalize_total_cycles = cycles_since(t_total);
}

static int sample_from_ball_triggered(quat_alg_elem_t *res,
                                      const quat_lattice_t *lattice,
                                      const quat_alg_t *alg,
                                      const ibz_t *radius)
{
    uint32_t t_total = response_cycles_read();
    uint32_t t = t_total;
    ibz_vec_4_t box;
    ibz_vec_4_init(&box);
    ibz_mat_4x4_t U, G;
    ibz_mat_4x4_init(&U);
    ibz_mat_4x4_init(&G);
    ibz_vec_4_t x;
    ibz_vec_4_init(&x);
    ibz_t rad, tmp;
    ibz_init(&rad);
    ibz_init(&tmp);

    quat_lattice_gram(&G, lattice, alg);
    ibz_mul(&rad, radius, &lattice->denom);
    ibz_mul(&rad, &rad, &lattice->denom);
    ibz_mul(&rad, &rad, &ibz_const_two);

    int ok = quat_lattice_bound_parallelogram(&box, &U, &G, &rad);
    if (!ok) {
        goto err;
    }
    response_profile.sample_setup_cycles = cycles_since(t);

    t = response_cycles_read();
    response_profile.rejection_iters = 0;
    do {
        response_profile.rejection_iters++;
        for (int i = 0; i < 4; i++) {
            if (ibz_is_zero(&box[i])) {
                ibz_copy(&x[i], &ibz_const_zero);
            } else {
                ibz_add(&tmp, &box[i], &box[i]);
                ok &= ibz_rand_interval(&x[i], &ibz_const_zero, &tmp);
                ibz_sub(&x[i], &x[i], &box[i]);
                if (!ok) {
                    goto err;
                }
            }
        }
        ibz_mat_4x4_eval_t(&x, &x, &U);
        quat_qf_eval(&tmp, &G, &x);
    } while (ibz_is_zero(&tmp) || (ibz_cmp(&tmp, &rad) > 0));
    response_profile.sample_rejection_cycles = cycles_since(t);

#if RESP_TRIGGER_MODE == 2
    trigger_high();
#endif
    t = response_cycles_read();
#if RESP_TRIGGER_MODE == 7
    /* Per-coordinate materialization: coord[i] = sum_j basis[i][j]*x[j], one trigger
       pulse per coordinate (res->coord is not aliased with basis or x here). Together
       with the per-coord division pulses in normalize, one capture yields 8 aligned
       segments: materialize coord0..3, then divide coord0..3. */
    {
        ibz_t prod;
        ibz_init(&prod);
        for (int i = 0; i < 4; i++) {
            trigger_high();
            ibz_set(&res->coord[i], 0);
            for (int j = 0; j < 4; j++) {
                ibz_mul(&prod, &lattice->basis[i][j], &x[j]);
                ibz_add(&res->coord[i], &res->coord[i], &prod);
            }
            trigger_low();
        }
        ibz_finalize(&prod);
    }
#else
    ibz_mat_4x4_eval(&res->coord, &lattice->basis, &x);
#endif
    ibz_copy(&res->denom, &lattice->denom);
    normalize_response_triggered(res);
    response_profile.sample_finalize_cycles = cycles_since(t);
#if RESP_TRIGGER_MODE == 2
    trigger_low();
#endif

err:
    ibz_finalize(&rad);
    ibz_finalize(&tmp);
    ibz_vec_4_finalize(&x);
    ibz_mat_4x4_finalize(&U);
    ibz_mat_4x4_finalize(&G);
    ibz_vec_4_finalize(&box);
    response_profile.sample_total_cycles = cycles_since(t_total);
    return ok;
}

static void sample_response(quat_alg_elem_t *x,
                            const quat_lattice_t *lattice,
                            const ibz_t *lattice_content)
{
    ibz_t bound;
    ibz_init(&bound);
    ibz_pow(&bound, &ibz_const_two, SQIsign_response_length);
    ibz_sub(&bound, &bound, &ibz_const_one);
    ibz_mul(&bound, &bound, lattice_content);

    int ok = sample_from_ball_triggered(x, lattice, &QUATALG_PINFTY, &bound);
    (void)ok;

    ibz_finalize(&bound);
}

void sqisign_response_quat_element_core(quat_alg_elem_t *resp_quat,
                                        ibz_t *lattice_content,
                                        const quat_left_ideal_t *secret_ideal,
                                        const quat_left_ideal_t *challenge_ideal,
                                        const quat_left_ideal_t *commit_ideal)
{
    uint32_t t;
    quat_left_ideal_t challenge_secret;
    quat_lattice_t hom_challenge_to_commit;
    quat_lattice_t commit_conj;

    quat_left_ideal_init(&challenge_secret);
    quat_lattice_init(&hom_challenge_to_commit);
    quat_lattice_init(&commit_conj);

    t = response_cycles_read();
    quat_lideal_inter(&challenge_secret, challenge_ideal, secret_ideal, &QUATALG_PINFTY);
    response_profile.inter_cycles = cycles_since(t);
    t = response_cycles_read();
    quat_lattice_conjugate_without_hnf(&commit_conj, &commit_ideal->lattice);
    response_profile.conjugate_cycles = cycles_since(t);
    t = response_cycles_read();
    quat_lattice_intersect(&hom_challenge_to_commit, &challenge_secret.lattice, &commit_conj);
    response_profile.lattice_intersect_cycles = cycles_since(t);

    t = response_cycles_read();
    ibz_mul(lattice_content, &challenge_secret.norm, &commit_ideal->norm);
    response_profile.content_mul_cycles = cycles_since(t);
    sample_response(resp_quat, &hom_challenge_to_commit, lattice_content);

    quat_lattice_finalize(&commit_conj);
    quat_lattice_finalize(&hom_challenge_to_commit);
    quat_left_ideal_finalize(&challenge_secret);
}
