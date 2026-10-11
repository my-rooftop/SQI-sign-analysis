/* Instrumented SQIsign Algorithm 3.18 (SampleQuaternionResponse) for side-channel
 * measurement of line 2: alpha_sk,chl <- Conjugate(alpha_sk,chl).
 *
 * Spec:  sqisign-20260901.pdf, Algorithm 3.18 (p.21), Conjugate = Algorithm 4.38 (p.30)
 * Upstream: the-sqisign.git f417ebd, src/quaternion/ref/lvlx/protocol.c
 *             quat_response_element()   <- this file mirrors it line for line
 *
 * Builds both for the CW308_STM32F4 target (triggers active) and for the host
 * vector/KAT tool (triggers compiled out, CONJ_HOST_BUILD defined).
 */
#ifndef CW_SQISIGN_CONJ_CORE_H
#define CW_SQISIGN_CONJ_CORE_H

#include <stdint.h>

#include <quaternion.h>
#include <prng.h>

/* ---------------------------------------------------------------------------
 * Trigger modes (firmware build time: make CONJ_TRIG_MODE=n)
 *
 *   0  whole Algorithm 3.18
 *   1  line 1  EquivalentCoprimeIdeal  (quat_ideal_small_equivalent_coprime)
 *   2  line 2  Conjugate               -- primary target, one pulse
 *   3  line 2  Conjugate, one pulse per negated coordinate: coord[1], coord[2],
 *              coord[3]  -> 3 aligned segments per capture
 *   4  line 2  Conjugate, one pulse per ibz_t touched: denom copy, coord[0] copy,
 *              coord[1..3] negate  -> 5 aligned segments per capture
 *   5  line 5  IdealMultiplication     (quat_ideal_mul_O0)
 *   6  line 5  RandomBoundedElement    (quat_lattice_sample_from_ball)
 *   7  limb-walk probe over the pre-conjugate alpha_sk,chl coordinates; repeat
 *              count set by CONJ_PROBE_REPS.  Artificial: establishes the
 *              leakage ceiling and a profiling baseline, not a real code path.
 *   8  line 2  Conjugate, one pulse per limb of coord[1]  -> limb-index alignment
 * ------------------------------------------------------------------------- */
#ifndef CONJ_TRIG_MODE
#define CONJ_TRIG_MODE 2
#endif

#ifndef CONJ_PROBE_REPS
#define CONJ_PROBE_REPS 1
#endif

/* 0: Conjugate in place, exactly as spec line 2 and upstream protocol.c:37.
 *    With the 2026 mp module ibz_t is fixed-width two's complement, so ibz_neg
 *    rewrites every limb even in place -- this already leaks.
 * 1: Conjugate into a separate destination, then copy back. Destination limbs
 *    differ from the source, so the Hamming distance is against unrelated
 *    previous content instead of the limb's own old value.
 */
#ifndef CONJ_OUT_OF_PLACE
#define CONJ_OUT_OF_PLACE 0
#endif

typedef struct
{
    uint32_t total_cycles;         /* whole Algorithm 3.18 */
    uint32_t intersect_cycles;     /* Intersect, Alg 3.17, before line 1 */
    uint32_t equiv_coprime_cycles; /* line 1 */
    uint32_t conjugate_cycles;     /* line 2  -- the target */
    uint32_t ideal_mul_cycles;     /* line 5, IdealMultiplication */
    uint32_t sample_cycles;        /* line 5, RandomBoundedElement */
    uint32_t conj_coord_cycles[4]; /* per-coordinate cost inside line 2 */
    uint32_t conj_limbs[4];        /* NUM_LIMBS(bitlen) actually negated */
    uint32_t conj_bitlen[4];       /* bitlen of each coordinate at line 2 */
    uint32_t equiv_ok;
    uint32_t sample_ok;
    uint32_t formation_cycles;      /* line 5: eval + denom copy + normalize */
    uint32_t eval_cycles;           /* line 5: ibz_mat_4x4_eval alone */
    uint32_t normalize_cycles;      /* line 5: quat_alg_normalize */
    uint32_t resp_bitlen[4];        /* bitlen of alpha_rsp's coordinates */
} conj_profile_t;

extern volatile conj_profile_t conj_profile;

void conj_profile_reset(void);
void conj_cycles_init(void);
uint32_t conj_cycles_read(void);

/* Algorithm 3.18 SampleQuaternionResponse, instrumented.
 *
 * Mirrors upstream quat_response_element(). Arguments and semantics are
 * unchanged; sk_chall_quat receives alpha_sk,chl *after* line 2, and
 * sk_chall_quat_pre receives the value *before* line 2 (ground truth for
 * side-channel labels; pass NULL to skip).
 */
/* Algorithm 3.18 line 2 on its own, with the trigger placed per CONJ_TRIG_MODE.
 * `pre` receives the value entering the step (the side-channel label); pass NULL
 * to skip. Lets a campaign drive line 2 directly with realistic alpha values
 * instead of paying for a full Algorithm 3.18 run per trace -- the measured
 * window is identical, only the surrounding code differs. */
void conj_run_conjugate(quat_alg_elem_t *x, quat_alg_elem_t *pre);

/* ---------------------------------------------------------------------------
 * Algorithm 3.18 line 5, RandomBoundedElement.
 *
 * quat_lattice_sample_from_ball (lattice.c:159) rejection-samples a vector x and
 * then forms the response quaternion from it:
 *
 *     ibz_mat_4x4_eval(&res->coord, &lat.basis, &x);   // coord[i] = sum_j b[i][j]*x[j]
 *     ibz_copy(&res->denom, &lat.denom);
 *     quat_alg_normalize(res);
 *
 * The rejection loop runs a data-dependent number of times, so the whole
 * subroutine does not align trace to trace. The formation runs exactly once
 * after acceptance, so triggering it alone gives the same deterministic
 * alignment line 2 had. Unlike Conjugate, alpha_rsp is *formed* by a
 * multiply-accumulate, which is the structure that leaks in Section 5.2.
 *
 * Trigger modes 15-17 (continuing conj_core's table):
 *   15  the whole formation: eval, denom copy, normalize
 *   16  ibz_mat_4x4_eval only -- the multiply-accumulate
 *   17  one pulse per output coordinate inside the eval -> 4 aligned segments
 *       (each ~9.7k cycles; needs a long capture window)
 *   18  one pulse per limb-wise write-out of a finished coordinate -> 4 short
 *       aligned segments carrying the complete formed value
 * ------------------------------------------------------------------------- */
void conj_run_response_formation(quat_alg_elem_t *res,
                                 const quat_lattice_t *lat,
                                 const ibz_vec_4_t *x,
                                 int radius_bound);

/* Faithful copy of quat_lattice_sample_from_ball that also returns the accepted
 * sample vector, so a campaign can replay the formation on the board with the
 * exact x the real routine would have used. Host build only. */
int conj_sample_from_ball_capture(quat_alg_elem_t *res,
                                  ibz_t *res_ideal_norm,
                                  ibz_vec_4_t *x_out,
                                  const quat_lattice_t *lattice,
                                  const ibz_t *radius,
                                  prng_domain_ctx_t *prng_domain);

void conj_response_element(quat_ideal_t *sk_chall_ideal,
                          quat_alg_elem_t *resp_quat,
                          ibz_t *norm,
                          quat_alg_elem_t *sk_chall_quat,
                          quat_alg_elem_t *sk_chall_quat_pre,
                          const quat_ideal_t *sk_ideal,
                          const quat_ideal_t *ideal_chall_two,
                          const quat_alg_elem_t *chall_split,
                          const quat_ideal_t *ideal_commit,
                          prng_domain_ctx_t *prng_domain,
                          int e);

#endif /* CW_SQISIGN_CONJ_CORE_H */
