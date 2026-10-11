/* Host-side vector generator and KAT checker for the Algorithm 3.18 line 2
 * Conjugate side-channel experiment.
 *
 * Builds the same quaternion inputs the signing flow feeds to
 * quat_response_element (Alg 3.18), runs the instrumented core, and reports the
 * value of alpha_sk,chl on both sides of line 2 together with its raw limbs --
 * the ground truth the trace analysis is scored against.
 *
 * Input generation is quaternion-only (no EC / isogeny code), following upstream
 * src/quaternion/ref/lvlx/test/protocol.c::quat_test_protocol_response_element_loop.
 *
 * Commands:
 *   gen  [seed_hex]            emit a board vector on stdout
 *   run  <vector_file>         run Alg 3.18 on a vector, print ground truth + KAT
 *   selftest [n]               generate+run n vectors, check upstream invariants
 *   alphas <n> [seed]          emit n real alpha_sk,chl values entering line 2
 *   responses <n> [seed]       emit n line 5 records (lattice, sample vector, alpha_rsp)
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#include "conj_core.h"

#include <internal.h>
#include <mp_internal.h>
#include <quaternion_constants.h>
#include <quaternion_data.h>
#include <torsion_constants.h>
#include <encoded_sizes.h>
#include <prng.h>
#include <quaternion_tests.h>
#include <id2iso.h>
#include <torsion_constants.h>

#define IBZ_STR_MAX (IBZ_NLIMBS * NUM_BITS_LIMB + 4)
#define NUM_LIMBS_OF(bits) ((int)(((bits) + (NUM_BITS_LIMB - 1)) / (NUM_BITS_LIMB)))

/* ------------------------------------------------------- deterministic rng */

/* The vendored pqm4 source set does not ship a randombytes(); pqm4 supplies one
 * on the target. Here we provide a deterministic one so host vectors and board
 * runs are reproducible and comparable. Same generator as cw_randombytes.c. */
static uint32_t rng_state = 0x13579bdfu;

int
randombytes(unsigned char *x, unsigned long long xlen)
{
    for (unsigned long long i = 0; i < xlen; i++) {
        rng_state = rng_state * 1664525u + 1013904223u;
        x[i] = (unsigned char)(rng_state >> 24);
    }
    return 0;
}

static void
rng_seed_from_hex(const char *hex)
{
    rng_state = 0x13579bdfu;
    if (hex == NULL) {
        return;
    }
    for (const char *p = hex; *p; p++) {
        rng_state = rng_state * 1664525u + 1013904223u + (uint32_t)(unsigned char)*p;
    }
}

/* ------------------------------------------------------------ text helpers */

static void
print_ibz(const char *label, const ibz_t *x)
{
    char buf[IBZ_STR_MAX];
    if (!ibz_convert_to_str(x, buf, 16)) {
        strcpy(buf, "?");
    }
    printf("%s %s\n", label, buf);
}

static void
print_ibz_bare(const ibz_t *x)
{
    char buf[IBZ_STR_MAX];
    if (!ibz_convert_to_str(x, buf, 16)) {
        strcpy(buf, "?");
    }
    printf("%s\n", buf);
}

/* Raw limbs exactly as ibz_neg walks them: index 0 first, NUM_LIMBS(bitlen)
 * of them. These are the 32-bit words whose Hamming weight the trace leaks. */
static void
print_limbs(const char *label, const ibz_t *x)
{
    int n = NUM_LIMBS_OF(x->bitlen);
    printf("%s bitlen=%d limbs=%d", label, x->bitlen, n);
    for (int i = 0; i < n; i++) {
        printf(" %08lx", (unsigned long)x->limbs[i]);
    }
    printf("\n");
}

static int
read_ibz_token(FILE *f, ibz_t *x)
{
    char tok[IBZ_STR_MAX];
    if (fscanf(f, "%2047s", tok) != 1) {
        return 0;
    }
    return ibz_set_from_str(x, tok, 16);
}

/* ------------------------------------------------------- vector generation */

typedef struct
{
    quat_ideal_t sk_ideal;
    quat_ideal_t ideal_chall;
    quat_alg_elem_t chall_split;
    quat_ideal_t ideal_commit;
} conj_vector_t;

static void
vector_init(conj_vector_t *v)
{
    quat_ideal_init(&v->sk_ideal);
    quat_ideal_init(&v->ideal_chall);
    quat_ideal_init(&v->ideal_commit);
    quat_alg_elem_init(&v->chall_split);
}

/* Mirrors the signing flow's quaternion inputs:
 *   sk_ideal     keygen.c:32-37  random prime-norm O0 ideal, reduced to a small
 *                                equivalent of prime norm
 *   ideal_commit sign.c:19-24    random prime-norm O0 ideal, reduced to a small
 *                                equivalent of odd norm
 *   ideal_chall  a 2^CHALLENGE_BITS-norm O0 ideal, with chall_split its split part
 */
static int
vector_generate(conj_vector_t *v, prng_domain_ctx_t *prng_key, prng_domain_ctx_t *prng_com)
{
    ibz_t p, n, tmp;
    ibz_init(&p);
    ibz_init(&n);
    ibz_init(&tmp);
    int ok = 1;

    /* Secret key ideal, exactly as keygen.c:32-37 builds it: a random O0-ideal of
     * norm SEC_DEGREE, replaced by a small equivalent of prime norm.
     *
     * An earlier version drew a random prime of bitsize(p)/2 for the norm, the
     * way upstream's quaternion test does. That reached a divide-by-zero inside
     * quat_ideal_small_equivalent_coprime_enumeration -> quat_ideal_create_O0_odd
     * -> ibz_mod -> mp_divmod_schoolbook (mp.c:581) on some draws, and it is not
     * what a real signer does anyway. */
    (void)p;
    ok = ok && quat_random_ideal_O0_given_prime_norm(&v->sk_ideal, &SEC_DEGREE, prng_key);
    ok = ok && quat_ideal_small_equivalent_coprime(NULL, &v->sk_ideal, &v->sk_ideal, &ibz_const_zero, prng_key);

    /* commitment ideal */
    ok = ok && quat_random_ideal_O0_given_prime_norm(&v->ideal_commit, &COM_DEGREE, prng_com);
    ok = ok && quat_ideal_small_equivalent_coprime(NULL, &v->ideal_commit, &v->ideal_commit, &ibz_const_two, prng_com);

    /* Challenge ideal and its split part, via the exact call the signing flow
     * makes (sign.c::compute_challenge_ideal_signature, line 48):
     *
     *   id2iso_kernel_dlogs_to_ideal_even(ideal_chall_two, chall_split, vec, CHALLENGE_BITS)
     *
     * vec = Msk * [1, chall_coeff]^T in the real flow. Here we sample vec
     * directly -- same distribution shape, and it needs no keygen or curve
     * arithmetic. This path is what guarantees chall_split.denom == 1, which
     * quat_ideal_intersect_O0 requires (ideal.c:435); the test-only
     * quat_represent_integer_even() route yields denom 2 and is rejected there.
     */
    if (ok) {
        ibz_vec_2_t vec;
        ibz_vec_2_init(&vec);
        ibz_mul_2exp(&n, &ibz_const_one, CHALLENGE_BITS);

        int guard = 0;
        int built = 0;
        while (!built && guard++ < 1000) {
            ibz_set(&vec.v[0], 1, 2);
            if (!ibz_rand_interval(&vec.v[1], &ibz_const_zero, &n)) {
                continue;
            }
            id2iso_kernel_dlogs_to_ideal_even(&v->ideal_chall, &v->chall_split, &vec, CHALLENGE_BITS);
            built = ibz_is_one(&v->chall_split.denom) && !quat_alg_elem_is_zero(&v->chall_split) &&
                    !ibz_is_zero(&v->ideal_chall.norm);
        }
        ok = ok && built;
    }

    return ok;
}

/* The 14 values a vector carries, in wire order. */
static void
vector_fields(const conj_vector_t *v, const ibz_t *out[14])
{
    int k = 0;
    out[k++] = &v->sk_ideal.x;
    out[k++] = &v->sk_ideal.y;
    out[k++] = &v->sk_ideal.norm;
    out[k++] = &v->ideal_chall.x;
    out[k++] = &v->ideal_chall.y;
    out[k++] = &v->ideal_chall.norm;
    out[k++] = &v->chall_split.denom;
    for (int i = 0; i < 4; i++) {
        out[k++] = &v->chall_split.coord.v[i];
    }
    out[k++] = &v->ideal_commit.x;
    out[k++] = &v->ideal_commit.y;
    out[k++] = &v->ideal_commit.norm;
}

static void
vector_emit(const conj_vector_t *v, const char *seed_hex)
{
    printf("# SQIsign Alg 3.18 (SampleQuaternionResponse) board vector\n");
    printf("# variant p324_3  RADIX=32  IBZ_NLIMBS=%d  CHALLENGE_BITS=%d  RESPONSE_BITS=%d\n",
           IBZ_NLIMBS,
           CHALLENGE_BITS,
           RESPONSE_BITS);
    printf("# field order: sk.x sk.y sk.norm chall.x chall.y chall.norm "
           "split.denom split.c0..c3 com.x com.y com.norm\n");
    printf("seed %s\n", seed_hex ? seed_hex : "default");

    /* ibz_t carries a declared bit bound, and every downstream ibz_t inherits it.
     * ibz_set_from_str() re-derives the bound from the parsed magnitude
     * (mp.c:311), which is tighter than the bound these values had in memory --
     * and that difference propagates all the way to alpha_sk,chl, changing how
     * many limbs ibz_neg walks at line 2. So the bounds travel with the values
     * and the reader restores them with ibz_set_bound(). Without this the board
     * measures a narrower window than a real signer would execute. */
    {
        const ibz_t *f[14];
        vector_fields(v, f);
        printf("bounds");
        for (int i = 0; i < 14; i++) {
            printf(" %d", f[i]->bitlen);
        }
        printf("\n");
    }

    printf("vector\n");
    print_ibz_bare(&v->sk_ideal.x);
    print_ibz_bare(&v->sk_ideal.y);
    print_ibz_bare(&v->sk_ideal.norm);
    print_ibz_bare(&v->ideal_chall.x);
    print_ibz_bare(&v->ideal_chall.y);
    print_ibz_bare(&v->ideal_chall.norm);
    print_ibz_bare(&v->chall_split.denom);
    for (int i = 0; i < 4; i++) {
        print_ibz_bare(&v->chall_split.coord.v[i]);
    }
    print_ibz_bare(&v->ideal_commit.x);
    print_ibz_bare(&v->ideal_commit.y);
    print_ibz_bare(&v->ideal_commit.norm);
}

static int
vector_load(conj_vector_t *v, const char *path, char *seed_out, size_t seed_len)
{
    FILE *f = fopen(path, "r");
    if (!f) {
        perror(path);
        return 0;
    }

    /* skip comments, pick up "seed <hex>" and "bounds ...", stop after "vector" */
    char line[4096];
    seed_out[0] = '\0';
    int saw_vector = 0;
    int bounds[14];
    int have_bounds = 0;
    long pos = 0;
    while (fgets(line, sizeof(line), f)) {
        if (line[0] == '#') {
            pos = ftell(f);
            continue;
        }
        if (strncmp(line, "bounds ", 7) == 0) {
            const char *q = line + 7;
            have_bounds = 1;
            for (int i = 0; i < 14; i++) {
                char *end = NULL;
                long b = strtol(q, &end, 10);
                if (end == q) {
                    have_bounds = 0;
                    break;
                }
                bounds[i] = (int)b;
                q = end;
            }
            pos = ftell(f);
            continue;
        }
        if (strncmp(line, "seed ", 5) == 0) {
            char *nl = strchr(line, '\n');
            if (nl) {
                *nl = '\0';
            }
            snprintf(seed_out, seed_len, "%s", line + 5);
            pos = ftell(f);
            continue;
        }
        if (strncmp(line, "vector", 6) == 0) {
            saw_vector = 1;
            pos = ftell(f);
            break;
        }
        pos = ftell(f);
    }
    if (!saw_vector) {
        fprintf(stderr, "%s: no 'vector' section\n", path);
        fclose(f);
        return 0;
    }
    fseek(f, pos, SEEK_SET);

    int ok = 1;
    ok = ok && read_ibz_token(f, &v->sk_ideal.x);
    ok = ok && read_ibz_token(f, &v->sk_ideal.y);
    ok = ok && read_ibz_token(f, &v->sk_ideal.norm);
    ok = ok && read_ibz_token(f, &v->ideal_chall.x);
    ok = ok && read_ibz_token(f, &v->ideal_chall.y);
    ok = ok && read_ibz_token(f, &v->ideal_chall.norm);
    ok = ok && read_ibz_token(f, &v->chall_split.denom);
    for (int i = 0; i < 4; i++) {
        ok = ok && read_ibz_token(f, &v->chall_split.coord.v[i]);
    }
    ok = ok && read_ibz_token(f, &v->ideal_commit.x);
    ok = ok && read_ibz_token(f, &v->ideal_commit.y);
    ok = ok && read_ibz_token(f, &v->ideal_commit.norm);
    fclose(f);
    if (!ok) {
        fprintf(stderr, "%s: malformed vector\n", path);
        return 0;
    }

    if (have_bounds) {
        const ibz_t *cf[14];
        vector_fields(v, cf);
        for (int i = 0; i < 14; i++) {
            ibz_t *x = (ibz_t *)cf[i];
            if (bounds[i] > x->bitlen) {
                ibz_set_bound(x, bounds[i]);
            }
        }
    } else {
        fprintf(stderr, "%s: no 'bounds' line; alpha_sk,chl will get a tighter "
                        "bit bound than a real signer produces\n",
                path);
    }
    return ok;
}

/* ------------------------------------------------------------------ runner */

typedef struct
{
    quat_ideal_t sk_chall_ideal;
    quat_alg_elem_t resp_quat;
    quat_alg_elem_t sk_chall_quat;     /* alpha_sk,chl after line 2 */
    quat_alg_elem_t sk_chall_quat_pre; /* alpha_sk,chl before line 2 */
    ibz_t norm;
} conj_result_t;

static void
result_init(conj_result_t *r)
{
    quat_ideal_init(&r->sk_chall_ideal);
    quat_alg_elem_init(&r->resp_quat);
    quat_alg_elem_init(&r->sk_chall_quat);
    quat_alg_elem_init(&r->sk_chall_quat_pre);
    ibz_init(&r->norm);
}

static void
run_alg318(conj_result_t *r, const conj_vector_t *v, prng_domain_ctx_t *prng_res)
{
    conj_profile_reset();
    conj_response_element(&r->sk_chall_ideal,
                          &r->resp_quat,
                          &r->norm,
                          &r->sk_chall_quat,
                          &r->sk_chall_quat_pre,
                          &v->sk_ideal,
                          &v->ideal_chall,
                          &v->chall_split,
                          &v->ideal_commit,
                          prng_res,
                          RESPONSE_BITS - 1);
}

/* Upstream invariants from quat_test_protocol_response_element_loop:
 *   Norm(resp_quat) == norm * n(sk_chall_ideal) * n(ideal_commit)
 *   bitsize(norm) <= e
 *   conj(alpha_sk,chl) == the value that entered line 2
 */
static int
check_invariants(const conj_result_t *r, const conj_vector_t *v)
{
    int bad = 0;
    ibz_t num, den, want;
    ibz_init(&num);
    ibz_init(&den);
    ibz_init(&want);

    if (ibz_bitsize(&r->norm) > RESPONSE_BITS - 1) {
        fprintf(stderr, "FAIL: bitsize(norm)=%d > e=%d\n", ibz_bitsize(&r->norm), RESPONSE_BITS - 1);
        bad = 1;
    }

    quat_alg_norm(&num, &den, &r->resp_quat, &QUATALG_PINFTY);
    if (!ibz_is_one(&den)) {
        fprintf(stderr, "FAIL: Norm(resp_quat) is not integral\n");
        bad = 1;
    }
    ibz_mul(&want, &r->norm, &r->sk_chall_ideal.norm);
    ibz_mul(&want, &want, &v->ideal_commit.norm);
    if (ibz_cmp(&want, &num) != 0) {
        fprintf(stderr, "FAIL: Norm(resp_quat) != norm * n(I'sk,chl) * n(Icom)\n");
        bad = 1;
    }

    /* line 2 is an involution: conjugating the result must restore the input */
    {
        quat_alg_elem_t back;
        quat_alg_elem_init(&back);
        quat_alg_conj(&back, &r->sk_chall_quat);
        for (int i = 0; i < 4; i++) {
            if (ibz_cmp(&back.coord.v[i], &r->sk_chall_quat_pre.coord.v[i]) != 0) {
                fprintf(stderr, "FAIL: conj(conj(alpha)) != alpha at coord[%d]\n", i);
                bad = 1;
            }
        }
    }

    /* the GCD coprimality check of Alg 3.18 line 3 (hoisted into sign.c upstream) */
    {
        ibz_t g;
        ibz_init(&g);
        ibz_gcd(&g, &r->sk_chall_ideal.norm, &v->ideal_commit.norm);
        if (!ibz_is_one(&g)) {
            fprintf(stderr, "FAIL: Alg 3.18 line 3, norms not coprime\n");
            bad = 1;
        }
    }

    return bad;
}

static void
report(const conj_result_t *r)
{
    printf("\n## Algorithm 3.18 line 2 target: alpha_sk,chl <- Conjugate(alpha_sk,chl)\n");
    printf("# value entering line 2 (side-channel ground truth)\n");
    print_ibz("pre.denom", &r->sk_chall_quat_pre.denom);
    for (int i = 0; i < 4; i++) {
        char lbl[24];
        snprintf(lbl, sizeof(lbl), "pre.coord[%d]", i);
        print_ibz(lbl, &r->sk_chall_quat_pre.coord.v[i]);
    }
    printf("\n# value leaving line 2\n");
    print_ibz("post.denom", &r->sk_chall_quat.denom);
    for (int i = 0; i < 4; i++) {
        char lbl[24];
        snprintf(lbl, sizeof(lbl), "post.coord[%d]", i);
        print_ibz(lbl, &r->sk_chall_quat.coord.v[i]);
    }

    printf("\n# recovery target size\n");
    {
        int total_limbs = 0, total_info = 0, total_mag = 0;
        for (int i = 1; i < 4; i++) {
            const ibz_t *x = &r->sk_chall_quat_pre.coord.v[i];
            int nl = NUM_LIMBS_OF(x->bitlen);
            int mag = ibz_bitsize(x);
            digit_t ext = (x->limbs[nl - 1] >> (NUM_BITS_LIMB - 1)) ? (digit_t)-1 : (digit_t)0;
            int sign_limbs = 0;
            for (int k = nl - 1; k >= 0 && x->limbs[k] == ext; k--) {
                sign_limbs++;
            }
            printf("coord[%d]  bound=%d bits  limbs=%d  magnitude=%d bits  "
                   "sign-ext limbs=%d  informative limbs=%d\n",
                   i,
                   x->bitlen,
                   nl,
                   mag,
                   sign_limbs,
                   nl - sign_limbs);
            total_limbs += nl;
            total_info += nl - sign_limbs;
            total_mag += mag;
        }
        printf("TOTAL     negated limbs=%d (%d bits of declared width)  "
               "informative limbs=%d (%d bits)  magnitude sum=%d bits\n",
               total_limbs,
               total_limbs * NUM_BITS_LIMB,
               total_info,
               total_info * NUM_BITS_LIMB,
               total_mag);
    }

    printf("\n# limbs walked by ibz_neg: coord[1..3] only (coord[0] and denom are copied)\n");
    for (int i = 1; i < 4; i++) {
        char lbl[24];
        snprintf(lbl, sizeof(lbl), "pre.limbs[%d]", i);
        print_limbs(lbl, &r->sk_chall_quat_pre.coord.v[i]);
        snprintf(lbl, sizeof(lbl), "post.limbs[%d]", i);
        print_limbs(lbl, &r->sk_chall_quat.coord.v[i]);
    }

    printf("\n# result of Algorithm 3.18\n");
    print_ibz("resp.norm", &r->norm);
    print_ibz("resp.denom", &r->resp_quat.denom);
    for (int i = 0; i < 4; i++) {
        char lbl[24];
        snprintf(lbl, sizeof(lbl), "resp.coord[%d]", i);
        print_ibz(lbl, &r->resp_quat.coord.v[i]);
    }
    print_ibz("skchall.norm", &r->sk_chall_ideal.norm);

    printf("\n# profile (cycle fields are zero on the host build)\n");
    printf("conj_bitlen  %u %u %u %u\n",
           conj_profile.conj_bitlen[0],
           conj_profile.conj_bitlen[1],
           conj_profile.conj_bitlen[2],
           conj_profile.conj_bitlen[3]);
    printf("conj_limbs   %u %u %u %u\n",
           conj_profile.conj_limbs[0],
           conj_profile.conj_limbs[1],
           conj_profile.conj_limbs[2],
           conj_profile.conj_limbs[3]);
    printf("equiv_ok=%u sample_ok=%u\n", conj_profile.equiv_ok, conj_profile.sample_ok);
}

/* -------------------------------------------------------------------- main */

static int
cmd_gen(const char *seed_hex)
{
    rng_seed_from_hex(seed_hex);
    if (prng_seed() != 0) {
        fprintf(stderr, "prng_seed failed\n");
        return 1;
    }
    prng_domain_ctx_t prng_key, prng_com;
    prng_domain_seed(&prng_key, DS_KEYGEN_DOMAIN);
    prng_domain_seed(&prng_com, DS_COMMITMENT_DOMAIN);

    conj_vector_t v;
    vector_init(&v);
    if (!vector_generate(&v, &prng_key, &prng_com)) {
        fprintf(stderr, "vector generation failed\n");
        return 1;
    }
    vector_emit(&v, seed_hex);
    return 0;
}

static int
cmd_run(const char *path)
{
    char seed[256];
    conj_vector_t v;
    vector_init(&v);
    if (!vector_load(&v, path, seed, sizeof(seed))) {
        return 1;
    }
    rng_seed_from_hex(seed[0] ? seed : NULL);
    if (prng_seed() != 0) {
        fprintf(stderr, "prng_seed failed\n");
        return 1;
    }
    prng_domain_ctx_t prng_res;
    prng_domain_seed(&prng_res, DS_RESPONSE_DOMAIN);

    conj_result_t r;
    result_init(&r);
    run_alg318(&r, &v, &prng_res);
    report(&r);
    int bad = check_invariants(&r, &v);
    printf("\n%s\n", bad ? "KAT: FAIL" : "KAT: OK");
    return bad;
}

/* Emit n real alpha_sk,chl values, each as it stands entering Algorithm 3.18
 * line 2. Running the whole algorithm on the board costs ~13 s per trace, which
 * makes a statistically useful MSB-recovery campaign impractical; generating the
 * values here and driving line 2 directly on the board costs milliseconds and
 * measures the identical window. */
static int
cmd_alphas(int n, const char *seed_hex)
{
    printf("# alpha_sk,chl values entering Algorithm 3.18 line 2\n");
    printf("# per record: 'a' <bitlen x5> then denom, coord0..3 as hex\n");
    printf("count %d\n", n);
    for (int i = 0; i < n; i++) {
        char seed[96];
        snprintf(seed, sizeof(seed), "%s_%d", seed_hex ? seed_hex : "alpha", i);
        rng_seed_from_hex(seed);
        if (prng_seed() != 0) {
            fprintf(stderr, "prng_seed failed at %d\n", i);
            return 1;
        }
        prng_domain_ctx_t prng_key, prng_com, prng_res;
        prng_domain_seed(&prng_key, DS_KEYGEN_DOMAIN);
        prng_domain_seed(&prng_com, DS_COMMITMENT_DOMAIN);
        prng_domain_seed(&prng_res, DS_RESPONSE_DOMAIN);

        conj_vector_t v;
        vector_init(&v);
        if (!vector_generate(&v, &prng_key, &prng_com)) {
            fprintf(stderr, "vector generation failed at %d\n", i);
            return 1;
        }
        conj_result_t r;
        result_init(&r);
        run_alg318(&r, &v, &prng_res);

        const quat_alg_elem_t *q = &r.sk_chall_quat_pre;
        printf("a %d %d %d %d %d\n", q->denom.bitlen, q->coord.v[0].bitlen,
               q->coord.v[1].bitlen, q->coord.v[2].bitlen, q->coord.v[3].bitlen);
        print_ibz_bare(&q->denom);
        for (int k = 0; k < 4; k++) {
            print_ibz_bare(&q->coord.v[k]);
        }
    }
    return 0;
}

/* Emit n line 5 records. Each is the product lattice, the accepted sample
 * vector and the radius bound that quat_lattice_sample_from_ball used, plus
 * alpha_rsp for scoring.
 *
 * The rejection loop runs a data-dependent number of times and costs seconds on
 * the board; the formation after it runs once. Capturing the formation with the
 * real x is both faster and better aligned than tracing the whole subroutine. */
static int
cmd_responses(int n, const char *seed_hex)
{
    printf("# Algorithm 3.18 line 5: product lattice, accepted sample vector, alpha_rsp\n");
    printf("# per record: 'r' <radius_bound> <bitlen x21>, then denom, basis[0..3][0..3],\n");
    printf("#             x[0..3], then (for scoring) alpha_rsp denom and coord0..3\n");
    printf("count %d\n", n);

    for (int i = 0; i < n; i++) {
        char seed[96];
        snprintf(seed, sizeof(seed), "%s_%d", seed_hex ? seed_hex : "resp", i);
        rng_seed_from_hex(seed);
        if (prng_seed() != 0) {
            fprintf(stderr, "prng_seed failed at %d\n", i);
            return 1;
        }
        prng_domain_ctx_t prng_key, prng_com, prng_res;
        prng_domain_seed(&prng_key, DS_KEYGEN_DOMAIN);
        prng_domain_seed(&prng_com, DS_COMMITMENT_DOMAIN);
        prng_domain_seed(&prng_res, DS_RESPONSE_DOMAIN);

        conj_vector_t v;
        vector_init(&v);
        if (!vector_generate(&v, &prng_key, &prng_com)) {
            fprintf(stderr, "vector generation failed at %d\n", i);
            return 1;
        }

        /* Reproduce Algorithm 3.18 up to the point line 5 samples. */
        quat_ideal_t sk_chall_ideal;
        quat_alg_elem_t gen;
        quat_lattice_t prod;
        ibz_t bound, lattice_content;
        quat_ideal_init(&sk_chall_ideal);
        quat_alg_elem_init(&gen);
        quat_lattice_init(&prod);
        ibz_init(&bound);
        ibz_init(&lattice_content);

        quat_ideal_intersect_O0(&sk_chall_ideal, &v.ideal_chall, &v.sk_ideal, &v.chall_split);
        if (!quat_ideal_small_equivalent_coprime(
                &gen, &sk_chall_ideal, &sk_chall_ideal, &ibz_const_zero, &PRNG_default_domain)) {
            fprintf(stderr, "EquivalentCoprimeIdeal failed at %d\n", i);
            return 1;
        }
        quat_alg_conj(&gen, &gen);
        quat_ideal_mul_O0(&prod, &sk_chall_ideal, &v.ideal_commit);
        ibz_div_2exp(&lattice_content, &prod.basis.m[0][0], 1);

        int e = RESPONSE_BITS - 1;
        ibz_mul_2exp(&bound, &ibz_const_one, e);
        ibz_sub(&bound, &bound, &ibz_const_one);
        ibz_set_bound(&bound, e + 1);

        quat_alg_elem_t rsp;
        ibz_vec_4_t xs;
        ibz_t norm;
        quat_alg_elem_init(&rsp);
        ibz_vec_4_init(&xs);
        ibz_init(&norm);
        if (!conj_sample_from_ball_capture(&rsp, &norm, &xs, &prod, &bound, &prng_res)) {
            fprintf(stderr, "sampling failed at %d\n", i);
            return 1;
        }

        printf("r %d %d", ibz_get_bound(&bound), prod.denom.bitlen);
        for (int a = 0; a < 4; a++) {
            for (int b = 0; b < 4; b++) {
                printf(" %d", prod.basis.m[a][b].bitlen);
            }
        }
        for (int a = 0; a < 4; a++) {
            printf(" %d", xs.v[a].bitlen);
        }
        printf("\n");
        print_ibz_bare(&prod.denom);
        for (int a = 0; a < 4; a++) {
            for (int b = 0; b < 4; b++) {
                print_ibz_bare(&prod.basis.m[a][b]);
            }
        }
        for (int a = 0; a < 4; a++) {
            print_ibz_bare(&xs.v[a]);
        }
        print_ibz_bare(&rsp.denom);
        for (int a = 0; a < 4; a++) {
            print_ibz_bare(&rsp.coord.v[a]);
        }
    }
    return 0;
}

static int
cmd_selftest(int n)
{
    int bad = 0;
    for (int i = 0; i < n; i++) {
        char seed[64];
        snprintf(seed, sizeof(seed), "selftest%d", i);
        rng_seed_from_hex(seed);
        if (prng_seed() != 0) {
            fprintf(stderr, "prng_seed failed\n");
            return 1;
        }
        prng_domain_ctx_t prng_key, prng_com, prng_res;
        prng_domain_seed(&prng_key, DS_KEYGEN_DOMAIN);
        prng_domain_seed(&prng_com, DS_COMMITMENT_DOMAIN);
        prng_domain_seed(&prng_res, DS_RESPONSE_DOMAIN);

        conj_vector_t v;
        vector_init(&v);
        if (!vector_generate(&v, &prng_key, &prng_com)) {
            fprintf(stderr, "[%d] vector generation failed\n", i);
            bad = 1;
            continue;
        }
        conj_result_t r;
        result_init(&r);
        run_alg318(&r, &v, &prng_res);
        int b = check_invariants(&r, &v);
        printf("[%d] %s\n", i, b ? "FAIL" : "ok");
        for (int c = 1; c < 4; c++) {
            const ibz_t *x = &r.sk_chall_quat_pre.coord.v[c];
            int nl = NUM_LIMBS_OF(x->bitlen);
            int mag = ibz_bitsize(x);
            /* count top limbs that are pure sign extension (0x00000000 / 0xffffffff) */
            digit_t ext = (x->limbs[nl - 1] >> (NUM_BITS_LIMB - 1)) ? (digit_t)-1 : (digit_t)0;
            int sign_limbs = 0;
            for (int k = nl - 1; k >= 0 && x->limbs[k] == ext; k--) {
                sign_limbs++;
            }
            printf("      coord[%d] bound=%d limbs=%d magnitude=%d bits  "
                   "sign-ext limbs=%d  informative limbs=%d\n",
                   c,
                   x->bitlen,
                   nl,
                   mag,
                   sign_limbs,
                   nl - sign_limbs);
        }
        bad |= b;
    }
    printf("\n%s\n", bad ? "selftest: FAIL" : "selftest: OK");
    return bad;
}

int
main(int argc, char **argv)
{
    if (argc < 2) {
        fprintf(stderr,
                "usage: %s gen [seed]\n"
                "       %s run <vector_file>\n"
                "       %s selftest [n]\n",
                argv[0],
                argv[0],
                argv[0]);
        return 2;
    }
    if (strcmp(argv[1], "gen") == 0) {
        return cmd_gen(argc > 2 ? argv[2] : NULL);
    }
    if (strcmp(argv[1], "run") == 0) {
        if (argc < 3) {
            fprintf(stderr, "run needs a vector file\n");
            return 2;
        }
        return cmd_run(argv[2]);
    }
    if (strcmp(argv[1], "alphas") == 0) {
        if (argc < 3) {
            fprintf(stderr, "alphas needs a count\n");
            return 2;
        }
        return cmd_alphas(atoi(argv[2]), argc > 3 ? argv[3] : NULL);
    }
    if (strcmp(argv[1], "responses") == 0) {
        if (argc < 3) {
            fprintf(stderr, "responses needs a count\n");
            return 2;
        }
        return cmd_responses(atoi(argv[2]), argc > 3 ? argv[3] : NULL);
    }
    if (strcmp(argv[1], "selftest") == 0) {
        return cmd_selftest(argc > 2 ? atoi(argv[2]) : 5);
    }
    fprintf(stderr, "unknown command: %s\n", argv[1]);
    return 2;
}
