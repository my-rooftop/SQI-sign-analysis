/* CW308_STM32F4 / STM32F415 simpleserial target for the SQIsign Algorithm 3.18
 * line 2 Conjugate side-channel experiment.
 *
 * The host sends the quaternion inputs of Alg 3.18 as ASCII hex (they are built
 * on the PC by conj_vector_host, which needs the EC/isogeny code the board does
 * not carry), the board runs the instrumented Algorithm 3.18 with the trigger
 * placed per CONJ_TRIG_MODE, and returns alpha_sk,chl as it stood on both sides
 * of line 2 plus the cycle profile.
 *
 * Commands (SS_VER_1_1, all reply 'r'):
 *   i  ->  info: struct sizes, trigger mode, out-of-place flag, IBZ_NLIMBS
 *   x  ->  clear input and output buffers
 *   a  ->  append CONJ_CHUNK bytes of vector text to the input buffer
 *   c  ->  parse the vector and run Algorithm 3.18 (trigger fires here)
 *   n  ->  seed the board RNG from 16 bytes
 *   l  ->  output length, big endian uint16
 *   g  ->  read output chunk at index (big endian uint16)
 *   p  ->  write the cycle/limb profile into the output buffer
 *
 * Algorithm 3.18 line 2 fast path (alpha injected, no full run):
 *   u  ->  load one alpha_sk,chl
 *   z  ->  run line 2 on it (trigger fires here)
 *   q  ->  read the alpha that entered line 2, for scoring
 *
 * Algorithm 3.18 line 5 fast path (lattice and sample vector injected):
 *   b  ->  load the product lattice, the sample vector and the radius bound
 *   m  ->  form alpha_rsp from them (trigger fires here)
 *   o  ->  read alpha_rsp back, for scoring
 *
 * 'v', 'w' and 'y' are reserved by simpleserial_init and must not be reused.
 *
 * Algorithm 3.23 line 12 target (secret matrix times known challenge coefficient):
 *   k  ->  load the fixed secret matrix Msk from the ASCII buffer
 *   e  ->  run the matrix-vector product on CHALLENGE_BYTES of challenge coefficient
 *   f  ->  read the result and the cycle profile of that run
 */
#include "hal.h"
#include "simpleserial.h"

#include <stddef.h>
#include <stdint.h>
#include <string.h>

#include "conj_core.h"
#include "chall_core.h"

#include <internal.h>
#include <mp_internal.h>
#include <quaternion_constants.h>
#include <quaternion_data.h>
#include <encoded_sizes.h>
#include <prng.h>
#include <rng.h>

#define INBUF_SIZE 4096u
#define OUTBUF_SIZE 6144u
#define CONJ_CHUNK 48u

#define CONJ_NUM_LIMBS(bits) ((int)(((bits) + (NUM_BITS_LIMB - 1)) / (NUM_BITS_LIMB)))

static char inbuf[INBUF_SIZE];
static char outbuf[OUTBUF_SIZE];
static size_t inlen;
static size_t outlen;

/* Alg 3.18 inputs, parsed from inbuf */
static quat_ideal_t in_sk_ideal;
static quat_ideal_t in_ideal_chall;
static quat_alg_elem_t in_chall_split;
static quat_ideal_t in_ideal_commit;

/* Alg 3.18 outputs */
static quat_ideal_t out_sk_chall_ideal;
static quat_alg_elem_t out_resp_quat;
static quat_alg_elem_t out_sk_chall_quat;     /* after line 2 */
static quat_alg_elem_t out_sk_chall_quat_pre; /* before line 2: the SCA label */
static ibz_t out_norm;
static uint8_t last_valid;

/* Algorithm 3.23 line 12 target: a fixed secret matrix, a known challenge
 * coefficient supplied per call. */
static ibz_mat_2x2_t chall_msk;
static ibz_vec_2_t chall_vec;
static ibz_vec_2_t chall_out;
static uint8_t chall_msk_valid;

/* line 2 fast path: an alpha injected directly, so a campaign does not pay for a
 * full Algorithm 3.18 run per trace. */
static quat_alg_elem_t alpha_work;
static quat_alg_elem_t alpha_pre;
static uint8_t alpha_valid;

/* Algorithm 3.18 line 5 fast path: the product lattice and the accepted sample
 * vector are injected, so the board runs only the formation of alpha_rsp and
 * not the seconds-long rejection sampling that precedes it. */
static quat_lattice_t resp_lat;
static ibz_vec_4_t resp_x;
static quat_alg_elem_t resp_out;
static int resp_radius_bound;
static uint8_t resp_valid;
static uint8_t last_bounds_applied;

/* Keep the core reachable even if the linker gets clever about a static call. */
static void (*volatile keep_conj_core)(quat_ideal_t *,
                                      quat_alg_elem_t *,
                                      ibz_t *,
                                      quat_alg_elem_t *,
                                      quat_alg_elem_t *,
                                      const quat_ideal_t *,
                                      const quat_ideal_t *,
                                      const quat_alg_elem_t *,
                                      const quat_ideal_t *,
                                      prng_domain_ctx_t *,
                                      int) = conj_response_element;

/* ------------------------------------------------------------ text helpers */

static int
next_token(const char **cursor, char *token, size_t token_size)
{
    const char *s = *cursor;
    while (*s == '\n' || *s == '\r' || *s == ' ' || *s == '\t') {
        s++;
    }
    if (*s == '\0') {
        *cursor = s;
        return 0;
    }
    const char *start = s;
    while (*s != '\0' && *s != '\n' && *s != '\r' && *s != ' ' && *s != '\t') {
        s++;
    }
    size_t len = (size_t)(s - start);
    if (len + 1 > token_size) {
        return 0;
    }
    memcpy(token, start, len);
    token[len] = '\0';
    *cursor = s;
    return 1;
}

static int
parse_ibz(ibz_t *x, const char **cursor)
{
    char token[560]; /* IBZ_NLIMBS=60 limbs -> at most 480 hex digits + sign */
    if (!next_token(cursor, token, sizeof(token))) {
        return 0;
    }
    return ibz_set_from_str(x, token, 16);
}

static int
parse_ideal(quat_ideal_t *ideal, const char **cursor)
{
    return parse_ibz(&ideal->x, cursor) && parse_ibz(&ideal->y, cursor) && parse_ibz(&ideal->norm, cursor);
}

static int
line_starts_with(const char *s, const char *kw)
{
    while (*kw) {
        if (*s++ != *kw++) {
            return 0;
        }
    }
    return 1;
}

static int
parse_int(const char **cursor, int *out)
{
    const char *s = *cursor;
    while (*s == ' ' || *s == '\t') {
        s++;
    }
    int neg = 0;
    if (*s == '-') {
        neg = 1;
        s++;
    }
    if (*s < '0' || *s > '9') {
        return 0;
    }
    int v = 0;
    while (*s >= '0' && *s <= '9') {
        v = v * 10 + (*s - '0');
        s++;
    }
    *cursor = s;
    *out = neg ? -v : v;
    return 1;
}

/* Walk the vector-file header, collecting the declared bit bounds and leaving the
 * cursor on the first hex token.
 *
 * The bounds matter: ibz_t carries a declared bit bound that every downstream
 * ibz_t inherits, and ibz_set_from_str() re-derives it from the parsed magnitude
 * (mp.c:311), which is tighter than what the host had in memory. That difference
 * propagates to alpha_sk,chl and changes how many limbs ibz_neg walks at line 2,
 * so without restoring the bounds the board would execute a narrower Conjugate
 * than a real signer.
 *
 * Returns 1 if a bounds line was found, 0 otherwise; the cursor is advanced past
 * the header either way. */
static int
parse_header(const char **cursor, int *bounds, int nbounds)
{
    const char *s = *cursor;
    int have = 0;
    for (;;) {
        while (*s == '\n' || *s == '\r') {
            s++;
        }
        if (*s == '\0') {
            break;
        }
        if (*s == '#' || line_starts_with(s, "seed")) {
            while (*s != '\0' && *s != '\n') {
                s++;
            }
            continue;
        }
        if (line_starts_with(s, "bounds")) {
            const char *q = s + 6;
            have = 1;
            for (int i = 0; i < nbounds; i++) {
                if (!parse_int(&q, &bounds[i])) {
                    have = 0;
                    break;
                }
            }
            while (*s != '\0' && *s != '\n') {
                s++;
            }
            continue;
        }
        if (line_starts_with(s, "vector")) {
            while (*s != '\0' && *s != '\n') {
                s++;
            }
            break;
        }
        /* not a keyword: the hex payload starts here */
        break;
    }
    *cursor = s;
    return have;
}

static int
parse_elem(quat_alg_elem_t *e, const char **cursor)
{
    if (!parse_ibz(&e->denom, cursor)) {
        return 0;
    }
    for (int i = 0; i < 4; i++) {
        if (!parse_ibz(&e->coord.v[i], cursor)) {
            return 0;
        }
    }
    return 1;
}

static void
append_str(const char *s)
{
    while (*s && outlen + 1 < OUTBUF_SIZE) {
        outbuf[outlen++] = *s++;
    }
    outbuf[outlen] = '\0';
}

static void
append_ibz(const ibz_t *x)
{
    char tmp[560];
    if (!ibz_convert_to_str(x, tmp, 16)) {
        tmp[0] = '?';
        tmp[1] = '\0';
    }
    append_str(tmp);
    append_str("\n");
}

static void
append_u32_dec(uint32_t x)
{
    char tmp[16];
    unsigned pos = sizeof(tmp);
    tmp[--pos] = '\0';
    if (x == 0) {
        tmp[--pos] = '0';
    } else {
        while (x && pos > 0) {
            tmp[--pos] = (char)('0' + (x % 10));
            x /= 10;
        }
    }
    append_str(&tmp[pos]);
}

static void
append_kv(const char *name, uint32_t value)
{
    append_str(name);
    append_str("=");
    append_u32_dec(value);
    append_str("\n");
}

/* Raw limbs, in the order ibz_neg walks them. These are the words whose Hamming
 * weight the trace is expected to leak, so the host wants them verbatim rather
 * than reconstructed from the decimal/hex value. */
static void
append_limbs(const ibz_t *x)
{
    int n = CONJ_NUM_LIMBS(x->bitlen);
    append_u32_dec((uint32_t)x->bitlen);
    append_str(" ");
    append_u32_dec((uint32_t)n);
    for (int i = 0; i < n; i++) {
        append_str(" ");
        char tmp[12];
        uint32_t v = (uint32_t)x->limbs[i];
        for (int k = 7; k >= 0; k--) {
            unsigned nib = (v >> (4 * k)) & 0xfu;
            tmp[7 - k] = (char)(nib < 10 ? ('0' + nib) : ('a' + nib - 10));
        }
        tmp[8] = '\0';
        append_str(tmp);
    }
    append_str("\n");
}

/* ----------------------------------------------------------------- results */

static void
serialize_result(void)
{
    outlen = 0;
    outbuf[0] = '\0';
    if (!last_valid) {
        append_str("0\n");
        return;
    }
    append_str("1\n");

    /* alpha_sk,chl entering line 2 -- the side-channel ground truth */
    append_ibz(&out_sk_chall_quat_pre.denom);
    for (int i = 0; i < 4; i++) {
        append_ibz(&out_sk_chall_quat_pre.coord.v[i]);
    }
    /* the same coordinates as raw limbs, coord[0..3] */
    for (int i = 0; i < 4; i++) {
        append_limbs(&out_sk_chall_quat_pre.coord.v[i]);
    }

    /* alpha_sk,chl leaving line 2 */
    append_ibz(&out_sk_chall_quat.denom);
    for (int i = 0; i < 4; i++) {
        append_ibz(&out_sk_chall_quat.coord.v[i]);
    }

    /* Algorithm 3.18 result */
    append_ibz(&out_norm);
    append_ibz(&out_sk_chall_ideal.norm);
    append_ibz(&out_resp_quat.denom);
    for (int i = 0; i < 4; i++) {
        append_ibz(&out_resp_quat.coord.v[i]);
    }
}

/* ---------------------------------------------------------------- commands */

#if SS_VER == SS_VER_2_1
#define SS_ARGS uint8_t cmd, uint8_t scmd, uint8_t len, uint8_t *buf
#define SS_UNUSED_HDR                                                                                                  \
    (void)cmd;                                                                                                         \
    (void)scmd
#else
#define SS_ARGS uint8_t *buf, uint8_t len
#define SS_UNUSED_HDR ((void)0)
#endif

uint8_t
conj_info(SS_ARGS)
{
    (void)buf;
    (void)len;
    SS_UNUSED_HDR;
    uint8_t out[10];
    uint16_t f[5] = {
        (uint16_t)sizeof(quat_alg_elem_t),
        (uint16_t)sizeof(quat_ideal_t),
        (uint16_t)IBZ_NLIMBS,
        (uint16_t)CONJ_TRIG_MODE,
        (uint16_t)((CONJ_OUT_OF_PLACE ? 0x100u : 0u) | (uint16_t)CONJ_PROBE_REPS),
    };
    for (unsigned i = 0; i < 5; i++) {
        out[2 * i] = (uint8_t)f[i];
        out[2 * i + 1] = (uint8_t)(f[i] >> 8);
    }
    simpleserial_put('r', sizeof(out), out);
    return 0;
}

uint8_t
conj_reset(SS_ARGS)
{
    (void)buf;
    (void)len;
    SS_UNUSED_HDR;
    inlen = 0;
    outlen = 0;
    last_valid = 0;
    inbuf[0] = '\0';
    outbuf[0] = '\0';
    uint8_t ok = 1;
    simpleserial_put('r', 1, &ok);
    return 0;
}

uint8_t
conj_write(SS_ARGS)
{
    SS_UNUSED_HDR;
    uint8_t ok = 0;
    if (inlen + len + 1 < INBUF_SIZE) {
        memcpy(inbuf + inlen, buf, len);
        inlen += len;
        inbuf[inlen] = '\0';
        ok = 1;
    }
    simpleserial_put('r', 1, &ok);
    return 0;
}

uint8_t
conj_compute(SS_ARGS)
{
    (void)buf;
    (void)len;
    SS_UNUSED_HDR;
    uint8_t ok = 0;
    const char *cursor = inbuf;
    int bounds[14];
    int have_bounds;

    conj_profile_reset();
    last_valid = 0;

    have_bounds = parse_header(&cursor, bounds, 14);

    if (parse_ideal(&in_sk_ideal, &cursor) && parse_ideal(&in_ideal_chall, &cursor) &&
        parse_elem(&in_chall_split, &cursor) && parse_ideal(&in_ideal_commit, &cursor)) {

        last_bounds_applied = (uint8_t)have_bounds;
        if (have_bounds) {
            ibz_t *f[14] = {
                &in_sk_ideal.x,          &in_sk_ideal.y,          &in_sk_ideal.norm,
                &in_ideal_chall.x,       &in_ideal_chall.y,       &in_ideal_chall.norm,
                &in_chall_split.denom,   &in_chall_split.coord.v[0],
                &in_chall_split.coord.v[1], &in_chall_split.coord.v[2],
                &in_chall_split.coord.v[3],
                &in_ideal_commit.x,      &in_ideal_commit.y,      &in_ideal_commit.norm,
            };
            for (int i = 0; i < 14; i++) {
                if (bounds[i] > f[i]->bitlen && bounds[i] <= IBZ_NLIMBS * NUM_BITS_LIMB) {
                    ibz_set_bound(f[i], bounds[i]);
                }
            }
        }

        prng_domain_ctx_t prng_res;
        if (prng_domain_seed(&prng_res, DS_RESPONSE_DOMAIN) == 0) {
            conj_response_element(&out_sk_chall_ideal,
                                  &out_resp_quat,
                                  &out_norm,
                                  &out_sk_chall_quat,
                                  &out_sk_chall_quat_pre,
                                  &in_sk_ideal,
                                  &in_ideal_chall,
                                  &in_chall_split,
                                  &in_ideal_commit,
                                  &prng_res,
                                  RESPONSE_BITS - 1);
            prng_domain_clear(&prng_res);
            last_valid = 1;
            ok = 1;
        }
    }

    serialize_result();
    simpleserial_put('r', 1, &ok);
    return 0;
}

uint8_t
conj_seed_rng(SS_ARGS)
{
    SS_UNUSED_HDR;
    uint8_t ok = 0;
    if (len == 16) {
        uint8_t seed[48];
        for (int i = 0; i < 48; i++) {
            seed[i] = buf[i & 15] ^ (uint8_t)(0x9dU + 17U * (uint8_t)i);
        }
        randombytes_init(seed, 0, 256);
        /* re-derive the PRNG root and the default domain (Alg 3.18 line 1 uses it) */
        ok = (prng_seed() == 0);
    }
    simpleserial_put('r', 1, &ok);
    return 0;
}

uint8_t
conj_len(SS_ARGS)
{
    (void)buf;
    (void)len;
    SS_UNUSED_HDR;
    uint8_t out[2] = { (uint8_t)(outlen >> 8), (uint8_t)outlen };
    simpleserial_put('r', sizeof(out), out);
    return 0;
}

uint8_t
conj_read(SS_ARGS)
{
    SS_UNUSED_HDR;
    if (len < 2) {
        return 0x01;
    }
    uint16_t idx = ((uint16_t)buf[0] << 8) | buf[1];
    size_t off = (size_t)idx * CONJ_CHUNK;
    uint8_t chunk[CONJ_CHUNK];
    memset(chunk, 0, sizeof(chunk));
    if (off < outlen) {
        size_t remain = outlen - off;
        uint8_t n = (uint8_t)(remain < CONJ_CHUNK ? remain : CONJ_CHUNK);
        memcpy(chunk, outbuf + off, n);
    }
    simpleserial_put('r', CONJ_CHUNK, chunk);
    return 0;
}

/* 'k': load the secret matrix Msk from the ASCII buffer (4 hex values, row
 * major) and give every entry the CHALLENGE_BITS+1 bound keygen leaves them
 * with. Sent once; it stays fixed for the whole campaign, which is the point --
 * unlike alpha_sk,chl it is the same in every signature, so traces average. */
uint8_t
chall_set_secret(SS_ARGS)
{
    (void)buf;
    (void)len;
    SS_UNUSED_HDR;
    const char *cursor = inbuf;
    uint8_t ok = 1;
    for (int i = 0; i < 2 && ok; i++) {
        for (int j = 0; j < 2 && ok; j++) {
            /* ibz_set_from_str returns -1 on success, 0 on failure (mp.c:365),
             * so normalise it rather than casting it straight to the ack byte. */
            ok = parse_ibz(&chall_msk.m[i][j], &cursor) ? 1 : 0;
            if (ok && chall_msk.m[i][j].bitlen < CHALLENGE_BITS + 1) {
                ibz_set_bound(&chall_msk.m[i][j], CHALLENGE_BITS + 1);
            }
        }
    }
    chall_msk_valid = ok;
    simpleserial_put('r', 1, &ok);
    return 0;
}

/* 'e': run Algorithm 3.23 line 12 on a known challenge coefficient.
 *
 * The coefficient arrives as CHALLENGE_BYTES little-endian bytes in the command
 * payload rather than through the ASCII buffer, so a trace costs one short
 * serial exchange instead of a multi-kilobyte transfer. */
uint8_t
chall_eval(SS_ARGS)
{
    SS_UNUSED_HDR;
    uint8_t ok = 0;
    if (len == CHALLENGE_BYTES && chall_msk_valid) {
        digit_t w[CHALLENGE_BYTES / sizeof(digit_t)];
        for (unsigned i = 0; i < CHALLENGE_BYTES / sizeof(digit_t); i++) {
            w[i] = (digit_t)buf[4 * i] | ((digit_t)buf[4 * i + 1] << 8) | ((digit_t)buf[4 * i + 2] << 16) |
                   ((digit_t)buf[4 * i + 3] << 24);
        }
        chall_profile_reset();
        ibz_set(&chall_vec.v[0], 1, CHALLENGE_BITS + 1);
        ibz_copy_digits(&chall_vec.v[1], w, CHALLENGE_BYTES / sizeof(digit_t));
        chall_mat_2x2_eval_triggered(&chall_out, &chall_msk, &chall_vec);
        ok = 1;
    }
    simpleserial_put('r', 1, &ok);
    return 0;
}

/* 'f': the result of the last 'e', plus the profile. Used to verify the board
 * against the host and to read the cycle counts; it reveals the secret, so it
 * exists for profiling runs only. */
uint8_t
chall_read(SS_ARGS)
{
    (void)buf;
    (void)len;
    SS_UNUSED_HDR;
    outlen = 0;
    outbuf[0] = '\0';
    append_ibz(&chall_out.v[0]);
    append_ibz(&chall_out.v[1]);
    append_kv("eval", chall_profile.eval_cycles);
    for (int i = 0; i < 4; i++) {
        char name[12] = "mul0";
        name[3] = (char)('0' + i);
        append_kv(name, chall_profile.mul_cycles[i]);
    }
    append_kv("secret_limbs", chall_profile.secret_limbs);
    append_kv("known_limbs", chall_profile.known_limbs);
    append_kv("secret_bitlen", chall_profile.secret_bitlen);
    append_kv("known_bitlen", chall_profile.known_bitlen);
    uint8_t ok = 1;
    simpleserial_put('r', 1, &ok);
    return 0;
}

/* 'u': load one alpha_sk,chl into the line 2 working register.
 *
 * Not 'y': simpleserial_init reserves 'v', 'w' and 'y' for check_version,
 * ss_get_commands and ss_num_commands. A second registration of a reserved
 * letter is shadowed by the built-in, which answers with the command count --
 * 0x11 here -- and silently desynchronises the link.
 *
 * Buffer layout matches conj_vector_host "alphas": a bounds line of five bitlen
 * values, then denom and coord0..3 as hex. The bounds are restored with
 * ibz_set_bound for the same reason the Algorithm 3.18 vectors carry them --
 * ibz_set_from_str re-derives a tighter bound from the parsed magnitude, which
 * would change how many limbs ibz_neg walks. */
uint8_t
conj_set_alpha(SS_ARGS)
{
    (void)buf;
    (void)len;
    SS_UNUSED_HDR;
    const char *cursor = inbuf;
    int bounds[5];
    int have = 0;

    /* Only a comment, a blank line, or an "a <bounds>" line is header. Matching a
     * bare letter would be wrong: a hex coordinate can legitimately begin with
     * 'a' or 'c', so the bounds line is recognised by "a " specifically. */
    while (*cursor == '#' || *cursor == '\n' || *cursor == '\r' ||
           (*cursor == 'a' && (cursor[1] == ' ' || cursor[1] == '\t'))) {
        if (*cursor == 'a') {
            const char *q = cursor + 1;
            have = 1;
            for (int i = 0; i < 5; i++) {
                if (!parse_int(&q, &bounds[i])) {
                    have = 0;
                    break;
                }
            }
        }
        while (*cursor != '\0' && *cursor != '\n') {
            cursor++;
        }
        if (*cursor == '\n') {
            cursor++;
        }
    }

    uint8_t ok = parse_elem(&alpha_work, &cursor) ? 1 : 0;
    if (ok && have) {
        ibz_t *f[5] = { &alpha_work.denom,      &alpha_work.coord.v[0], &alpha_work.coord.v[1],
                        &alpha_work.coord.v[2], &alpha_work.coord.v[3] };
        for (int i = 0; i < 5; i++) {
            if (bounds[i] > f[i]->bitlen && bounds[i] <= IBZ_NLIMBS * NUM_BITS_LIMB) {
                ibz_set_bound(f[i], bounds[i]);
            }
        }
    }
    alpha_valid = ok;
    simpleserial_put('r', 1, &ok);
    return 0;
}

/* 'z': run Algorithm 3.18 line 2 on the loaded alpha. The trigger fires per
 * CONJ_TRIG_MODE exactly as it does inside the full algorithm. */
uint8_t
conj_run_line2(SS_ARGS)
{
    (void)buf;
    (void)len;
    SS_UNUSED_HDR;
    uint8_t ok = 0;
    if (alpha_valid) {
        conj_profile_reset();
        quat_alg_elem_copy(&alpha_pre, &alpha_work);
        conj_run_conjugate(&alpha_work, NULL);
        ok = 1;
    }
    simpleserial_put('r', 1, &ok);
    return 0;
}

/* 'q': the limbs of the value that entered line 2, for scoring the attack.
 * Not 'w' -- simpleserial_init reserves that for ss_get_commands, and a second
 * registration is shadowed by the built-in. */
uint8_t
conj_read_alpha(SS_ARGS)
{
    (void)buf;
    (void)len;
    SS_UNUSED_HDR;
    outlen = 0;
    outbuf[0] = '\0';
    for (int i = 0; i < 4; i++) {
        append_limbs(&alpha_pre.coord.v[i]);
    }
    for (int i = 0; i < 4; i++) {
        append_ibz(&alpha_work.coord.v[i]);
    }
    append_kv("conjugate", conj_profile.conjugate_cycles);
    uint8_t ok = 1;
    simpleserial_put('r', 1, &ok);
    return 0;
}

/* 'b': load the product lattice, the accepted sample vector and the radius bound
 * for one line 5 formation. Buffer layout, hex unless noted:
 *   r <radius_bound> <denom bitlen> <16 basis bitlens> <4 x bitlens>   (decimal)
 *   denom, basis[0][0..3], basis[1][0..3], basis[2][0..3], basis[3][0..3],
 *   x[0..3] */
uint8_t
resp_set_input(SS_ARGS)
{
    (void)buf;
    (void)len;
    SS_UNUSED_HDR;
    const char *cursor = inbuf;
    int bounds[21];
    int have = 0;

    while (*cursor == '#' || *cursor == '\n' || *cursor == '\r' ||
           (*cursor == 'r' && (cursor[1] == ' ' || cursor[1] == '\t'))) {
        if (*cursor == 'r') {
            const char *q = cursor + 1;
            have = parse_int(&q, &resp_radius_bound);
            for (int i = 0; i < 21 && have; i++) {
                have = parse_int(&q, &bounds[i]);
            }
        }
        while (*cursor != '\0' && *cursor != '\n') {
            cursor++;
        }
        if (*cursor == '\n') {
            cursor++;
        }
    }

    uint8_t ok = parse_ibz(&resp_lat.denom, &cursor) ? 1 : 0;
    for (int i = 0; i < 4 && ok; i++) {
        for (int j = 0; j < 4 && ok; j++) {
            ok = parse_ibz(&resp_lat.basis.m[i][j], &cursor) ? 1 : 0;
        }
    }
    for (int i = 0; i < 4 && ok; i++) {
        ok = parse_ibz(&resp_x.v[i], &cursor) ? 1 : 0;
    }

    if (ok && have) {
        ibz_t *f[21];
        int k = 0;
        f[k++] = &resp_lat.denom;
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                f[k++] = &resp_lat.basis.m[i][j];
            }
        }
        for (int i = 0; i < 4; i++) {
            f[k++] = &resp_x.v[i];
        }
        for (int i = 0; i < 21; i++) {
            if (bounds[i] > f[i]->bitlen && bounds[i] <= IBZ_NLIMBS * NUM_BITS_LIMB) {
                ibz_set_bound(f[i], bounds[i]);
            }
        }
    }
    resp_valid = ok;
    simpleserial_put('r', 1, &ok);
    return 0;
}

/* 'm': form alpha_rsp from the loaded lattice and sample vector. The trigger
 * fires per CONJ_TRIG_MODE exactly as it does inside the full algorithm. */
uint8_t
resp_run(SS_ARGS)
{
    (void)buf;
    (void)len;
    SS_UNUSED_HDR;
    uint8_t ok = 0;
    if (resp_valid) {
        conj_profile_reset();
        conj_run_response_formation(&resp_out, &resp_lat, &resp_x, resp_radius_bound);
        ok = 1;
    }
    simpleserial_put('r', 1, &ok);
    return 0;
}

/* 'o': alpha_rsp as limbs and as values, for scoring. */
uint8_t
resp_read(SS_ARGS)
{
    (void)buf;
    (void)len;
    SS_UNUSED_HDR;
    outlen = 0;
    outbuf[0] = '\0';
    for (int i = 0; i < 4; i++) {
        append_limbs(&resp_out.coord.v[i]);
    }
    append_ibz(&resp_out.denom);
    for (int i = 0; i < 4; i++) {
        append_ibz(&resp_out.coord.v[i]);
    }
    append_kv("formation", conj_profile.formation_cycles);
    append_kv("eval", conj_profile.eval_cycles);
    append_kv("normalize", conj_profile.normalize_cycles);
    uint8_t ok = 1;
    simpleserial_put('r', 1, &ok);
    return 0;
}

uint8_t
conj_profile_read(SS_ARGS)
{
    (void)buf;
    (void)len;
    SS_UNUSED_HDR;
    outlen = 0;
    outbuf[0] = '\0';
    append_kv("trig_mode", (uint32_t)CONJ_TRIG_MODE);
    append_kv("out_of_place", (uint32_t)CONJ_OUT_OF_PLACE);
    append_kv("probe_reps", (uint32_t)CONJ_PROBE_REPS);
    append_kv("total", conj_profile.total_cycles);
    append_kv("intersect", conj_profile.intersect_cycles);
    append_kv("equiv_coprime", conj_profile.equiv_coprime_cycles);
    append_kv("conjugate", conj_profile.conjugate_cycles);
    append_kv("ideal_mul", conj_profile.ideal_mul_cycles);
    append_kv("sample", conj_profile.sample_cycles);
    for (int i = 0; i < 4; i++) {
        char name[20] = "conj_coord0";
        name[10] = (char)('0' + i);
        append_kv(name, conj_profile.conj_coord_cycles[i]);
    }
    for (int i = 0; i < 4; i++) {
        char name[20] = "conj_limbs0";
        name[10] = (char)('0' + i);
        append_kv(name, conj_profile.conj_limbs[i]);
    }
    for (int i = 0; i < 4; i++) {
        char name[20] = "conj_bitlen0";
        name[11] = (char)('0' + i);
        append_kv(name, conj_profile.conj_bitlen[i]);
    }
    append_kv("equiv_ok", conj_profile.equiv_ok);
    append_kv("sample_ok", conj_profile.sample_ok);
    append_kv("bounds_applied", (uint32_t)last_bounds_applied);
    uint8_t ok = 1;
    simpleserial_put('r', 1, &ok);
    return 0;
}

int
main(void)
{
    platform_init();
    init_uart();
    trigger_setup();
    conj_cycles_init();

    quat_ideal_init(&in_sk_ideal);
    quat_ideal_init(&in_ideal_chall);
    quat_ideal_init(&in_ideal_commit);
    quat_alg_elem_init(&in_chall_split);
    quat_lattice_init(&resp_lat);
    ibz_vec_4_init(&resp_x);
    quat_alg_elem_init(&resp_out);
    resp_valid = 0;
    quat_alg_elem_init(&alpha_work);
    quat_alg_elem_init(&alpha_pre);
    alpha_valid = 0;
    ibz_mat_2x2_init(&chall_msk);
    ibz_vec_2_init(&chall_vec);
    ibz_vec_2_init(&chall_out);
    chall_msk_valid = 0;
    quat_ideal_init(&out_sk_chall_ideal);
    quat_alg_elem_init(&out_resp_quat);
    quat_alg_elem_init(&out_sk_chall_quat);
    quat_alg_elem_init(&out_sk_chall_quat_pre);
    ibz_init(&out_norm);
    last_valid = 0;

    /* The board RNG must be seeded before Alg 3.18 line 1, which draws from
     * PRNG_default_domain. The host normally overrides this with 'n'. */
    randombytes_init(0, 0, 256);
    prng_seed();

    putch('S');
    putch('Q');
    putch('I');
    putch('c');
    putch('o');
    putch('n');
    putch('j');
    putch('\n');

    simpleserial_init();
    simpleserial_addcmd('i', 0, conj_info);
    simpleserial_addcmd('x', 0, conj_reset);
    simpleserial_addcmd('a', CONJ_CHUNK, conj_write);
    simpleserial_addcmd('c', 0, conj_compute);
    simpleserial_addcmd('n', 16, conj_seed_rng);
    simpleserial_addcmd('l', 0, conj_len);
    simpleserial_addcmd('g', 2, conj_read);
    simpleserial_addcmd('p', 0, conj_profile_read);
    simpleserial_addcmd('k', 0, chall_set_secret);
    simpleserial_addcmd('e', CHALLENGE_BYTES, chall_eval);
    simpleserial_addcmd('f', 0, chall_read);
    simpleserial_addcmd('u', 0, conj_set_alpha);
    simpleserial_addcmd('z', 0, conj_run_line2);
    simpleserial_addcmd('q', 0, conj_read_alpha);
    simpleserial_addcmd('b', 0, resp_set_input);
    simpleserial_addcmd('m', 0, resp_run);
    simpleserial_addcmd('o', 0, resp_read);

    (void)keep_conj_core;
    while (1) {
        simpleserial_get();
    }
    return 0;
}
