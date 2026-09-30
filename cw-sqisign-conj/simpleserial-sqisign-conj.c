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
 */
#include "hal.h"
#include "simpleserial.h"

#include <stddef.h>
#include <stdint.h>
#include <string.h>

#include "conj_core.h"

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

    conj_profile_reset();
    last_valid = 0;

    if (parse_ideal(&in_sk_ideal, &cursor) && parse_ideal(&in_ideal_chall, &cursor) &&
        parse_elem(&in_chall_split, &cursor) && parse_ideal(&in_ideal_commit, &cursor)) {

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

    (void)keep_conj_core;
    while (1) {
        simpleserial_get();
    }
    return 0;
}
