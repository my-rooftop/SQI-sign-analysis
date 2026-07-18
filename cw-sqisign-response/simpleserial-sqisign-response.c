#include "hal.h"
#include "simpleserial.h"

#include <stddef.h>
#include <stdint.h>
#include <string.h>

#include "response_core.h"
#include "intbig_internal.h"
#include <quaternion_data.h>
#include <rng.h>

#define INBUF_SIZE 8192u
#define OUTBUF_SIZE 4096u
#define RESP_CHUNK 48u

#ifndef RESP_TRIGGER_MODE
#define RESP_TRIGGER_MODE 0
#endif

static char inbuf[INBUF_SIZE];
static char outbuf[OUTBUF_SIZE];
static size_t inlen;
static size_t outlen;
static quat_alg_elem_t last_resp;
static ibz_t last_lattice_content;
static uint8_t last_valid;

static void (*volatile keep_response_core)(quat_alg_elem_t *,
                                           ibz_t *,
                                           const quat_left_ideal_t *,
                                           const quat_left_ideal_t *,
                                           const quat_left_ideal_t *) =
    sqisign_response_quat_element_core;

#if SS_VER == SS_VER_2_1
uint8_t response_info(uint8_t cmd, uint8_t scmd, uint8_t len, uint8_t *buf)
#else
uint8_t response_info(uint8_t *buf, uint8_t len)
#endif
{
    (void)buf;
    (void)len;
#if SS_VER == SS_VER_2_1
    (void)cmd;
    (void)scmd;
#endif
    uint8_t out[8];
    uint16_t sizes[4] = {
        (uint16_t)sizeof(quat_alg_elem_t),
        (uint16_t)sizeof(quat_left_ideal_t),
        (uint16_t)sizeof(quat_lattice_t),
        (uint16_t)((uintptr_t)keep_response_core & 0xffffu),
    };
    for (unsigned i = 0; i < 4; i++) {
        out[2 * i] = (uint8_t)sizes[i];
        out[2 * i + 1] = (uint8_t)(sizes[i] >> 8);
    }
    simpleserial_put('r', sizeof(out), out);
    return 0;
}

static int next_token(const char **cursor, char *token, size_t token_size)
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

static int parse_ibz(ibz_t *x, const char **cursor)
{
    char token[1024];
    if (!next_token(cursor, token, sizeof(token))) {
        return 0;
    }
    return ibz_set_from_str(x, token, 16);
}

static int parse_ideal(quat_left_ideal_t *ideal, const char **cursor)
{
    ideal->parent_order = &MAXORD_O0;
    if (!parse_ibz(&ideal->norm, cursor)) {
        return 0;
    }
    if (!parse_ibz(&ideal->lattice.denom, cursor)) {
        return 0;
    }
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            if (!parse_ibz(&ideal->lattice.basis[i][j], cursor)) {
                return 0;
            }
        }
    }
    return 1;
}

static void append_str(const char *s)
{
    while (*s && outlen + 1 < OUTBUF_SIZE) {
        outbuf[outlen++] = *s++;
    }
    outbuf[outlen] = '\0';
}

static void append_ibz(const ibz_t *x)
{
    char tmp[1024];
    ibz_convert_to_str(x, tmp, 16);
    append_str(tmp);
    append_str("\n");
}

static void append_u32_dec(uint32_t x)
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

static void append_profile_line(const char *name, uint32_t value)
{
    append_str(name);
    append_str("=");
    append_u32_dec(value);
    append_str("\n");
}

static void serialize_response(const ibz_t *lattice_content,
                               const quat_alg_elem_t *resp,
                               uint8_t do_trigger)
{
    outlen = 0;
    outbuf[0] = '\0';
    if (do_trigger) {
        trigger_high();
    }
    append_str("1\n");
    append_ibz(lattice_content);
    append_ibz(&resp->denom);
    for (int i = 0; i < 4; i++) {
        append_ibz(&resp->coord[i]);
    }
    if (do_trigger) {
        trigger_low();
    }
}

#if SS_VER == SS_VER_2_1
uint8_t response_reset(uint8_t cmd, uint8_t scmd, uint8_t len, uint8_t *buf)
#else
uint8_t response_reset(uint8_t *buf, uint8_t len)
#endif
{
    (void)buf;
    (void)len;
#if SS_VER == SS_VER_2_1
    (void)cmd;
    (void)scmd;
#endif
    inlen = 0;
    outlen = 0;
    last_valid = 0;
    inbuf[0] = '\0';
    outbuf[0] = '\0';
    uint8_t ok = 1;
    simpleserial_put('r', 1, &ok);
    return 0;
}

#if SS_VER == SS_VER_2_1
uint8_t response_write(uint8_t cmd, uint8_t scmd, uint8_t len, uint8_t *buf)
#else
uint8_t response_write(uint8_t *buf, uint8_t len)
#endif
{
#if SS_VER == SS_VER_2_1
    (void)cmd;
    (void)scmd;
#endif
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

#if SS_VER == SS_VER_2_1
uint8_t response_compute(uint8_t cmd, uint8_t scmd, uint8_t len, uint8_t *buf)
#else
uint8_t response_compute(uint8_t *buf, uint8_t len)
#endif
{
    (void)buf;
    (void)len;
#if SS_VER == SS_VER_2_1
    (void)cmd;
    (void)scmd;
#endif
    uint32_t t_handler = response_cycles_read();
    uint32_t t = t_handler;
    uint8_t ok = 0;
    response_profile_reset();
    const char *cursor = inbuf;
    quat_left_ideal_t secret_ideal;
    quat_left_ideal_t challenge_ideal;
    quat_left_ideal_t commit_ideal;
    quat_alg_elem_t resp;
    ibz_t lattice_content;

    quat_left_ideal_init(&secret_ideal);
    quat_left_ideal_init(&challenge_ideal);
    quat_left_ideal_init(&commit_ideal);
    quat_alg_elem_init(&resp);
    ibz_init(&lattice_content);

    outlen = 0;
    outbuf[0] = '\0';

    if (parse_ideal(&secret_ideal, &cursor) &&
        parse_ideal(&challenge_ideal, &cursor) &&
        parse_ideal(&commit_ideal, &cursor)) {
        response_profile.parse_cycles = response_cycles_read() - t;
#if RESP_TRIGGER_MODE == 0
        trigger_high();
#endif
        t = response_cycles_read();
        sqisign_response_quat_element_core(&resp,
                                           &lattice_content,
                                           &secret_ideal,
                                           &challenge_ideal,
                                           &commit_ideal);
        response_profile.core_cycles = response_cycles_read() - t;
#if RESP_TRIGGER_MODE == 0
        trigger_low();
#endif
        ok = 1;
        quat_alg_elem_copy(&last_resp, &resp);
        ibz_copy(&last_lattice_content, &lattice_content);
        last_valid = 1;
        t = response_cycles_read();
        serialize_response(&lattice_content, &resp, RESP_TRIGGER_MODE == 1);
        response_profile.serialize_cycles = response_cycles_read() - t;
    } else {
        append_str("0\n");
    }

    quat_alg_elem_finalize(&resp);
    ibz_finalize(&lattice_content);
    quat_left_ideal_finalize(&commit_ideal);
    quat_left_ideal_finalize(&challenge_ideal);
    quat_left_ideal_finalize(&secret_ideal);

    response_profile.handler_cycles = response_cycles_read() - t_handler;
    simpleserial_put('r', 1, &ok);
    return 0;
}

#if SS_VER == SS_VER_2_1
uint8_t response_serialize_last(uint8_t cmd, uint8_t scmd, uint8_t len, uint8_t *buf)
#else
uint8_t response_serialize_last(uint8_t *buf, uint8_t len)
#endif
{
    (void)buf;
    (void)len;
#if SS_VER == SS_VER_2_1
    (void)cmd;
    (void)scmd;
#endif
    uint8_t ok = 0;
    if (last_valid) {
        serialize_response(&last_lattice_content, &last_resp, 1);
        ok = 1;
    } else {
        outlen = 0;
        outbuf[0] = '\0';
        append_str("0\n");
    }
    simpleserial_put('r', 1, &ok);
    return 0;
}

#if SS_VER == SS_VER_2_1
uint8_t response_set_calibration(uint8_t cmd, uint8_t scmd, uint8_t len, uint8_t *buf)
#else
uint8_t response_set_calibration(uint8_t *buf, uint8_t len)
#endif
{
#if SS_VER == SS_VER_2_1
    (void)cmd;
    (void)scmd;
#endif
    uint8_t ok = 0;
    if (len == 4) {
        digit_t limb[1];
        limb[0] = ((uint32_t)buf[0] << 24) |
                  ((uint32_t)buf[1] << 16) |
                  ((uint32_t)buf[2] << 8) |
                  (uint32_t)buf[3];
        ibz_set(&last_lattice_content, 1);
        ibz_set(&last_resp.denom, 1);
        ibz_copy_digits(&last_resp.coord[0], limb, 1);
        ibz_set(&last_resp.coord[1], 0);
        ibz_set(&last_resp.coord[2], 0);
        ibz_set(&last_resp.coord[3], 0);
        last_valid = 1;
        ok = 1;
    }
    simpleserial_put('r', 1, &ok);
    return 0;
}

#if SS_VER == SS_VER_2_1
uint8_t response_seed_rng(uint8_t cmd, uint8_t scmd, uint8_t len, uint8_t *buf)
#else
uint8_t response_seed_rng(uint8_t *buf, uint8_t len)
#endif
{
#if SS_VER == SS_VER_2_1
    (void)cmd;
    (void)scmd;
#endif
    uint8_t seed[48];
    uint8_t ok = 0;
    if (len == 16) {
        for (int i = 0; i < 48; i++) {
            seed[i] = buf[i & 15] ^ (uint8_t)(0x9dU + 17U * (uint8_t)i);
        }
        randombytes_init(seed, 0, 256);
        ok = 1;
    }
    simpleserial_put('r', 1, &ok);
    return 0;
}

#if SS_VER == SS_VER_2_1
uint8_t response_read(uint8_t cmd, uint8_t scmd, uint8_t len, uint8_t *buf)
#else
uint8_t response_read(uint8_t *buf, uint8_t len)
#endif
{
#if SS_VER == SS_VER_2_1
    (void)cmd;
    (void)scmd;
#endif
    if (len < 2) {
        return 0x01;
    }
    uint16_t idx = ((uint16_t)buf[0] << 8) | buf[1];
    size_t off = (size_t)idx * RESP_CHUNK;
    uint8_t chunk[RESP_CHUNK];
    memset(chunk, 0, sizeof(chunk));
    if (off < outlen) {
        size_t remain = outlen - off;
        uint8_t n = (uint8_t)(remain < RESP_CHUNK ? remain : RESP_CHUNK);
        memcpy(chunk, outbuf + off, n);
    }
    simpleserial_put('r', RESP_CHUNK, chunk);
    return 0;
}

#if SS_VER == SS_VER_2_1
uint8_t response_len(uint8_t cmd, uint8_t scmd, uint8_t len, uint8_t *buf)
#else
uint8_t response_len(uint8_t *buf, uint8_t len)
#endif
{
    (void)buf;
    (void)len;
#if SS_VER == SS_VER_2_1
    (void)cmd;
    (void)scmd;
#endif
    uint8_t out[2] = {
        (uint8_t)(outlen >> 8),
        (uint8_t)outlen,
    };
    simpleserial_put('r', sizeof(out), out);
    return 0;
}

#if SS_VER == SS_VER_2_1
uint8_t response_profile_read(uint8_t cmd, uint8_t scmd, uint8_t len, uint8_t *buf)
#else
uint8_t response_profile_read(uint8_t *buf, uint8_t len)
#endif
{
    (void)buf;
    (void)len;
#if SS_VER == SS_VER_2_1
    (void)cmd;
    (void)scmd;
#endif
    outlen = 0;
    outbuf[0] = '\0';
    append_profile_line("parse", response_profile.parse_cycles);
    append_profile_line("core", response_profile.core_cycles);
    append_profile_line("serialize", response_profile.serialize_cycles);
    append_profile_line("handler", response_profile.handler_cycles);
    append_profile_line("inter", response_profile.inter_cycles);
    append_profile_line("conjugate", response_profile.conjugate_cycles);
    append_profile_line("lattice_intersect", response_profile.lattice_intersect_cycles);
    append_profile_line("content_mul", response_profile.content_mul_cycles);
    append_profile_line("sample_total", response_profile.sample_total_cycles);
    append_profile_line("sample_setup", response_profile.sample_setup_cycles);
    append_profile_line("sample_rejection", response_profile.sample_rejection_cycles);
    append_profile_line("sample_finalize", response_profile.sample_finalize_cycles);
    append_profile_line("normalize_total", response_profile.normalize_total_cycles);
    append_profile_line("normalize_gcd", response_profile.normalize_gcd_cycles);
    append_profile_line("normalize_coord0_div", response_profile.normalize_coord0_div_cycles);
    append_profile_line("normalize_coord_rest_div", response_profile.normalize_coord_rest_div_cycles);
    append_profile_line("normalize_sign", response_profile.normalize_sign_cycles);
    append_profile_line("rejection_iters", response_profile.rejection_iters);
    uint8_t ok = 1;
    simpleserial_put('r', 1, &ok);
    return 0;
}

int main(void)
{
    platform_init();
    init_uart();
    trigger_setup();
    response_cycles_init();
    quat_alg_elem_init(&last_resp);
    ibz_init(&last_lattice_content);
    last_valid = 0;

    putch('S');
    putch('Q');
    putch('I');
    putch('r');
    putch('e');
    putch('s');
    putch('p');
    putch('\n');

    simpleserial_init();
    simpleserial_addcmd('i', 0, response_info);
    simpleserial_addcmd('x', 0, response_reset);
    simpleserial_addcmd('a', RESP_CHUNK, response_write);
    simpleserial_addcmd('c', 0, response_compute);
    simpleserial_addcmd('q', 4, response_set_calibration);
    simpleserial_addcmd('n', 16, response_seed_rng);
    simpleserial_addcmd('s', 0, response_serialize_last);
    simpleserial_addcmd('g', 2, response_read);
    simpleserial_addcmd('l', 0, response_len);
    simpleserial_addcmd('p', 0, response_profile_read);

    while (1) {
        simpleserial_get();
    }
}
