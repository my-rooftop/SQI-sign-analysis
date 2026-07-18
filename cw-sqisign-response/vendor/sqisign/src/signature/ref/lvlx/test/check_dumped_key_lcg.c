#include <api.h>
#include <rng.h>
#include <sig.h>

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MSG_LEN 32

static uint32_t rng_state = 0x12345678u;
static unsigned long progress_counts[256];

void
sqisign_cw_progress(unsigned char tag)
{
    progress_counts[tag]++;
}

void
randombytes_init(unsigned char *entropy_input, unsigned char *personalization_string, int security_strength)
{
    (void)personalization_string;
    (void)security_strength;
    rng_state = 0x12345678u;
    if (entropy_input) {
        for (int i = 0; i < 48; i++) {
            rng_state = (rng_state * 1664525u + 1013904223u) ^ entropy_input[i];
        }
    }
}

int
randombytes_select(unsigned char *x, unsigned long long xlen)
{
    return randombytes(x, xlen);
}

int
randombytes(unsigned char *x, unsigned long long xlen)
{
    for (unsigned long long i = 0; i < xlen; i++) {
        rng_state = rng_state * 1664525u + 1013904223u;
        x[i] = (unsigned char)(rng_state >> 24);
    }
    return 0;
}

static int
read_exact(const char *path, unsigned char *buf, size_t len)
{
    FILE *f = fopen(path, "rb");
    if (!f) {
        perror(path);
        return 0;
    }
    size_t got = fread(buf, 1, len, f);
    int extra = fgetc(f);
    fclose(f);
    if (got != len || extra != EOF) {
        fprintf(stderr, "bad key dump size: got %zu bytes, expected %zu\n", got + (extra == EOF ? 0 : 1), len);
        return 0;
    }
    return 1;
}

static void
seed_board_rng(void)
{
    unsigned char seed[48];
    for (unsigned i = 0; i < sizeof(seed); i++) {
        seed[i] = (unsigned char)i;
    }
    randombytes_init(seed, NULL, 256);
}

static void
print_count(unsigned char tag, const char *label)
{
    printf("%c %-24s %lu\n", tag, label, progress_counts[tag]);
}

int
main(int argc, char **argv)
{
    const char *path = argc > 1 ? argv[1] : "/tmp/sqisign_board_key.bin";
    unsigned char dumped_pk[CRYPTO_PUBLICKEYBYTES];
    unsigned char dumped_sk[CRYPTO_SECRETKEYBYTES];
    unsigned char msg[MSG_LEN] = { 0 };
    unsigned char opened[MSG_LEN];
    unsigned char sm[CRYPTO_BYTES + MSG_LEN];
    unsigned long long smlen = 0;
    unsigned long long opened_len = 0;

    unsigned char keybuf[CRYPTO_PUBLICKEYBYTES + CRYPTO_SECRETKEYBYTES];
    if (!read_exact(path, keybuf, sizeof(keybuf))) {
        return 1;
    }
    memcpy(dumped_pk, keybuf, CRYPTO_PUBLICKEYBYTES);
    memcpy(dumped_sk, keybuf + CRYPTO_PUBLICKEYBYTES, CRYPTO_SECRETKEYBYTES);

    memset(progress_counts, 0, sizeof(progress_counts));
    seed_board_rng();
    int sign_ret = sqisign_sign(sm, &smlen, msg, sizeof(msg), dumped_sk);

    printf("sign_ret %d smlen %llu\n", sign_ret, smlen);
    print_count('a', "sign attempts");
    print_count('P', "aux sample start");
    print_count('Q', "aux sample done");
    print_count('N', "represent loop heartbeat");
    print_count('M', "final split reject");
    print_count('Y', "theta chain success");
    print_count('m', "sign finished");

    if (sign_ret != 0) {
        return 4;
    }
    if (sqisign_open(opened, &opened_len, sm, smlen, dumped_pk) != 0) {
        fprintf(stderr, "open with dumped pk failed\n");
        return 5;
    }
    if (opened_len != sizeof(msg) || memcmp(opened, msg, sizeof(msg)) != 0) {
        fprintf(stderr, "opened message mismatch\n");
        return 6;
    }

    printf("open ok\n");
    return 0;
}
