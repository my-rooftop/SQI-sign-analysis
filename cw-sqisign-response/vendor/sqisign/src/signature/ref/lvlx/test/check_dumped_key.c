#include <api.h>
#include <rng.h>
#include <sig.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MSG_LEN 32

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
seed_rng(void)
{
    unsigned char seed[48];
    for (unsigned i = 0; i < sizeof(seed); i++) {
        seed[i] = (unsigned char)i;
    }
    randombytes_init(seed, NULL, 256);
}

int
main(int argc, char **argv)
{
    const char *path = argc > 1 ? argv[1] : "/tmp/sqisign_board_key.bin";
    unsigned char dumped_pk[CRYPTO_PUBLICKEYBYTES];
    unsigned char dumped_sk[CRYPTO_SECRETKEYBYTES];
    unsigned char pc_pk[CRYPTO_PUBLICKEYBYTES];
    unsigned char pc_sk[CRYPTO_SECRETKEYBYTES];
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

    seed_rng();
    if (sqisign_keypair(pc_pk, pc_sk) != 0) {
        fprintf(stderr, "PC keygen failed\n");
        return 1;
    }
    int matches_pc = 1;
    if (memcmp(dumped_pk, pc_pk, CRYPTO_PUBLICKEYBYTES) != 0) {
        fprintf(stderr, "dumped pk differs from PC keygen with same seed\n");
        matches_pc = 0;
    }
    if (memcmp(dumped_sk, pc_sk, CRYPTO_SECRETKEYBYTES) != 0) {
        fprintf(stderr, "dumped sk differs from PC keygen with same seed\n");
        matches_pc = 0;
    }

    seed_rng();
    if (sqisign_sign(sm, &smlen, msg, sizeof(msg), dumped_sk) != 0) {
        fprintf(stderr, "sign with dumped sk failed\n");
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

    if (matches_pc) {
        printf("dumped key matches PC keygen and signs/verifies\n");
    } else {
        printf("dumped key differs from PC deterministic keygen, but signs/verifies\n");
    }
    return 0;
}
