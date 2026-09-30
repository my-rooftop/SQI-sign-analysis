#include "rng.h"
#include <stdint.h>

static uint32_t rng_state = 0x13579bdfu;

void randombytes_init(unsigned char *entropy_input,
                      unsigned char *personalization_string,
                      int security_strength)
{
    (void)personalization_string;
    (void)security_strength;
    rng_state = 0x13579bdfu;
    if (entropy_input) {
        for (int i = 0; i < 48; i++) {
            rng_state = (rng_state * 1664525u + 1013904223u) ^ entropy_input[i];
        }
    }
}

int randombytes_select(unsigned char *x, unsigned long long xlen)
{
    return randombytes(x, xlen);
}

int randombytes(unsigned char *x, unsigned long long xlen)
{
    for (unsigned long long i = 0; i < xlen; i++) {
        rng_state = rng_state * 1664525u + 1013904223u;
        x[i] = (unsigned char)(rng_state >> 24);
    }
    return 0;
}
