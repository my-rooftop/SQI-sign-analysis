/* Not shipped by scripts/gen_pqm4_sources.sh: pqm4 supplies randombytes() on the
 * target. Here the firmware provides it in cw_randombytes.c and the host tool in
 * conj_vector_host.c. */
#ifndef RANDOMBYTES_H
#define RANDOMBYTES_H

int randombytes(unsigned char *x, unsigned long long xlen);
int randombytes_select(unsigned char *x, unsigned long long xlen);
void randombytes_init(unsigned char *entropy_input, unsigned char *personalization_string, int security_strength);

#endif
