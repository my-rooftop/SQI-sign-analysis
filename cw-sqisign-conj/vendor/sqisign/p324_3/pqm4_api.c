// SPDX-License-Identifier: Apache-2.0

#include <api.h>
#include <encoded_sizes.h>
#include <sig.h>

// The CRYPTO_* macros in api.h and the sizes in encoded_sizes.h are generated separately. Callers allocate by
// the former, while the library checks lengths against the latter, so a mismatch would mean undersized buffers
// at runtime.
_Static_assert(CRYPTO_BYTES == SIGNATURE_BYTES, "api.h and encoded_sizes.h disagree on the signature size");
_Static_assert(CRYPTO_PUBLICKEYBYTES == PUBLICKEY_BYTES,
               "api.h and encoded_sizes.h disagree on the public key size");
_Static_assert(CRYPTO_SECRETKEYBYTES == SECRETKEY_BYTES,
               "api.h and encoded_sizes.h disagree on the secret key size");

int crypto_sign_keypair(unsigned char *pk, unsigned char *sk) {
  return sqisign_keypair(pk, sk);
}

int crypto_sign(unsigned char *sm, size_t *smlen, const unsigned char *m,
                size_t mlen, const unsigned char *sk) {
  unsigned long long smlen_ull = 0;
  int ret = sqisign_sign(sm, &smlen_ull, m, mlen, sk);
  if (smlen) {
    *smlen = smlen_ull;
  }
  return ret;
}

int crypto_sign_open(unsigned char *m, size_t *mlen, const unsigned char *sm,
                     size_t smlen, const unsigned char *pk) {
  unsigned long long mlen_ull = 0;
  int ret = sqisign_open(m, &mlen_ull, sm, smlen, pk);
  if (mlen) {
    *mlen = mlen_ull;
  }
  return ret;
}
