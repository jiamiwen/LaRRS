#ifndef API_H
#define API_H

#include "config.h"

#if LARRS_MODE == 2
#define CRYPTO_PUBLICKEYBYTES 1312
#define CRYPTO_SECRETKEYBYTES 2544
#define CRYPTO_BYTES 2420

#elif LARRS_MODE == 1
#define CRYPTO_PUBLICKEYBYTES 1952
#define CRYPTO_SECRETKEYBYTES 4016
#define CRYPTO_BYTES 3293

#endif

#define crypto_sign_keypair LARRS_NAMESPACE(_keypair)
int crypto_sign_keypair(unsigned char *pk, unsigned char *sk);

#define crypto_sign LARRS_NAMESPACE()
int crypto_sign(unsigned char *sm, unsigned long long *smlen,
                const unsigned char *msg, unsigned long long len,
                const unsigned char *sk);

#define crypto_sign_open LARRS_NAMESPACE(_open)
int crypto_sign_open(unsigned char *m, unsigned long long *mlen,
                     const unsigned char *sm, unsigned long long smlen,
                     const unsigned char *pk);

#endif
