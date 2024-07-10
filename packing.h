#ifndef PACKING_H
#define PACKING_H

#include <stdint.h>
#include "params.h"
#include "polyvec.h"
#include "fips202.h"

#define pack_y LARRS_NAMESPACE(_pack_y)
void pack_y(uint8_t pk[POLYQ_PACKEDBYTES],
             const polyvech *t1);


#define pack_sk LARRS_NAMESPACE(_pack_sk)
void pack_sk(uint8_t sk[CRYPTO_SECRETKEYBYTES],
             const polyvecv *x,
             const polyvech *xprime);

#define pack_ctag LARRS_NAMESPACE(_pack_ctag)
void pack_ctag(uint8_t ctag[CTAG_BYTES],
               const polyvecv *c1,
               const poly *c2);

#define pack_epsilon LARRS_NAMESPACE(_pack_epsilon)
void pack_epsilon(uint8_t sig[EPSILON_BYTES],
                  const polyvecv *epsilon);

#define pack_respR LARRS_NAMESPACE(_pack_respR)
void pack_respR(uint8_t resp[R_BYTES],
               const polyvecv *r1,
               const polyvech *r2);
#define pack_respT LARRS_NAMESPACE(_pack_respT)
void pack_respT(uint8_t resp[T_BYTES],
                const polyvech *t);
#define pack_respRho LARRS_NAMESPACE(_pack_respRho)
void pack_respRho(uint8_t resp[POLYRHO_PACKEDBYTES],
                const poly *rho);

#define pack_bytes LARRS_NAMESPACE(_pack_bytes)
void pack_bytes(uint8_t *sig,
             const uint8_t *bytes, unsigned int len);

//#define pack_sig LARRS_NAMESPACE(_pack_sig)
//void pack_sig(uint8_t sig[CRYPTO_BYTES],
//              const uint8_t c[SEEDBYTES], const polyvecv *z, const polyvech *h);

#define unpack_y LARRS_NAMESPACE(_unpack_y)
void unpack_y(polyvech *t1,
               const uint8_t pk[CRYPTO_PUBLICKEYBYTES]);
//#define unpack_TApk LARRS_NAMESPACE(_unpack_TApk)
//void unpack_TApk(polyvech *t1,
//                 const uint8_t pk[CRYPTO_PUBLICKEYBYTES]);
//
//#define unpack_pk1 LARRS_NAMESPACE(_unpack_pk1)
//void unpack_pk1(polyvech *t1,
//               const uint8_t pk[CRYPTO_PUBLICKEYBYTES]);

#define unpack_sk LARRS_NAMESPACE(_unpack_sk)
void unpack_sk(polyvecv *s1,
               polyvech *s2,
               const uint8_t sk[CRYPTO_SECRETKEYBYTES]);
#define unpack_ctag LARRS_NAMESPACE(_unpack_ctag)
void unpack_ctag(polyvecv *c1,
                poly *c2,
                const uint8_t ctag[CTAG_BYTES]);

#define unpack_epsilon LARRS_NAMESPACE(_unpack_epsilon)
void unpack_epsilon(polyvecv *epsilon,
                    const uint8_t sig[EPSILON_BYTES]);

#define unpack_respR LARRS_NAMESPACE(_unpack_respR)
void unpack_respR(polyvecv *r1,
               polyvech *r2,
               const uint8_t resp[RESPONSE_BYTES]);
#define unpack_respT LARRS_NAMESPACE(_unpack_respT)
void unpack_respT(polyvech *t,
                 const uint8_t resp[RESPONSE_BYTES]);
#define unpack_respRho LARRS_NAMESPACE(_unpack_respRho)
void unpack_respRho(poly *rho,
                 const uint8_t resp[RESPONSE_BYTES]);

#define unpack_bytes LARRS_NAMESPACE(_unpack_bytes)
void unpack_bytes(uint8_t *bytes,
             const uint8_t *sig, unsigned int len);
//#define unpack_sig LARRS_NAMESPACE(_unpack_sig)
//int unpack_sig(uint8_t c[SEEDBYTES], polyvecv *z, polyvech *h,
//               const uint8_t sig[CRYPTO_BYTES]);


#endif
