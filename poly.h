#ifndef POLY_H
#define POLY_H

#include <stdint.h>
#include "params.h"


typedef struct {
  int32_t coeffs[N];
} poly;


#define poly_reduce LARRS_NAMESPACE(_poly_reduce)
void poly_reduce(poly *a);
#define poly_caddq LARRS_NAMESPACE(_poly_caddq)
void poly_caddq(poly *a);
#define poly_freeze LARRS_NAMESPACE(_poly_freeze)
void poly_freeze(poly *a);

#define poly_add LARRS_NAMESPACE(_poly_add)
void poly_add(poly *c, const poly *a, const poly *b);
#define poly_sub LARRS_NAMESPACE(_poly_sub)
void poly_sub(poly *c, const poly *a, const poly *b);

#define poly_cmp LARRS_NAMESPACE(_poly_cmp)
int poly_cmp(poly *a, poly *b);

#define poly_ntt LARRS_NAMESPACE(_poly_ntt)
void poly_ntt(poly *a);
#define poly_invntt_tomont LARRS_NAMESPACE(_poly_invntt_tomont)
void poly_invntt_tomont(poly *a);
#define poly_pointwise_montgomery LARRS_NAMESPACE(_poly_pointwise_montgomery)
void poly_pointwise_montgomery(poly *c, const poly *a, const poly *b);

#define poly_inv LARRS_NAMESPACE(_poly_inv)
int poly_inv(poly *a);

#define poly_chknorm LARRS_NAMESPACE(_poly_chknorm)
int poly_chknorm(const poly *a, int32_t B);
#define poly_uniform LARRS_NAMESPACE(_poly_uniform)
void poly_uniform(poly *a,
                  const uint8_t seed[SEEDBYTES],
                  uint16_t nonce);
#define poly_uniform_beta LARRS_NAMESPACE(_poly_uniform_beta)
void poly_uniform_beta(poly *a,
                      const uint8_t seed[SEEDBYTES],
                      uint16_t nonce);

#define poly_uniform_gamma LARRS_NAMESPACE(_poly_uniform_gamma)
void poly_uniform_gamma(poly *a,
                         const uint8_t seed[CRHBYTES],
                         uint16_t nonce);
#define poly_uniform_eta LARRS_NAMESPACE(_poly_uniform_eta)
void poly_uniform_eta(poly*a, const uint8_t seed[SEEDBYTES], uint16_t nonce);

//#define poly_gen_epsilon LARRS_NAMESPACE(_poly_gen_epsilon)
//uint16_t poly_gen_epsilon2(poly *epsilon2,
//                           poly *i_epsilon2,
//                           const uint8_t seed[SEEDBYTES],
//                           uint16_t nonce);

#define poly_challenge LARRS_NAMESPACE(_poly_challenge)
void poly_challenge(poly *c, const uint8_t seed[SEEDBYTES]);

#define poly_frommsg LARRS_NAMESPACE(_poly_frommsg)
void poly_frommsg(poly *a, const uint8_t msg[SEEDBYTES]);
#define poly_round LARRS_NAMESPACE(_poly_round)
void poly_round(poly *a);

#define polybeta_pack LARRS_NAMESPACE(_polybeta_pack)
void polybeta_pack(uint8_t *r, const poly *a);
#define polybeta_unpack LARRS_NAMESPACE(_polybeta_unpack)
void polybeta_unpack(poly *r, const uint8_t *a);

#define polyeta_pack LARRS_NAMESPACE(_polyeta_pack)
void polyeta_pack(uint8_t *r, const poly *a);
#define polyeta_unpack LARRS_NAMESPACE(_polyeta_unpack)
void polyeta_unpack(poly *r, const uint8_t *a);

#define polyq_pack LARRS_NAMESPACE(_polyq_pack)
void polyq_pack(uint8_t *r, const poly *a);
#define polyq_unpack LARRS_NAMESPACE(_polyq_unpack)
void polyq_unpack(poly *r, const uint8_t *a);

#define polyrho_pack LARRS_NAMESPACE(_polyrho_pack)
void polyrho_pack(uint8_t *r, const poly *a);
#define polyrho_unpack LARRS_NAMESPACE(_polyrho_unpack)
void polyrho_unpack(poly *r, const uint8_t *a);
#endif
