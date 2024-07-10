#ifndef POLYVEC_H
#define POLYVEC_H

#include <stdint.h>
#include "params.h"
#include "poly.h"

/* Vectors of polynomials of length L */
typedef struct {
  poly vec[V];
} polyvecv;

#define polyvecv_uniform_beta LARRS_NAMESPACE(_polyvecv_uniform_beta)
void polyvecv_uniform_beta(polyvecv *v, const uint8_t seed[SEEDBYTES], uint16_t nonce);

#define polyvecv_uniform_eta LARRS_NAMESPACE(_polyvecv_uniform_eta)
void polyvecv_uniform_eta(polyvecv *v, const uint8_t seed[SEEDBYTES], uint16_t nonce);

#define polyvecv_uniform_gamma LARRS_NAMESPACE(_polyvecv_uniform_gamma)
void polyvecv_uniform_gamma(polyvecv *v, const uint8_t seed[SEEDBYTES], uint16_t nonce);

#define polyvecv_reduce LARRS_NAMESPACE(_polyvecv_reduce)
void polyvecv_reduce(polyvecv *v);

#define polyvecv_freeze LARRS_NAMESPACE(_polyvecv_freeze)
void polyvecv_freeze(polyvecv *v);

#define polyvecv_add LARRS_NAMESPACE(_polyvecv_add)
void polyvecv_add(polyvecv *w, const polyvecv *u, const polyvecv *v);
#define polyvecv_sub LARRS_NAMESPACE(_polyvecv_sub)
void polyvecv_sub(polyvecv *w, const polyvecv *u, const polyvecv *v);

#define polyvecv_ntt LARRS_NAMESPACE(_polyvecv_ntt)
void polyvecv_ntt(polyvecv *v);
#define polyvecv_invntt_tomont LARRS_NAMESPACE(_polyvecv_invntt_tomont)
void polyvecv_invntt_tomont(polyvecv *v);
#define polyvecv_pointwise_poly_montgomery LARRS_NAMESPACE(_polyvecv_pointwise_poly_montgomery)
void polyvecv_pointwise_poly_montgomery(polyvecv *r, const poly *a, const polyvecv *v);
#define polyvecv_pointwise_acc_montgomery \
        LARRS_NAMESPACE(_polyvecv_pointwise_acc_montgomery)
void polyvecv_pointwise_acc_montgomery(poly *w,
                                       const polyvecv *u,
                                       const polyvecv *v);


#define polyvecv_chknorm LARRS_NAMESPACE(_polyvecv_chknorm)
int polyvecv_chknorm(const polyvecv *v, int32_t B);



/* Vectors of polynomials of length K */
typedef struct {
  poly vec[H];
} polyvech;

#define polyvech_uniform_beta LARRS_NAMESPACE(_polyvech_uniform_beta)
void polyvech_uniform_beta(polyvech *v, const uint8_t seed[SEEDBYTES], uint16_t nonce);
#define polyvech_uniform_eta LARRS_NAMESPACE(_polyvech_uniform_eta)
void polyvech_uniform_eta(polyvech *v, const uint8_t seed[SEEDBYTES], uint16_t nonce);
#define polyvech_uniform_gamma LARRS_NAMESPACE(_polyvech_uniform_gamma)
void polyvech_uniform_gamma(polyvech *v, const uint8_t seed[SEEDBYTES], uint16_t nonce);

#define polyvech_reduce LARRS_NAMESPACE(_polyvech_reduce)
void polyvech_reduce(polyvech *v);
#define polyvech_caddq LARRS_NAMESPACE(_polyvech_caddq)
void polyvech_caddq(polyvech *v);
#define polyvech_freeze LARRS_NAMESPACE(_polyvech_freeze)
void polyvech_freeze(polyvech *v);

#define polyvech_add LARRS_NAMESPACE(_polyvech_add)
void polyvech_add(polyvech *w, const polyvech *u, const polyvech *v);
#define polyvech_sub LARRS_NAMESPACE(_polyvech_sub)
void polyvech_sub(polyvech *w, const polyvech *u, const polyvech *v);

#define polyvech_ntt LARRS_NAMESPACE(_polyvech_ntt)
void polyvech_ntt(polyvech *v);
#define polyvech_invntt_tomont LARRS_NAMESPACE(_polyvech_invntt_tomont)
void polyvech_invntt_tomont(polyvech *v);
#define polyvech_pointwise_poly_montgomery LARRS_NAMESPACE(_polyvech_pointwise_poly_montgomery)
void polyvech_pointwise_poly_montgomery(polyvech *r, const poly *a, const polyvech *v);

#define polyvech_pointwise_acc_montgomery \
        LARRS_NAMESPACE(_polyvech_pointwise_acc_montgomery)
void polyvech_pointwise_acc_montgomery(poly *w,
                                       const polyvech *u,
                                       const polyvech *v);
#define polyvech_chknorm LARRS_NAMESPACE(_polyvech_chknorm)
int polyvech_chknorm(const polyvech *v, int32_t B);


#define polyvec_matrix_expand LARRS_NAMESPACE(_polyvec_matrix_expand)
void polyvec_matrix_expand(polyvecv mat[H], const uint8_t rho[SEEDBYTES]);

#define polyvec_matrix_pointwise_montgomery LARRS_NAMESPACE(_polyvec_matrix_pointwise_montgomery)
void polyvec_matrix_pointwise_montgomery(polyvech *t, const polyvecv mat[H], const polyvecv *v);

#define polyvec_matrix_trans_pointwise_montgomery LARRS_NAMESPACE(_polyvec_matrix_trans_pointwise_montgomery)
void polyvec_matrix_trans_pointwise_montgomery(polyvecv *t, const polyvecv mat[H], const polyvech *v);
#endif
