#include <stdint.h>
#include "params.h"
#include "polyvec.h"
#include "poly.h"

/*************************************************
* Name:        expand_mat
*
* Description: Implementation of ExpandA. Generates matrix A with uniformly
*              random coefficients a_{i,j} by performing rejection
*              sampling on the output stream of SHAKE128(rho|j|i)
*              or AES256CTR(rho,j|i).
*
* Arguments:   - polyvecv mat[K]: output matrix
*              - const uint8_t rho[]: byte array containing seed rho
**************************************************/
void polyvec_matrix_expand(polyvecv mat[H], const uint8_t rho[SEEDBYTES]) {
  unsigned int i, j;

  for(i = 0; i < H; ++i)
    for(j = 0; j < V; ++j)
      poly_uniform(&mat[i].vec[j], rho, (i << 8) + j);
}

void polyvec_matrix_pointwise_montgomery(polyvech *t, const polyvecv mat[H], const polyvecv *v) {
  unsigned int i;

  for(i = 0; i < H; ++i)
    polyvecv_pointwise_acc_montgomery(&t->vec[i], &mat[i], v);
}

/*************************************************
* Name:        polyvec_matrix_trans_pointwise_montgomery
*
* Description: compute t = mat^T * v
*
* Arguments:   - poly *t: output polynomial
*              - const polyvecv mat[H]: pointer to first input matrix
*              - const polyvech *v: pointer to second input vector
**************************************************/
void polyvec_matrix_trans_pointwise_montgomery(polyvecv *t, const polyvecv mat[H], const polyvech *v){
    unsigned int i, j;
    poly w;
    for(i = 0; i < V; ++i){
        poly_pointwise_montgomery(&t->vec[i], &mat[0].vec[i], &v->vec[0]);
        for(j = 1; j < H; ++j){
            poly_pointwise_montgomery(&w, &mat[j].vec[i], &v->vec[j]);
            poly_add(&t->vec[i], &t->vec[i], &w);
        }
    }
}

/**************************************************************/
/************ Vectors of polynomials of length V **************/
/**************************************************************/

void polyvecv_uniform_beta(polyvecv *v, const uint8_t seed[SEEDBYTES], uint16_t nonce) {
  unsigned int i;

  for(i = 0; i < V; ++i)
    poly_uniform_beta(&v->vec[i], seed, nonce++);
}

void polyvecv_uniform_gamma(polyvecv *v, const uint8_t seed[SEEDBYTES], uint16_t nonce) {
  unsigned int i;

  for(i = 0; i < V; ++i)
    poly_uniform_gamma(&v->vec[i], seed, nonce++);
}

void polyvecv_uniform_eta(polyvecv *v, const uint8_t seed[SEEDBYTES], uint16_t nonce){
    unsigned int i;
    for(i = 0; i < V; ++i){
        poly_uniform_eta(&v->vec[i], seed, nonce++);
    }
}

void polyvecv_reduce(polyvecv *v) {
  unsigned int i;

  for(i = 0; i < V; ++i)
    poly_reduce(&v->vec[i]);
}

/*************************************************
* Name:        polyvecv_freeze
*
* Description: Reduce coefficients of polynomials in vector of length V
*              to standard representatives.
*
* Arguments:   - polyvecv *v: pointer to input/output vector
**************************************************/
void polyvecv_freeze(polyvecv *v) {
  unsigned int i;

  for(i = 0; i < V; ++i)
    poly_freeze(&v->vec[i]);
}

/*************************************************
* Name:        polyvecv_add
*
* Description: Add vectors of polynomials of length V.
*              No modular reduction is performed.
*
* Arguments:   - polyvecv *w: pointer to output vector
*              - const polyvecv *u: pointer to first summand
*              - const polyvecv *v: pointer to second summand
**************************************************/
void polyvecv_add(polyvecv *w, const polyvecv *u, const polyvecv *v) {
  unsigned int i;

  for(i = 0; i < V; ++i)
    poly_add(&w->vec[i], &u->vec[i], &v->vec[i]);
}

void polyvecv_sub(polyvecv *w, const polyvecv *u, const polyvecv *v){
    unsigned int i;

    for(i = 0; i < V; ++i)
        poly_sub(&w->vec[i], &u->vec[i], &v->vec[i]);
}

/*************************************************
* Name:        polyvecv_ntt
*
* Description: Forward NTT of all polynomials in vector of length L. Output
*              coefficients can be up to 16*Q larger than input coefficients.
*
* Arguments:   - polyvecv *v: pointer to input/output vector
**************************************************/
void polyvecv_ntt(polyvecv *v) {
  unsigned int i;

  for(i = 0; i < V; ++i)
    poly_ntt(&v->vec[i]);
}

void polyvecv_invntt_tomont(polyvecv *v) {
  unsigned int i;

  for(i = 0; i < V; ++i)
    poly_invntt_tomont(&v->vec[i]);
}

void polyvecv_pointwise_poly_montgomery(polyvecv *r, const poly *a, const polyvecv *v) {
  unsigned int i;

  for(i = 0; i < V; ++i)
    poly_pointwise_montgomery(&r->vec[i], a, &v->vec[i]);
}

/*************************************************
* Name:        polyvecv_pointwise_acc_montgomery
*
* Description: Pointwise multiply vectors of polynomials of length L, multiply
*              resulting vector by 2^{-32} and add (accumulate) polynomials
*              in it. Input/output vectors are in NTT domain representation.
*
* Arguments:   - poly *w: output polynomial
*              - const polyvecv *u: pointer to first input vector
*              - const polyvecv *v: pointer to second input vector
**************************************************/
void polyvecv_pointwise_acc_montgomery(poly *w,
                                       const polyvecv *u,
                                       const polyvecv *v)
{
  unsigned int i;
  poly t;

  poly_pointwise_montgomery(w, &u->vec[0], &v->vec[0]);
  for(i = 1; i < V; ++i) {
    poly_pointwise_montgomery(&t, &u->vec[i], &v->vec[i]);
    poly_add(w, w, &t);
  }
}

/*************************************************
* Name:        polyvecv_chknorm
*
* Description: Check infinity norm of polynomials in vector of length L.
*              Assumes input polyvecv to be reduced by polyvecv_reduce().
*
* Arguments:   - const polyvecv *v: pointer to vector
*              - int32_t B: norm bound
*
* Returns 0 if norm of all polynomials is strictly smaller than B <= (Q-1)/8
* and 1 otherwise.
**************************************************/
int polyvecv_chknorm(const polyvecv *v, int32_t bound)  {
  unsigned int i;

  for(i = 0; i < V; ++i)
    if(poly_chknorm(&v->vec[i], bound))
      return 1;

  return 0;
}

/**************************************************************/
/************ Vectors of polynomials of length K **************/
/**************************************************************/

void polyvech_uniform_beta(polyvech *v, const uint8_t seed[SEEDBYTES], uint16_t nonce) {
  unsigned int i;

  for(i = 0; i < H; ++i)
    poly_uniform_beta(&v->vec[i], seed, nonce++);
}

void polyvech_uniform_gamma(polyvech *v, const uint8_t seed[SEEDBYTES], uint16_t nonce) {
    unsigned int i;

    for(i = 0; i < H; ++i)
        poly_uniform_gamma(&v->vec[i], seed, nonce++);
}

void polyvech_uniform_eta(polyvech *v, const uint8_t seed[SEEDBYTES], uint16_t nonce){
    unsigned int i;
    for(i = 0; i < H; ++i){
        poly_uniform_eta(&v->vec[i], seed, nonce++);
    }
}


/*************************************************
* Name:        polyvech_reduce
*
* Description: Reduce coefficients of polynomials in vector of length H
*              to representatives in [-6283009,6283007].
*
* Arguments:   - polyvech *v: pointer to input/output vector
**************************************************/
void polyvech_reduce(polyvech *v) {
  unsigned int i;

  for(i = 0; i < H; ++i)
    poly_reduce(&v->vec[i]);
}

/*************************************************
* Name:        polyvech_caddq
*
* Description: For all coefficients of polynomials in vector of length H
*              add Q if coefficient is negative.
*
* Arguments:   - polyvech *v: pointer to input/output vector
**************************************************/
void polyvech_caddq(polyvech *v) {
  unsigned int i;

  for(i = 0; i < H; ++i)
    poly_caddq(&v->vec[i]);
}

/*************************************************
* Name:        polyvech_freeze
*
* Description: Reduce coefficients of polynomials in vector of length K
*              to standard representatives.
*
* Arguments:   - polyvech *v: pointer to input/output vector
**************************************************/
void polyvech_freeze(polyvech *v)  {
  unsigned int i;

  for(i = 0; i < H; ++i)
    poly_freeze(&v->vec[i]);
}

/*************************************************
* Name:        polyvech_add
*
* Description: Add vectors of polynomials of length H.
*              No modular reduction is performed.
*
* Arguments:   - polyvech *w: pointer to output vector
*              - const polyvech *u: pointer to first summand
*              - const polyvech *v: pointer to second summand
**************************************************/
void polyvech_add(polyvech *w, const polyvech *u, const polyvech *v) {
  unsigned int i;

  for(i = 0; i < H; ++i)
    poly_add(&w->vec[i], &u->vec[i], &v->vec[i]);
}

/*************************************************
* Name:        polyvech_sub
*
* Description: Subtract vectors of polynomials of length H.
*              No modular reduction is performed.
*
* Arguments:   - polyvech *w: pointer to output vector
*              - const polyvech *u: pointer to first input vector
*              - const polyvech *v: pointer to second input vector to be
*                                   subtracted from first input vector
**************************************************/
void polyvech_sub(polyvech *w, const polyvech *u, const polyvech *v) {
  unsigned int i;

  for(i = 0; i < H; ++i)
    poly_sub(&w->vec[i], &u->vec[i], &v->vec[i]);
}


/*************************************************
* Name:        polyvech_ntt
*
* Description: Forward NTT of all polynomials in vector of length H. Output
*              coefficients can be up to 16*Q larger than input coefficients.
*
* Arguments:   - polyvech *v: pointer to input/output vector
**************************************************/
void polyvech_ntt(polyvech *v) {
  unsigned int i;

  for(i = 0; i < H; ++i)
    poly_ntt(&v->vec[i]);
}

/*************************************************
* Name:        polyvech_invntt_tomont
*
* Description: Inverse NTT and multiplication by 2^{32} of polynomials
*              in vector of length H. Input coefficients need to be less
*              than 2*Q.
*
* Arguments:   - polyvech *v: pointer to input/output vector
**************************************************/
void polyvech_invntt_tomont(polyvech *v) {
  unsigned int i;

  for(i = 0; i < H; ++i)
    poly_invntt_tomont(&v->vec[i]);
}

/*************************************************
* Name:        polyvech_pointwise_acc_montgomery
*
* Description: Pointwise multiply vectors of polynomials of length H, multiply
*              resulting vector by 2^{-32} and add (accumulate) polynomials
*              in it. Input/output vectors are in NTT domain representation.
*
* Arguments:   - poly *w: output polynomial
*              - const polyvech *u: pointer to first input vector
*              - const polyvech *v: pointer to second input vector
**************************************************/
void polyvech_pointwise_acc_montgomery(poly *w,
                                       const polyvech *u,
                                       const polyvech *v)
{
    unsigned int i;
    poly t;

    poly_pointwise_montgomery(w, &u->vec[0], &v->vec[0]);
    for(i = 1; i < H; ++i) {
        poly_pointwise_montgomery(&t, &u->vec[i], &v->vec[i]);
        poly_add(w, w, &t);
    }
}

void polyvech_pointwise_poly_montgomery(polyvech *r, const poly *a, const polyvech *v) {
  unsigned int i;

  for(i = 0; i < H; ++i)
    poly_pointwise_montgomery(&r->vec[i], a, &v->vec[i]);
}


/*************************************************
* Name:        polyvech_chknorm
*
* Description: Check infinity norm of polynomials in vector of length K.
*              Assumes input polyvech to be reduced by polyvech_reduce().
*
* Arguments:   - const polyvech *v: pointer to vector
*              - int32_t B: norm bound
*
* Returns 0 if norm of all polynomials are strictly smaller than B <= (Q-1)/8
* and 1 otherwise.
**************************************************/
int polyvech_chknorm(const polyvech *v, int32_t bound) {
  unsigned int i;

  for(i = 0; i < H; ++i)
    if(poly_chknorm(&v->vec[i], bound))
      return 1;

  return 0;
}



