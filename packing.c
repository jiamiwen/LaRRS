#include "params.h"
#include "packing.h"
#include "polyvec.h"
#include "poly.h"

/*************************************************
* Name:        packy
*
* Description: Bit-pack part of public key
*
* Arguments:   - uint8_t pk[]: output byte array
*              - const polyvech *y: pointer to vector yi
*              - int TA: if 1, do not generate id
**************************************************/
void pack_y(uint8_t pk_y[POLYQ_PACKEDBYTES],
             const polyvech *y) {

    unsigned int i;
    for (i = 0; i < H; ++i)
        polyq_pack(pk_y + i * POLYQ_PACKEDBYTES, &y->vec[i]);


}

/*************************************************
* Name:        unpack_y
*
* Description: Unpack public key pk = (yi, H'(yi),).
*
*
* Arguments:
*              - const polyvech *y: pointer to output vector y
*              - uint8_t pk[]: byte array containing bit-packed pk
**************************************************/
void unpack_y(polyvech *y,
               const uint8_t pk[CRYPTO_PUBLICKEYBYTES]) {
    unsigned int i;

    for (i = 0; i < H; ++i)
        polyq_unpack(&y->vec[i], pk + i * POLYQ_PACKEDBYTES);

}

//void unpack_TApk(polyvech *y,
//                 const uint8_t pk[CRYPTO_PUBLICKEYBYTES]){
//    unsigned int i;
//
//    for (i = 0; i < H; ++i)
//        polybeta_unpack(&y->vec[i], pk + i * POLYBETA_PACKEDBYTES);
//}

/*************************************************
* Name:        unpack_pk1
*
* Description: Unpack public key pk = (yi).
*
* Arguments:
*              - const polyvech *y: pointer to output vector y
*              - uint8_t pk[]: byte array containing bit-packed pk
**************************************************/
//void unpack_pk1(polyvech *y,
//               const uint8_t pk[CRYPTO_PUBLICKEYBYTES]) {
//    unsigned int i;
//
//    for (i = 0; i < H; ++i)
//        polybeta_unpack(&y->vec[i], pk + i * POLYBETA_PACKEDBYTES);
//
//}

/*************************************************
* Name:        pack_sk
*
* Description: Bit-pack secret key sk = (x, x').
*
* Arguments:   - uint8_t sk[]: output byte array
*              - const polyvecv *x: pointer to vector x
*              - const polyvech *xprime: pointer to vector x'
**************************************************/
void pack_sk(uint8_t sk[CRYPTO_SECRETKEYBYTES],
             const polyvecv *x,
             const polyvech *xprime) {
    unsigned int i;

    for (i = 0; i < V; ++i)
        polybeta_pack(sk + i * POLYBETA_PACKEDBYTES, &x->vec[i]);
    sk += V * POLYBETA_PACKEDBYTES;

    for (i = 0; i < H; ++i)
        polybeta_pack(sk + i * POLYBETA_PACKEDBYTES, &xprime->vec[i]);

}

/*************************************************
* Name:        unpack_sk
*
* Description: Unpack secret key sk = (x,xprime).
*
* Arguments:
*              - const polyvecv *x: pointer to output vector x
*              - const polyvech *xprime: pointer to output vector xprime
*              - uint8_t sk[]: byte array containing bit-packed sk
**************************************************/
void unpack_sk(polyvecv *x,
               polyvech *xprime,
               const uint8_t sk[CRYPTO_SECRETKEYBYTES]) {
    unsigned int i;

    for (i = 0; i < V; ++i)
        polybeta_unpack(&x->vec[i], sk + i * POLYBETA_PACKEDBYTES);
    sk += V * POLYBETA_PACKEDBYTES;

    for (i = 0; i < H; ++i)
        polybeta_unpack(&xprime->vec[i], sk + i * POLYBETA_PACKEDBYTES);
}

/*************************************************
* Name:        pack_sig
*
* Description: Bit-pack signature sig = (c, z, h).
*
* Arguments:   - uint8_t sig[]: output byte array
*              - const uint8_t *c: pointer to challenge hash length SEEDBYTES
*              - const polyvecv *z: pointer to vector z
*              - const polyvech *h: pointer to hint vector h
**************************************************/
//void pack_sig(uint8_t sig[CRYPTO_BYTES],
//              const uint8_t c[SEEDBYTES],
//              const polyvecv *z,
//              const polyvech *h)
//{
//  unsigned int i, j, k;
//
//  for(i=0; i < SEEDBYTES; ++i)
//    sig[i] = c[i];
//  sig += SEEDBYTES;
//
//  for(i = 0; i < L; ++i)
//    polyz_pack(sig + i*POLYZ_PACKEDBYTES, &z->vec[i]);
//  sig += L*POLYZ_PACKEDBYTES;
//
//  /* Encode h */
//  for(i = 0; i < OMEGA + K; ++i)
//    sig[i] = 0;
//
//  k = 0;
//  for(i = 0; i < K; ++i) {
//    for(j = 0; j < N; ++j)
//      if(h->vec[i].coeffs[j] != 0)
//        sig[k++] = j;
//
//    sig[OMEGA + i] = k;
//  }
//}

/*************************************************
* Name:        unpack_sig
*
* Description: Unpack signature sig = (c, z, h).
*
* Arguments:   - uint8_t *c: pointer to output challenge hash
*              - polyvecv *z: pointer to output vector z
*              - polyvech *h: pointer to output hint vector h
*              - const uint8_t sig[]: byte array containing
*                bit-packed signature
*
* Returns 1 in case of malformed signature; otherwise 0.
**************************************************/
//int unpack_sig(uint8_t c[SEEDBYTES],
//               polyvecv *z,
//               polyvech *h,
//               const uint8_t sig[CRYPTO_BYTES])
//{
//  unsigned int i, j, k;
//
//  for(i = 0; i < SEEDBYTES; ++i)
//    c[i] = sig[i];
//  sig += SEEDBYTES;
//
//  for(i = 0; i < L; ++i)
//    polyz_unpack(&z->vec[i], sig + i*POLYZ_PACKEDBYTES);
//  sig += L*POLYZ_PACKEDBYTES;
//
//  /* Decode h */
//  k = 0;
//  for(i = 0; i < K; ++i) {
//    for(j = 0; j < N; ++j)
//      h->vec[i].coeffs[j] = 0;
//
//    if(sig[OMEGA + i] < k || sig[OMEGA + i] > OMEGA)
//      return 1;
//
//    for(j = k; j < sig[OMEGA + i]; ++j) {
//      /* Coefficients are ordered for strong unforgeability */
//      if(j > k && sig[j] <= sig[j-1]) return 1;
//      h->vec[i].coeffs[sig[j]] = 1;
//    }
//
//    k = sig[OMEGA + i];
//  }
//
//  /* Extra indices are zero for strong unforgeability */
//  for(j = k; j < OMEGA; ++j)
//    if(sig[j])
//      return 1;
//
//  return 0;
//}
/*************************************************
* Name:        pack_ctag
*
* Description: Bit-pack signature ctag = (c1,c2) \in R_q^v \times R_q.
*
* Arguments:   - uint8_t ctag[]: output byte array
*              - const polyvecv *c1: pointer to vector c1
*              - const poly *c2: pointer to hint vector c2
**************************************************/
void pack_ctag(uint8_t ctag[CTAG_BYTES],
               const polyvecv *c1,
               const poly *c2) {
    unsigned int i;
    for (i = 0; i < V; ++i) {
        polyq_pack(ctag + i * POLYQ_PACKEDBYTES, &c1->vec[i]);
    }
    ctag += V * POLYQ_PACKEDBYTES;

    polyq_pack(ctag, c2);

}

/*************************************************
* Name:        unpack_ctag
*
* Description: Unpack revoke tag ctag = (c1, c2).
*
* Arguments:
*              - polyvecv *c1: pointer to output vector c1
*              - poly *c2: pointer to output vector c2
*              - const uint8_t ctag[]: byte array containing
*                bit-packed tag
*
**************************************************/
void unpack_ctag(polyvecv *c1,
                 poly *c2,
                 const uint8_t ctag[CTAG_BYTES]) {
    unsigned int i;
    for (i = 0; i < V; ++i) {
        polyq_unpack(&c1->vec[i], ctag + i * POLYQ_PACKEDBYTES);
    }
    ctag += V * POLYQ_PACKEDBYTES;

    polyq_unpack(c2, ctag);
}

/*************************************************
* Name:        pack_resp
*
* Description: Bit-pack responses resp = (r_i,t_i,rho_i) \in S_{\eta}^{v+h}
*              \times S_{\eta}^h  \times S_{\BETA * \KAPPA}
*
*
* Arguments:   - uint8_t resp[]: output byte array
*              - const polyvecv *r1: pointer to vector r1
*              - const polyvech *r2: pointer to vector r2
*              - const polyvech *t: pointer to vector t
*              - const poly *rho: pointer to vector rho
**************************************************/
void pack_respR(uint8_t resp[R_BYTES],
                const polyvecv *r1,
                const polyvech *r2) {
    unsigned int i;
    for (i = 0; i < V; ++i) {
        polyeta_pack(resp + i * POLYETA_PACKEDBYTES, &r1->vec[i]);
    }
    resp += V * POLYETA_PACKEDBYTES;

    for (i = 0; i < H; ++i) {
        polyeta_pack(resp + i * POLYETA_PACKEDBYTES, &r2->vec[i]);
    }
}
void pack_respT(uint8_t resp[T_BYTES],
                const polyvech *t) {
    unsigned int i;
    for (i = 0; i < H; ++i) {
        polyeta_pack(resp + i * POLYETA_PACKEDBYTES, &t->vec[i]);
    }
}

void pack_respRho(uint8_t resp[POLYRHO_PACKEDBYTES],
                  const poly *rho){
    polyrho_pack(resp, rho);

}

/*************************************************
* Name:        unpack_resp
*
* Description: Unpack responses resp = (r_i,t_i,rho_i).
*
* Arguments:
*              - polyvecv *r1: pointer to output vector r1
*              - polyvech *r2: pointer to output vector r2
*              - polyvecv *t: pointer to output vector t
*              - poly *rho: pointer to output vector rho
*              - const uint8_t resp[]: byte array containing
*                bit-packed tag
*
**************************************************/
void unpack_respR(polyvecv *r1,
                  polyvech *r2,
                  const uint8_t resp[RESPONSE_BYTES]) {
    unsigned int i;
    for (i = 0; i < V; ++i) {
        polyeta_unpack(&r1->vec[i], resp + i * POLYETA_PACKEDBYTES);
    }
    resp += V * POLYETA_PACKEDBYTES;

    for (i = 0; i < H; ++i) {
        polyeta_unpack(&r2->vec[i], resp + i * POLYETA_PACKEDBYTES);
    }
}
void unpack_respT(polyvech *t,
                  const uint8_t resp[RESPONSE_BYTES]) {
    unsigned int i;
    for (i = 0; i < H; ++i) {
        polyeta_unpack(&t->vec[i], resp + i * POLYETA_PACKEDBYTES);
    }
}
void unpack_respRho(poly *rho,
                    const uint8_t resp[RESPONSE_BYTES]){

    polyrho_unpack(rho, resp);
}

/*************************************************
* Name:        pack_e1
*
* Description: Bit-pack bytes
*
*
* Arguments:   - uint8_t sig[]: output byte array
*              - const uint8_t e1_seed[]: byte array of e1 seed
**************************************************/
void pack_bytes(uint8_t *sig,
                const uint8_t *bytes, unsigned int len) {
    unsigned int i;
    for (i = 0; i < len; ++i) {
        sig[i] = bytes[i];
    }
}

/*************************************************
* Name:        unpack_bytes
*
* Description: Bit-pack e1 seed
*
*
* Arguments:   - uint8_t e1_seed[]: output byte array
*              - const uint8_t sig[]: byte array of sig
**************************************************/
void unpack_bytes(uint8_t *bytes,
                  const uint8_t *sig, unsigned int len) {
    unsigned int i;
    for (i = 0; i < len; ++i) {
        bytes[i] = sig[i];
    }
}

void pack_epsilon(uint8_t sig[EPSILON_BYTES],
                  const polyvecv *epsilon){

    unsigned int i;
    for(i = 0; i < V; ++i){
        polyq_pack(sig + i * POLYQ_PACKEDBYTES, &epsilon->vec[i]);
    }
}

void unpack_epsilon(polyvecv *epsilon,
                    const uint8_t sig[EPSILON_BYTES]){

    unsigned int i;
    for(i = 0; i < V; ++i){
        polyq_unpack( &epsilon->vec[i], sig + i * POLYQ_PACKEDBYTES);
    }
}
