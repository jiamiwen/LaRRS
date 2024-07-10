#include <stdint.h>

#include "params.h"
#include "sign.h"
#include "packing.h"
#include "polyvec.h"
#include "poly.h"
#include "randombytes.h"

#include "fips202.h"


static polyvecv mat[H];

static inline uint64_t cpucycles(void) {
    uint64_t result;

    __asm__ volatile ("rdtsc; shlq $32,%%rdx; orq %%rdx,%%rax"
            : "=a" (result) : : "%rdx");

    return result;
}
//extern uint64_t t0, t1, epoch, o, rep;
//#define TIME_START() t0 = cpucycles();
//#define TIME_PAUSE() t1= cpucycles(); \
//                    epoch += t1 - t0 - o;

/*************************************************
* Name:        crypto_setup
*
* Description: Generates common parameters A.
*
* Arguments:
*
**************************************************/
void crypto_setup(void) {
    uint8_t seedbuf[SEEDBYTES];
    randombytes(seedbuf, SEEDBYTES);
    polyvec_matrix_expand(mat, seedbuf);
//    for (int i = 0; i < H; ++i) {
//        polyvecv_ntt(&mat[i]);
//    }
}

/*************************************************
* Name:        crypto_sign_keypair
*
* Description: Generates and packs public and private key
*
* Arguments:   - uint8_t *pk: pointer to output public key (allocated
*                             array of CRYPTO_PUBLICKEYBYTES bytes)
*              - uint8_t *sk: pointer to output private key (allocated
*                             array of CRYPTO_SECRETKEYBYTES bytes)
*              - uint8_t seedbuf[] : seed, can be replaced by random bytes
*              - int TA: if 1, do not generate id
*
* Returns 0 (success)
**************************************************/
int crypto_sign_keypair(uint8_t *pk_y, uint8_t *pk_id, uint8_t *sk) {
    polyvecv x, xhat;
    polyvech xprime, y;
    keccak_state state;
    uint8_t seedbuf[SEEDBYTES];

    randombytes(seedbuf, SEEDBYTES);

    polyvecv_uniform_beta(&x, seedbuf, 0);
    polyvech_uniform_beta(&xprime, seedbuf, V);

    xhat = x;

    //y = A*x+x'
    polyvecv_ntt(&xhat);
    polyvec_matrix_pointwise_montgomery(&y, mat, &xhat);
    polyvech_reduce(&y);
    polyvech_invntt_tomont(&y);
    polyvech_add(&y, &y, &xprime);

    pack_sk(sk, &x, &xprime);
    pack_y(pk_y, &y);

    //id = H'(y)
    shake256_init(&state);
    shake256_absorb(&state, pk_y, H * POLYQ_PACKEDBYTES);
    shake256_finalize(&state);
    shake256_squeeze(pk_id, SEEDBYTES, &state);

    return 0;
}

/*************************************************
* Name:        crypto_sign_keypair1
*
* Description: if TA, generates and packs y and sk; else, generates and packs y and id
*
* Arguments:   - uint8_t *var: pointer to output public key (allocated
*                             array of CRYPTO_PUBLICKEYBYTES bytes)
*              - uint8_t seedbuf[] : seed, can be replaced by random bytes
*
* Returns 0 (success)
**************************************************/
int crypto_sign_keypair1(uint8_t *pk_y, uint8_t *var, int TA) {
    polyvecv x, xhat;
    polyvech xprime, y;
    uint8_t seedbuf[SEEDBYTES];

    randombytes(seedbuf, SEEDBYTES);

    polyvecv_uniform_beta(&x, seedbuf, 0);
    polyvech_uniform_beta(&xprime, seedbuf, V);

    xhat = x;

    //y = A*x+x'
    polyvecv_ntt(&xhat);
    polyvec_matrix_pointwise_montgomery(&y, mat, &xhat);
    polyvech_reduce(&y);
    polyvech_invntt_tomont(&y);
    polyvech_add(&y, &y, &xprime);

    pack_y(pk_y, &y);
    if(TA){
        pack_sk(var, &x, &xprime);
    }else{
        keccak_state state;
        shake256_init(&state);
        shake256_absorb(&state, pk_y, H * POLYQ_PACKEDBYTES);
        shake256_finalize(&state);
        shake256_squeeze(var, SEEDBYTES, &state);
    }

    return 0;
}
/*************************************************
* Name:        crypto_sign_signature
*
* Description: Computes signature.
*
* Arguments:   - uint8_t *sig:   pointer to output signature (of length CRYPTO_BYTES)
*                                sig = epsilon||c1||c2||T||e1||r1||t1||rho1||...||rn||tn||rhon
*              - size_t *siglen: pointer to output length of signature
*              - uint8_t *m:     pointer to message to be signed
*              - size_t mlen:    length of message
*              - uint8_t *sk:    pointer to bit-packed secret key
*
* Returns 0 (success)
**************************************************/
int crypto_sign_signature(uint8_t *sig,
                          size_t *siglen,
                          const uint8_t *m,
                          size_t mlen,
                          const SIGN_IN_PARAMS *inParams) {
    unsigned int i, j;
    uint8_t *pt;
    uint8_t seedbuf[SEEDBYTES], seed_beta[SEEDBYTES], seed_eta[SEEDBYTES], seed_gamma[CRHBYTES];
    uint8_t bytes[POLYQ_PACKEDBYTES];
    uint16_t nonce_beta, nonce_eta, nonce_gamma;
    polyvech s, y_TA, u2, w, xprime,t_vech, t_vech1;
    polyvecv epsilon, c1, u1, x, t_vecv, t_vecv1;
    poly c2, epsilon2, e,t_poly, t_poly1;
    keccak_state state, state0;
    time_t T;
//    rep = 1;
#define PACKQ_ABSORB(p)    polyq_pack(bytes, &p);\
            shake256_absorb(&state, bytes, POLYQ_PACKEDBYTES);

    // generates seed
    nonce_beta = nonce_eta = nonce_gamma = 0;
    randombytes(seedbuf, SEEDBYTES);
    shake256_init(&state);
    shake256_absorb(&state, seedbuf, SEEDBYTES);
    shake256_finalize(&state);
    shake256_squeeze(seed_beta, SEEDBYTES, &state);
    shake256_squeeze(seed_gamma, CRHBYTES, &state);
    shake256_squeeze(seed_eta, SEEDBYTES, &state);

    // 1. revocable tag C
    // sample s->S_beta^h, epsilon1 = t_vecv->S_beta^v
//    TIME_START()
    polyvech_uniform_beta(&s, seed_beta, nonce_beta);
    nonce_beta += H;
    polyvecv_uniform_beta(&t_vecv, seed_beta, nonce_beta);
    nonce_beta += V;
//    TIME_PAUSE()

    //sample invertible polynomial epsilon2, t_poly = epsilon2^-1
gen_inverse:
//    TIME_START()
    poly_uniform_beta(&epsilon2, seed_beta, nonce_beta++);
//    TIME_PAUSE()
    t_poly = epsilon2;
    if(!poly_inv(&t_poly))
        goto gen_inverse;

    //C1=A^T*s + \epsilon_1
    polyvech_ntt(&s);//note: s in NTT
    polyvec_matrix_trans_pointwise_montgomery(&c1, mat, &s);
    polyvecv_reduce(&c1);
    polyvecv_invntt_tomont(&c1);
    polyvecv_add(&c1, &c1, &t_vecv);

    // epsilon = epsilon1 * epsilon2^-1
    polyvecv_ntt(&t_vecv);
    poly_ntt(&t_poly);
    polyvecv_pointwise_poly_montgomery(&epsilon, &t_poly, &t_vecv);
    polyvecv_invntt_tomont(&epsilon);

    //pack epsilon
    pack_epsilon(sig, &epsilon);


    //unpack TA's public key
    unpack_y(&y_TA, inParams->y_TA);

    //c2 = y^T*s + \epsilon_2 +  [q/2]*id_pi
    polyvech_ntt(&y_TA);//note:TA_y in NTT
    polyvech_pointwise_acc_montgomery(&t_poly, &y_TA, &s);
//    poly_reduce(&t_poly);
    poly_invntt_tomont(&t_poly);

    poly_add(&c2, &t_poly, &epsilon2);
    poly_frommsg(&t_poly, inParams->id_list + PI * SEEDBYTES);
    poly_add(&c2, &c2, &t_poly);
    poly_reduce(&c2);

    //pack c1 & c2
    pack_ctag(sig + CTAG_POS, &c1, &c2);


    time(&T);
    pt = (uint8_t *)&T;
    //pack T
    pack_bytes(sig + TIMESTAMP_POS, pt, TIMESTAMP_BYTES);

    //keccak_absorb(T,L,M)
    shake256_init(&state0);
    shake256_absorb(&state0, pt, TIMESTAMP_BYTES);
    shake256_absorb(&state0, inParams->y_list, USER_NUM * Y_BYTES);
    shake256_absorb(&state0, m, mlen);

    // prepare epsion, epsilon2, C1 and sk
    polyvecv_ntt(&epsilon);//note: epsilon in NTT form
    polyvecv_ntt(&c1);//note: C1 in NTT form
    poly_ntt(&epsilon2);//note: epsilon2 in NTT form
    unpack_sk(&x, &xprime, inParams->sk);
    polyvecv_ntt(&x);
    polyvech_ntt(&xprime);

rej:
    // 2. compute the next of user pi in the Signing Ring
    // sample u=(u1,u2) and w
//    TIME_START()
    polyvecv_uniform_gamma(&u1, seed_gamma, nonce_gamma);
    nonce_gamma += V;
    polyvech_uniform_gamma(&u2, seed_gamma, nonce_gamma);
    nonce_gamma += H;
    polyvech_uniform_gamma(&w, seed_gamma, nonce_gamma);
    nonce_gamma += H;
//    TIME_PAUSE()

    // t_vech = Ahat*u = A*u1 + u2
    t_vecv = u1;
    polyvecv_ntt(&t_vecv);
    polyvec_matrix_pointwise_montgomery(&t_vech, mat, &t_vecv);
    polyvech_invntt_tomont(&t_vech);
    polyvech_add(&t_vech, &t_vech, &u2);
    polyvech_reduce(&t_vech);

    state = state0;
    //absorb Ahat*u
    for(i = 0; i < H; ++i){
        PACKQ_ABSORB(t_vech.vec[i])
    }

    //t_vecv = A^T * w
    t_vech = w;
    polyvech_ntt(&t_vech);
    polyvec_matrix_trans_pointwise_montgomery(&t_vecv, mat, &t_vech);
    polyvecv_reduce(&t_vecv);
    polyvecv_invntt_tomont(&t_vecv);

    //absorb
    for(i = 0; i < V; ++i){
        PACKQ_ABSORB(t_vecv.vec[i])
    }

    //t_poly = yhat ^ T * w
    polyvech_pointwise_acc_montgomery(&t_poly, &y_TA, &t_vech);
    poly_invntt_tomont(&t_poly);

    //absorb
    PACKQ_ABSORB(t_poly)


    //compute e_{pi+1}
    shake256_finalize(&state);
    shake256_squeeze(bytes, SEEDBYTES, &state);
    poly_challenge(&e, bytes);


    // 3.for i=pi+1,...n,1,...,pi-1, compute e_{i+1}
    i = (PI + 1) % USER_NUM;
    while(i != PI){
        if(i == 0){
            //store e1
            pack_bytes(sig + E1_POS, bytes, SEEDBYTES);
        }
        state = state0;
//

        // 1) \alpha_i = \bar{A} * ri - ei * y_i
        // sample ri = (S_eta^v, S_eta^h) = (t_vecv, t_vech1)
//        TIME_START()
        polyvecv_uniform_eta(&t_vecv,  seed_eta, nonce_eta);
        nonce_eta += V;
        polyvech_uniform_eta(&t_vech1, seed_eta, nonce_eta);
        nonce_eta += H;
//        TIME_PAUSE()
        // pack ri
        pack_respR(sig + RESPONSE_POS + RESPONSE_BYTES * i, &t_vecv, &t_vech1);

        // t_vech = \bar{A} * ri = A * r1 + r2
        polyvecv_ntt(&t_vecv);
        polyvec_matrix_pointwise_montgomery(&t_vech, mat, &t_vecv);
        polyvech_invntt_tomont(&t_vech);
        polyvech_add(&t_vech, &t_vech, &t_vech1);

        // unpack yi in t_vech1
        unpack_y(&t_vech1, inParams->y_list + i * Y_BYTES);
        // t_vech1 = ei * yi
        poly_ntt(&e); //note: ei in NTT form
        polyvech_ntt(&t_vech1);
        polyvech_pointwise_poly_montgomery(&t_vech1, &e, &t_vech1);
        polyvech_invntt_tomont(&t_vech1);

        //t_vech = alpha_i = \bar{A} * ri - ei * yi
        polyvech_sub(&t_vech, &t_vech, &t_vech1);
        polyvech_reduce(&t_vech);
        // pack & absorb alpha_i
        for(j = 0; j < H; ++j){
            PACKQ_ABSORB(t_vech.vec[j])
        }


        // 2) omega_i = A^T*t_i + epsilon * pi - ei * C1
        // sample t_vech = ti <- S_eta^h
//        TIME_START()
        polyvech_uniform_eta(&t_vech, seed_eta, nonce_eta);
//        TIME_PAUSE()
        nonce_eta += H;
        // store ti
        pack_respT(sig + RESPONSE_POS + i * RESPONSE_BYTES + T_POS, &t_vech);

        // t_vecv = A^T * t_i
        polyvech_ntt(&t_vech);// ti in NTT form
        polyvec_matrix_trans_pointwise_montgomery(&t_vecv, mat, &t_vech);
        polyvecv_invntt_tomont(&t_vecv);

        //sample t_poly = vi <- S_beta
//        TIME_START()
        poly_uniform_beta(&t_poly, seed_beta, nonce_beta++);
//        TIME_PAUSE()
        // t_poly = rho_i = ei * vi
        poly_ntt(&t_poly);
        poly_pointwise_montgomery(&t_poly, &t_poly, &e);
        poly_invntt_tomont(&t_poly);
//        poly_reduce(&t_poly);
        pack_respRho(sig + RESPONSE_POS + i * RESPONSE_BYTES + RHO_POS, &t_poly);

        // t_vecv1 = epsilon(v*1) * rho_i,  v * 1v
        t_poly1 = t_poly;
        poly_ntt(&t_poly1);
        polyvecv_pointwise_poly_montgomery(&t_vecv1, &t_poly1, &epsilon);
        polyvecv_invntt_tomont(&t_vecv1);

        // t_vecv = A^T*t_i + epsilon * pi
        polyvecv_add(&t_vecv, &t_vecv, &t_vecv1);

        // t_vecv1 = ei * C_1
        polyvecv_pointwise_poly_montgomery(&t_vecv1, &e, &c1);
        polyvecv_invntt_tomont(&t_vecv1);

        polyvecv_sub(&t_vecv, &t_vecv, &t_vecv1);
        polyvecv_reduce(&t_vecv);

        // pack & absorb
        for(j = 0; j < V; ++j){
            PACKQ_ABSORB(t_vecv.vec[j]);
        }

        // 3) delta_i = y_TA^T * t_i + rho_i - ei * (C_2 - [q/2] * id_i)
        // t_poly1 = y_TA^T * t_i, ti = t_vech
        polyvech_pointwise_acc_montgomery(&t_poly1, &y_TA, &t_vech);
        poly_invntt_tomont(&t_poly1);

        // t_poly = y_TA^T * t_i + rho_i, rho_i = t_poly
        poly_add(&t_poly, &t_poly, &t_poly1);

        // t_poly1 = C_2 - [q/2] * id_i
        poly_frommsg(&t_poly1, inParams->id_list + i * SEEDBYTES);
        poly_sub(&t_poly1, &c2, &t_poly1);

        // t_poly1 = ei * (C_2 - [q/2] * id_i)
        poly_ntt(&t_poly1);
        poly_pointwise_montgomery(&t_poly1, &e, &t_poly1);
        poly_invntt_tomont(&t_poly1);

        // t_poly = delta_i
        poly_sub(&t_poly, &t_poly, &t_poly1);
        poly_reduce(&t_poly);

        // pack & absorb
        PACKQ_ABSORB(t_poly)

        // e_{i+1}
        shake256_finalize(&state);
        shake256_squeeze(bytes, SEEDBYTES, &state);

        poly_challenge(&e, bytes);

        i = (i + 1) % USER_NUM;

    }
    // 4. compute r_pi, t_pi, rho_pi, now e1 is e_{pi}
    // (t_vecv, t_vech) = (e_pi * x_pi, e_pi * x_pi')
    poly_ntt(&e);
    polyvecv_pointwise_poly_montgomery(&t_vecv, &e, &x);
    polyvecv_invntt_tomont(&t_vecv);
    polyvech_pointwise_poly_montgomery(&t_vech, &e, &xprime);
    polyvech_invntt_tomont(&t_vech);

    // (t_vecv, t_vech) = r_pi = (u1 + e_pi * x_pi, u2 + e_pi * x_pi')
    polyvecv_add(&t_vecv, &t_vecv, &u1);
    polyvecv_reduce(&t_vecv);
    polyvech_add(&t_vech, &t_vech, &u2);
    polyvech_reduce(&t_vech);

    // reject sampling
    if(polyvecv_chknorm(&t_vecv,ETA) || polyvech_chknorm(&t_vech, ETA)) {
//        rep++;
        goto rej;
    }
    // t_vech1 = e_pi * s
    polyvech_pointwise_poly_montgomery(&t_vech1, &e, &s);
    polyvech_invntt_tomont(&t_vech1);

    // t_vech1 = t_pi = w + e_pi * s
    polyvech_add(&t_vech1, &w, &t_vech1);
    polyvech_reduce(&t_vech1);

    if(polyvech_chknorm(&t_vech1, ETA)){
//        rep++;
        goto rej;
    }


    // t_poly = rho_pi = e_pi * epsilon_2
    poly_pointwise_montgomery(&t_poly, &e, &epsilon2);
    poly_invntt_tomont(&t_poly);
    poly_reduce(&t_poly);

    pack_respR(sig + RESPONSE_POS + PI * RESPONSE_BYTES, &t_vecv, &t_vech);
    pack_respT(sig + RESPONSE_POS + PI * RESPONSE_BYTES + T_POS, &t_vech1);
    pack_respRho(sig + RESPONSE_POS + PI * RESPONSE_BYTES + RHO_POS, &t_poly);
    if(PI == 0){
        pack_bytes(sig + E1_POS, bytes, SEEDBYTES);
    }
    *siglen = CRYPTO_BYTES;
    return 0;

#undef PACKQ_ABSORB
}

/*************************************************
* Name:        crypto_sign
*
* Description: Compute signed message.
*
* Arguments:   - uint8_t *sm: pointer to output signed message (allocated
*                             array with CRYPTO_BYTES + mlen bytes),
*                             can be equal to m
*              - size_t *smlen: pointer to output length of signed
*                               message
*              - const uint8_t *m: pointer to message to be signed
*              - size_t mlen: length of message
*              - const uint8_t *sk: pointer to bit-packed secret key
*              - const uint8_t *pk_list: pk of users in the ring
*              - const uint8_t pi: index of signer in the ring
*              - uint8_t seedbuf[SEEDBYTES]:
*
* Returns 0 (success)
**************************************************/
int crypto_sign(uint8_t *sm,
                size_t *smlen,
                const uint8_t *m,
                size_t mlen,
                const SIGN_IN_PARAMS *inParams) {
    size_t i;

    for (i = 0; i < mlen; ++i)
        sm[CRYPTO_BYTES + mlen - 1 - i] = m[mlen - 1 - i];
    crypto_sign_signature(sm, smlen, sm + CRYPTO_BYTES, mlen, inParams);
    *smlen += mlen;
    return 0;
}

/*************************************************
* Name:        crypto_sign_verify
*
* Description: Verifies signature.
*
* Arguments:   - uint8_t *m: pointer to input signature
*              - size_t siglen: length of signature
*              - const uint8_t *m: pointer to message
*              - size_t mlen: length of message
*              - const uint8_t *pk: pointer to bit-packed public key
*
* Returns 0 if signature could be verified correctly and -1 otherwise
**************************************************/

int crypto_sign_verify(const uint8_t *sig,
                       size_t siglen,
                       const uint8_t *m,
                       size_t mlen,
                       const VERIFY_IN_PARAMS *inParams) {
    unsigned int i, j;
    uint8_t bytes[POLYQ_PACKEDBYTES];
    polyvecv c1, epsilon, t_vecv, t_vecv1;
    polyvech y_TA, t_vech, t_vech1;
    poly c2, t_poly, t_poly1, e;
    keccak_state state, state0;

#define PACKQ_ABSORB(p)    polyq_pack(bytes, &p);\
            shake256_absorb(&state, bytes, POLYQ_PACKEDBYTES);
    // restore
    poly_challenge(&e, sig + E1_POS);
    unpack_ctag(&c1, &c2, sig +  CTAG_POS);
    unpack_epsilon(&epsilon, sig);
    unpack_y(&y_TA, inParams->y_TA);

    //ntt prepare
    polyvecv_ntt(&c1);//note: c1 in NTT form
    polyvecv_ntt(&epsilon);//note: epsilon in NTT form
    polyvech_ntt(&y_TA);//note: y_TA in NTT form


    // keccak absorb T L M
    shake256_init(&state0);
    shake256_absorb(&state0, sig + TIMESTAMP_POS, TIMESTAMP_BYTES);
    shake256_absorb(&state0, inParams->y_list, USER_NUM * Y_BYTES);
    shake256_absorb(&state0, m, mlen);

    // compute e_{i+1}, i = 1,..., n - 1
    for(i = 0; i < USER_NUM; ++i){
        state = state0;

        // 1) \alpha_i = \bar{A} * ri - ei * y_i
        // unpack ri = (t_vecv, t_vech1)
        unpack_respR(&t_vecv, &t_vech1, sig + RESPONSE_POS + i * RESPONSE_BYTES);

        // t_vech = \bar{A} * ri = A * r1 + r2
        polyvecv_ntt(&t_vecv);
        polyvec_matrix_pointwise_montgomery(&t_vech, mat, &t_vecv);
        polyvech_invntt_tomont(&t_vech);
        polyvech_add(&t_vech, &t_vech, &t_vech1);

        // unpack yi in t_vech1
        unpack_y(&t_vech1, inParams->y_list + i * Y_BYTES);
        // t_vech1 = ei * yi
        poly_ntt(&e); //note: ei in NTT form
        polyvech_ntt(&t_vech1);
        polyvech_pointwise_poly_montgomery(&t_vech1, &e, &t_vech1);
        polyvech_invntt_tomont(&t_vech1);

        //t_vech = alpha_i = \bar{A} * ri - ei * yi
        polyvech_sub(&t_vech, &t_vech, &t_vech1);
        polyvech_reduce(&t_vech);
        // pack & absorb alpha_i
        for(j = 0; j < H; ++j){
            PACKQ_ABSORB(t_vech.vec[j])
        }


        // 2) omega_i = A^T*t_i + epsilon * pi - ei * C1
        // unpack t_vech = ti
        unpack_respT(&t_vech, sig + RESPONSE_POS + i * RESPONSE_BYTES + T_POS);

        // t_vecv = A^T * t_i
        polyvech_ntt(&t_vech);// ti in NTT form
        polyvec_matrix_trans_pointwise_montgomery(&t_vecv, mat, &t_vech);
        polyvecv_invntt_tomont(&t_vecv);

        // unpack t_poly = rho_i
        unpack_respRho(&t_poly, sig + RESPONSE_POS + i * RESPONSE_BYTES + RHO_POS);

        // t_vecv1 = epsilon(v*1) * rho_i,  v * 1v
        t_poly1 = t_poly;
        poly_ntt(&t_poly1);
        polyvecv_pointwise_poly_montgomery(&t_vecv1, &t_poly1, &epsilon);
        polyvecv_invntt_tomont(&t_vecv1);

        // t_vecv = A^T*t_i + epsilon * pi
        polyvecv_add(&t_vecv, &t_vecv, &t_vecv1);

        // t_vecv1 = ei * C_1
        polyvecv_pointwise_poly_montgomery(&t_vecv1, &e, &c1);
        polyvecv_invntt_tomont(&t_vecv1);

        polyvecv_sub(&t_vecv, &t_vecv, &t_vecv1);
        polyvecv_reduce(&t_vecv);

        // pack & absorb
        for(j = 0; j < V; ++j){
            PACKQ_ABSORB(t_vecv.vec[j]);
        }

        // 3) delta_i = y_TA^T * t_i + rho_i - ei * (C_2 - [q/2] * id_i)
        // t_poly1 = y_TA^T * t_i, ti = t_vech
        polyvech_pointwise_acc_montgomery(&t_poly1, &y_TA, &t_vech);
        poly_invntt_tomont(&t_poly1);

        // t_poly = y_TA^T * t_i + rho_i, rho_i = t_poly
        poly_add(&t_poly, &t_poly, &t_poly1);

        // t_poly1 = C_2 - [q/2] * id_i
        poly_frommsg(&t_poly1, inParams->id_list + i * SEEDBYTES);
        poly_sub(&t_poly1, &c2, &t_poly1);

        // t_poly1 = ei * (C_2 - [q/2] * id_i)
        poly_ntt(&t_poly1);
        poly_pointwise_montgomery(&t_poly1, &e, &t_poly1);
        poly_invntt_tomont(&t_poly1);

        // t_poly = delta_i
        poly_sub(&t_poly, &t_poly, &t_poly1);
        poly_reduce(&t_poly);

        // pack & absorb
        PACKQ_ABSORB(t_poly)

        // e_{i+1}
        shake256_finalize(&state);
        shake256_squeeze(bytes, SEEDBYTES, &state);

        poly_challenge(&e, bytes);

    }

    for(j = 0; j < SEEDBYTES; ++j){
        if(bytes[j] != sig[E1_POS + j]){
            return 1;

        }
    }

    return 0;


#undef PACKQ_ABSORB
}

/*************************************************
* Name:        crypto_sign_open
*
* Description: Verify signed message.
*
* Arguments:   - uint8_t *m: pointer to output message (allocated
*                            array with smlen bytes), can be equal to sm
*              - size_t *mlen: pointer to output length of message
*              - const uint8_t *sm: pointer to signed message
*              - size_t smlen: length of signed message
*              - const uint8_t *pk: pointer to bit-packed public key
*
* Returns 0 if signed message could be verified correctly and -1 otherwise
**************************************************/
int crypto_sign_open(size_t mlen,
                     const uint8_t *sm,
                     size_t smlen,
                     const VERIFY_IN_PARAMS *inParams) {
    size_t i;

    if (smlen < CRYPTO_BYTES)
        goto badsig;

    if (crypto_sign_verify(sm, CRYPTO_BYTES, sm + CRYPTO_BYTES, mlen, inParams))
        goto badsig;
    else {
        /* All good, return 0 */
        return 0;
    }

    badsig:
    /* Signature verification failed */
    return -1;

}
/*************************************************
* Name:        crypto_sign_open
*
* Description: TA revoke signer's identity.
*
* Arguments:   -
*
* Returns signer's index
**************************************************/
int crypto_revoke(REVOKE_IN_PARMAS *inParmas){
    unsigned int i;
    polyvecv x_TA, c1;
    polyvech xprime_TA;
    poly c2, id_pi, id_i;


    //restore
    unpack_ctag(&c1, &c2, inParmas->ctag);
    unpack_sk(&x_TA, &xprime_TA, inParmas->sk_TA);

    polyvecv_ntt(&x_TA);
    polyvecv_ntt(&c1);
    polyvecv_pointwise_acc_montgomery(&id_pi, &x_TA, &c1);
    poly_invntt_tomont(&id_pi);
    poly_sub(&id_pi, &c2, &id_pi);
    poly_reduce(&id_pi);
    poly_round(&id_pi);

    for(i = 0; i < USER_NUM; ++i){
        poly_frommsg(&id_i, inParmas->id_list + i * SEEDBYTES);
        if(poly_cmp(&id_pi, &id_i)){
            return i;
        }
    }


    return -1;

}
