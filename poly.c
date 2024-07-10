#include <stdint.h>
#include "params.h"
#include "poly.h"
#include "ntt.h"
#include "reduce.h"

#include "symmetric.h"

#ifdef DBENCH
#include "test/cpucycles.h"
extern const uint64_t timing_overhead;
extern uint64_t *tred, *tadd, *tmul, *tround, *tsample, *tpack;
#define DBENCH_START() uint64_t time = cpucycles()
#define DBENCH_STOP(t) t += cpucycles() - time - timing_overhead
#else
#define DBENCH_START()
#define DBENCH_STOP(t)
#endif

/*************************************************
* Name:        poly_reduce
*
* Description: Inplace reduction of all coefficients of polynomial to
*              representative in [-6283009,6283007].
*
* Arguments:   - poly *a: pointer to input/output polynomial
**************************************************/
void poly_reduce(poly *a) {
  unsigned int i;
  DBENCH_START();

  for(i = 0; i < N; ++i)
    a->coeffs[i] = reduce32(a->coeffs[i]);

  DBENCH_STOP(*tred);
}

/*************************************************
* Name:        poly_caddq
*
* Description: For all coefficients of in/out polynomial add Q if
*              coefficient is negative.
*
* Arguments:   - poly *a: pointer to input/output polynomial
**************************************************/
void poly_caddq(poly *a) {
  unsigned int i;
  DBENCH_START();

  for(i = 0; i < N; ++i)
    a->coeffs[i] = caddq(a->coeffs[i]);

  DBENCH_STOP(*tred);
}

/*************************************************
* Name:        poly_freeze
*
* Description: Inplace reduction of all coefficients of polynomial to
*              standard representatives.
*
* Arguments:   - poly *a: pointer to input/output polynomial
**************************************************/
void poly_freeze(poly *a) {
  unsigned int i;
  DBENCH_START();

  for(i = 0; i < N; ++i)
    a->coeffs[i] = freeze(a->coeffs[i]);

  DBENCH_STOP(*tred);
}

/*************************************************
* Name:        poly_add
*
* Description: Add polynomials. No modular reduction is performed.
*
* Arguments:   - poly *c: pointer to output polynomial
*              - const poly *a: pointer to first summand
*              - const poly *b: pointer to second summand
**************************************************/
void poly_add(poly *c, const poly *a, const poly *b)  {
  unsigned int i;
  DBENCH_START();

  for(i = 0; i < N; ++i)
    c->coeffs[i] = a->coeffs[i] + b->coeffs[i];

  DBENCH_STOP(*tadd);
}

/*************************************************
* Name:        poly_sub
*
* Description: Subtract polynomials. No modular reduction is
*              performed.
*
* Arguments:   - poly *c: pointer to output polynomial
*              - const poly *a: pointer to first input polynomial
*              - const poly *b: pointer to second input polynomial to be
*                               subtraced from first input polynomial
**************************************************/
void poly_sub(poly *c, const poly *a, const poly *b) {
  unsigned int i;
  DBENCH_START();

  for(i = 0; i < N; ++i)
    c->coeffs[i] = a->coeffs[i] - b->coeffs[i];

  DBENCH_STOP(*tadd);
}

/*************************************************
* Name:        poly_cmp
*
* Description: Compare two polynomials.
*
* Arguments:   - poly *c: pointer to output polynomial
*              - const poly *a: pointer to first input polynomial
*              - const poly *b: pointer to second input polynomial to be
*                               compared with first input polynomial
*
* Return: 1:equal/0:else
**************************************************/
int poly_cmp(poly *a, poly *b){
    unsigned int i;
    for(i = 0; i < N; ++i) {
        if(a->coeffs[i] != b->coeffs[i]) {
            return 0;
        }
    }
    return 1;

}


/*************************************************
* Name:        poly_ntt
*
* Description: Inplace forward NTT. Coefficients can grow by
*              8*Q in absolute value.
*
* Arguments:   - poly *a: pointer to input/output polynomial
**************************************************/
void poly_ntt(poly *a) {
  DBENCH_START();

  ntt(a->coeffs);

  DBENCH_STOP(*tmul);
}

/*************************************************
* Name:        poly_invntt_tomont
*
* Description: Inplace inverse NTT and multiplication by 2^{32}.
*              Input coefficients need to be less than Q in absolute
*              value and output coefficients are again bounded by Q.
*
* Arguments:   - poly *a: pointer to input/output polynomial
**************************************************/
void poly_invntt_tomont(poly *a) {
  DBENCH_START();

  invntt_tomont(a->coeffs, 41978);

  DBENCH_STOP(*tmul);
}

/*************************************************
* Name:        poly_pointwise_montgomery
*
* Description: Pointwise multiplication of polynomials in NTT domain
*              representation and multiplication of resulting polynomial
*              by 2^{-32}.
*
* Arguments:   - poly *c: pointer to output polynomial
*              - const poly *a: pointer to first input polynomial
*              - const poly *b: pointer to second input polynomial
**************************************************/
void poly_pointwise_montgomery(poly *c, const poly *a, const poly *b) {
  unsigned int i;
  DBENCH_START();

  for(i = 0; i < N; ++i)
    c->coeffs[i] = montgomery_reduce((int64_t)a->coeffs[i] * b->coeffs[i]);

  DBENCH_STOP(*tmul);
}

#define mont_mul(x, y) montgomery_reduce((int64_t)(x)*(y))
#define mont_sqr(x) montgomery_reduce((int64_t)(x)*(x))
/*************************************************
* Name:        modinv
*
* Description: compute a^-1 mod Q = a^{QD-2} mod QD
*              required: a > 0
*
* Arguments:   - int32_t: finite field element a
*
* Returns r in montgomery form
**************************************************/
static int32_t modinv(int32_t a){
    //Q-2 = 0x7fdfff
    int32_t t0, t1, t2;
    a = mont_mul(a, R2); //mont epr:a*R^2
    t0 = a;
    t1 = mont_sqr(t0);
    t2 = mont_mul(t0, t1);//a^3
    t0 = mont_sqr(t2);
    t0 = mont_sqr(t0);//a^{0xc}
    t1 = mont_mul(t0, t2);//a^{0xf}
    t0 = mont_sqr(t1);
    t0 = mont_sqr(t0);
    t0 = mont_sqr(t0);
    t0 = mont_sqr(t0);//a^{0xf0}
    t2 = mont_mul(t0, t1);//a^{0xff}
    t0 = mont_sqr(t2);//a^{0x1fe}
    t0 = mont_mul(t0, a);//a^{0x1ff}
    t2 = mont_sqr(t0);
    t2 = mont_sqr(t2);
    t2 = mont_sqr(t2);
    t2 = mont_sqr(t2);
    t2 = mont_sqr(t2);
    t2 = mont_sqr(t2);
    t2 = mont_sqr(t2);
    t2 = mont_sqr(t2);
    t2 = mont_sqr(t2);
    t2 = mont_sqr(t2);//a^{0x7fc000}
    t2 = mont_mul(t2, t0);//a^{0x7fdff}
    t0 = mont_sqr(t2);
    t0 = mont_sqr(t0);
    t0 = mont_sqr(t0);
    t0 = mont_sqr(t0);//a^{0x7fdff0}
    t2 = mont_mul(t0, t1);//a^{0x7fdfff]

    return t2;
}

/*************************************************
* Name:        poly_inv
*
* Description: compute r = a^-1 mod qd where a is in NTT form
*
* Arguments:
*
* Return: 0 not invertible/1
**************************************************/
int poly_inv(poly* a){
    poly_ntt(a);
    DBENCH_START();
    for(unsigned int i = 0; i < N; ++i){
        if(a->coeffs[i] == 0)
            // not invertible
            return 0;
//        a[i] = reduce32(a[i]);
        a->coeffs[i] = modinv(a->coeffs[i]);
    }
    poly_reduce(a);
    invntt_tomont(a->coeffs, 8347681);
    DBENCH_STOP(*tsample);
    return 1;

}

///*************************************************
//* Name:        poly_gen_epsilon2
//*
//* Description: Sample invertible polynomial epsilon2 with uniformly random coefficients
//*              in [-BETA, BETA] by performing rejection sampling on the
//*              output stream of SHAKE256(seed|nonce), and compute its inverse i_epsilon2
//*
//* Arguments:   - poly *epsilon2: pointer to output invertible polynomial
//*              - poly *i_epsilon2: pointer to output inverse of epsilon2
//*              - const uint8_t seed[]: byte array with seed of length SEEDBYTES
//*              - uint16_t nonce: 2-byte nonce
//* Return:  retry times
//**************************************************/
//uint16_t poly_gen_epsilon2(poly *epsilon2,
//                           poly *i_epsilon2,
//                           const uint8_t seed[SEEDBYTES],
//                           uint16_t nonce){
//    uint16_t T = 0;
//
//retry:
//    poly_uniform_beta(epsilon2, seed, nonce + T);
//
//
//
//}

/*************************************************
* Name:        poly_chknorm
*
* Description: Check infinity norm of polynomial against given bound.
*              Assumes input coefficients were reduced by reduce32().
*
* Arguments:   - const poly *a: pointer to polynomial
*              - int32_t B: norm bound
*
* Returns 0 if norm is strictly smaller than B <= (Q-1)/8 and 1 otherwise.
**************************************************/
int poly_chknorm(const poly *a, int32_t B) {
  unsigned int i;
  int32_t t;
  DBENCH_START();

  if(B > (Q-1)/8)
    return 1;

  /* It is ok to leak which coefficient violates the bound since
     the probability for each coefficient is independent of secret
     data but we must not leak the sign of the centralized representative. */
  for(i = 0; i < N; ++i) {
    /* Absolute value */
    t = a->coeffs[i] >> 31;
    t = a->coeffs[i] - (t & 2*a->coeffs[i]);

    if(t >= B) {
      DBENCH_STOP(*tsample);
      return 1;
    }
  }

  DBENCH_STOP(*tsample);
  return 0;
}

/*************************************************
* Name:        poly_frommsg
*
* Description: Convert 32-byte message to polynomial with coefficients in {0, q/2}
*
* Arguments:   - const poly *a: pointer to polynomial
*              - int32_t msg: message
*
* Returns
**************************************************/
void poly_frommsg(poly *a, const uint8_t msg[SEEDBYTES]) {
    unsigned int i,j;
    int16_t mask;
    for(i = 0; i < SEEDBYTES; i++) {
        for(j = 0; j <  8; j++) {
            mask = -(int16_t)((msg[i] >> j)&1);
            a->coeffs[8*i+j] = mask & (Q_HALF);
        }
    }
}
/*************************************************
* Name:        poly_round
*
* Description: In-place round the coefficients of p, where c <= |q/4| to 0, else to 1
*
* Arguments:   - const poly *a: pointer to polynomial
*              - int32_t msg: message
*
* Returns
**************************************************/
void poly_round(poly *a) {
    unsigned int i;
    int32_t t;
    for(i = 0; i < N; ++i) {
        t  = a->coeffs[i];
        t += (t >> 31) & Q;
        t = (((t << 1) + Q_HALF) / Q) & 1;
        a->coeffs[i] = (-t >> 31) & Q_HALF;
    }
}

/*************************************************
* Name:        rej_uniform
*
* Description: Sample uniformly random coefficients in [0, Q-1] by
*              performing rejection sampling on array of random bytes.
*
* Arguments:   - int32_t *a: pointer to output array (allocated)
*              - unsigned int len: number of coefficients to be sampled
*              - const uint8_t *buf: array of random bytes
*              - unsigned int buflen: length of array of random bytes
*
* Returns number of sampled coefficients. Can be smaller than len if not enough
* random bytes were given.
**************************************************/
static unsigned int rej_uniform(int32_t *a,
                                unsigned int len,
                                const uint8_t *buf,
                                unsigned int buflen)
{
  unsigned int ctr, pos;
  uint32_t t;
  DBENCH_START();

  ctr = pos = 0;
  while(ctr < len && pos + 3 <= buflen) {
    t  = buf[pos++];
    t |= (uint32_t)buf[pos++] << 8;
    t |= (uint32_t)buf[pos++] << 16;
    t &= 0x7FFFFF;
    if(t < Q)
      a[ctr++] = t;
  }

  DBENCH_STOP(*tsample);
  return ctr;
}

/*************************************************
* Name:        poly_uniform
*
* Description: Sample polynomial with uniformly random coefficients
*              in [0,Q-1] by performing rejection sampling on the
*              output stream of SHAKE256(seed|nonce) or AES256CTR(seed,nonce).
*
* Arguments:   - poly *a: pointer to output polynomial
*              - const uint8_t seed[]: byte array with seed of length SEEDBYTES
*              - uint16_t nonce: 2-byte nonce
**************************************************/
#define POLY_UNIFORM_NBLOCKS ((768 + STREAM128_BLOCKBYTES - 1)/STREAM128_BLOCKBYTES)
void poly_uniform(poly *a,
                  const uint8_t seed[SEEDBYTES],
                  uint16_t nonce)
{
  unsigned int i, ctr, off;
  unsigned int buflen = POLY_UNIFORM_NBLOCKS*STREAM128_BLOCKBYTES;
  uint8_t buf[POLY_UNIFORM_NBLOCKS*STREAM128_BLOCKBYTES + 2];
  stream128_state state;

  stream128_init(&state, seed, nonce);
  stream128_squeezeblocks(buf, POLY_UNIFORM_NBLOCKS, &state);

  ctr = rej_uniform(a->coeffs, N, buf, buflen);

  while(ctr < N) {
    off = buflen % 3;//3是p需要3个字节？
    for(i = 0; i < off; ++i)//交换最后off和前off个字节？
      buf[i] = buf[buflen - off + i];

    stream128_squeezeblocks(buf + off, 1, &state);
    buflen = STREAM128_BLOCKBYTES + off;
    ctr += rej_uniform(a->coeffs + ctr, N - ctr, buf, buflen);
  }
}

/*************************************************
* Name:        poly_uniform_gamma1m1
*
* Description: Sample polynomial with uniformly random coefficients
*              in [-(GAMMA1 - 1), GAMMA1] by unpacking output stream
*              of SHAKE256(seed|nonce) or AES256CTR(seed,nonce).
*
* Arguments:   - poly *a: pointer to output polynomial
*              - const uint8_t seed[]: byte array with seed of length CRHBYTES
*              - uint16_t nonce: 16-bit nonce
**************************************************/
#define POLY_UNIFORM_GAMMA1_NBLOCKS ((576 + STREAM256_BLOCKBYTES - 1)/STREAM256_BLOCKBYTES)
void poly_uniform_gamma(poly *a,
                        const uint8_t seed[CRHBYTES],
                        uint16_t nonce){
    uint8_t buf[POLY_UNIFORM_GAMMA1_NBLOCKS*STREAM256_BLOCKBYTES];
    stream256_state state;

    stream256_init(&state, seed, nonce);
    stream256_squeezeblocks(buf, POLY_UNIFORM_GAMMA1_NBLOCKS, &state);
    polyeta_unpack(a, buf);

}

/*************************************************
* Name:        rej_beta
*
* Description: Sample uniformly random coefficients in [-BETA, BETA] by
*              performing rejection sampling on array of random bytes.
*
* Arguments:   - int32_t *a: pointer to output array (allocated)
*              - unsigned int len: number of coefficients to be sampled
*              - const uint8_t *buf: array of random bytes
*              - unsigned int buflen: length of array of random bytes
*
* Returns number of sampled coefficients. Can be smaller than len if not enough
* random bytes were given.
**************************************************/
static unsigned int rej_beta(int32_t *a,
                            unsigned int len,
                            const uint8_t *buf,
                            unsigned int buflen)
{
  unsigned int ctr, pos;
  uint32_t t0, t1;
  DBENCH_START();

  ctr = pos = 0;
  while(ctr < len && pos < buflen) {
    t0 = buf[pos] & 0x0F;
    t1 = buf[pos++] >> 4;

#if BETA == 2
    if(t0 < 15) {
      t0 = t0 - (205*t0 >> 10)*5;
      a[ctr++] = 2 - t0;
    }
    if(t1 < 15 && ctr < len) {
      t1 = t1 - (205*t1 >> 10)*5;
      a[ctr++] = 2 - t1;
    }
#elif BETA == 3
    if(t0 < 12) {
        t0 = t0 - (171 * t0 >> 10) * 6;
        a[ctr++] = 3 - t0;
    }

    if(t1 < 12) {
        t1 = t1 - (171 * t1 >> 10) * 6;
        a[ctr++] = 3 - t1;
    }
#endif
  }

  DBENCH_STOP(*tsample);
  return ctr;
}

/*************************************************
* Name:        poly_uniform_beta
*
* Description: Sample polynomial with uniformly random coefficients
*              in [-BETA,BETA] by performing rejection sampling on the
*              output stream from SHAKE256(seed|nonce) or AES256CTR(seed,nonce).
*
* Arguments:   - poly *a: pointer to output polynomial
*              - const uint8_t seed[]: byte array with seed of length SEEDBYTES
*              - uint16_t nonce: 2-byte nonce
**************************************************/
// 推测是beta = 2时每个字节有15/16的概率被使用，beta=4时有9/16概率被使用，尽量保证总共使用的字节相同（大概为128字节）
#if BETA == 2
#define POLY_UNIFORM_BETA_NBLOCKS ((136 + STREAM128_BLOCKBYTES - 1)/STREAM128_BLOCKBYTES)
#elif BETA == 3
#define POLY_UNIFORM_BETA_NBLOCKS ((170 + STREAM128_BLOCKBYTES - 1)/STREAM128_BLOCKBYTES)
#endif
void poly_uniform_beta(poly *a,
                      const uint8_t seed[SEEDBYTES],
                      uint16_t nonce)
{
  unsigned int ctr;
  unsigned int buflen = POLY_UNIFORM_BETA_NBLOCKS*STREAM128_BLOCKBYTES;
  uint8_t buf[POLY_UNIFORM_BETA_NBLOCKS*STREAM128_BLOCKBYTES];
  stream128_state state;

  stream128_init(&state, seed, nonce);
  stream128_squeezeblocks(buf, POLY_UNIFORM_BETA_NBLOCKS, &state);

  ctr = rej_beta(a->coeffs, N, buf, buflen);

  while(ctr < N) {
    stream128_squeezeblocks(buf, 1, &state);
    ctr += rej_beta(a->coeffs + ctr, N - ctr, buf, STREAM128_BLOCKBYTES);
  }
}

/*************************************************
* Name:        rej_eta
*
* Description: Sample uniformly random coefficients in [-ETA, ETA] by
*              performing rejection sampling on array of random bytes.
*
* Arguments:   - int32_t *a: pointer to output array (allocated)
*              - unsigned int len: number of coefficients to be sampled
*              - const uint8_t *buf: array of random bytes
*              - unsigned int buflen: length of array of random bytes
*
* Returns number of sampled coefficients. Can be smaller than len if not enough
* random bytes were given.
**************************************************/
#define ETA_D  (2 * ETA - 1)
static unsigned int rej_eta(int32_t *a,
                            unsigned int len,
                            const uint8_t *buf,
                            unsigned int buflen){
    unsigned int ctr, pos;
    uint32_t t;
    ctr = pos = 0;
    DBENCH_START();
    while(ctr < len && pos + 3 <= buflen){
        t = (buf[pos] | (buf[pos + 1] << 8) | (buf[pos + 2] << 16)) & 0x3FFFF;
        pos += 2;
        if(t < ETA_D){
            a[ctr++] = ETA - t;
        }

        if(ctr < len && pos + 3 <= buflen){
            t = ((buf[pos] >> 2) | (buf[pos + 1] << 6) | (buf[pos + 2] << 14)) & 0x3FFFF;
            if(t < ETA_D){
                a[ctr++] = ETA - t;
            }
        }
        pos += 2;

        if(ctr < len && pos + 3 <= buflen){
            t = ((buf[pos] >> 4) | (buf[pos + 1] << 4) | (buf[pos + 2] << 12)) & 0x3FFFF;
            if(t < ETA_D){
                a[ctr++] = ETA - t;
            }
        }
        pos += 2;

        if(ctr < len && pos + 3 <= buflen){
            t = ((buf[pos] >> 6) | (buf[pos + 1] << 2) | (buf[pos + 2] << 10)) & 0x3FFFF;
            if(t < ETA_D){
                a[ctr++] = ETA - t;
            }
        }
        pos += 3;

    }
    DBENCH_STOP(*tsample);
    return ctr;

}

/*************************************************
* Name:        poly_uniform_eta
*
* Description: Sample polynomial with uniformly random coefficients
*              in [-ETA,ETA] by performing rejection sampling on the
*              output stream from SHAKE256(seed|nonce)
*
* Arguments:   - poly *a: pointer to output polynomial
*              - const uint8_t seed[]: byte array with seed of length SEEDBYTES
*              - uint16_t nonce: 2-byte nonce
**************************************************/
#define POLY_UNIFORM_ETA_NBLOCKS ((POLYETA_PACKEDBYTES + STREAM128_BLOCKBYTES - 1)/STREAM128_BLOCKBYTES)
void poly_uniform_eta(poly *a,
                       const uint8_t seed[SEEDBYTES],
                       uint16_t nonce) {
    unsigned int ctr;
    unsigned int buflen = POLY_UNIFORM_ETA_NBLOCKS*STREAM128_BLOCKBYTES;
    uint8_t buf[POLY_UNIFORM_ETA_NBLOCKS*STREAM128_BLOCKBYTES];
    stream128_state state;

    stream128_init(&state, seed, nonce);
    stream128_squeezeblocks(buf, POLY_UNIFORM_ETA_NBLOCKS, &state);

    ctr = rej_eta(a->coeffs, N, buf, buflen);

    while(ctr < N) {
        stream128_squeezeblocks(buf, 1, &state);
        ctr += rej_eta(a->coeffs + ctr, N - ctr, buf, STREAM128_BLOCKBYTES);
    }
}


/*************************************************
* Name:        challenge
*
* Description: Implementation of H. Samples polynomial with KAPPA nonzero
*              coefficients in {-1,1} using the output stream of
*              SHAKE256(seed).
*
* Arguments:   - poly *c: pointer to output polynomial
*              - const uint8_t mu[]: byte array containing seed of length SEEDBYTES
**************************************************/
void poly_challenge(poly *c, const uint8_t seed[SEEDBYTES]) {
  unsigned int i, b, pos;
  uint64_t signs;
  uint8_t buf[SHAKE256_RATE];
  keccak_state state;

  shake256_init(&state);
  shake256_absorb(&state, seed, SEEDBYTES);
  shake256_finalize(&state);
  shake256_squeezeblocks(buf, 1, &state);

  signs = 0;
  for(i = 0; i < 8; ++i)
    signs |= (uint64_t)buf[i] << 8*i;
  pos = 8;

  for(i = 0; i < N; ++i)
    c->coeffs[i] = 0;
  for(i = N-KAPPA; i < N; ++i) {
    do {
      if(pos >= SHAKE256_RATE) {
        shake256_squeezeblocks(buf, 1, &state);
        pos = 0;
      }

      b = buf[pos++];
    } while(b > i); //需要b小于i，交换b和i的值

    c->coeffs[i] = c->coeffs[b];
    c->coeffs[b] = 1 - 2*(signs & 1);
    signs >>= 1;
  }
}

/*************************************************
* Name:        polybeta_pack
*
* Description: Bit-pack polynomial with coefficients in [-BETA,BETA].
*
* Arguments:   - uint8_t *r: pointer to output byte array with at least
*                            POLYBETA_PACKEDBYTES bytes
*              - const poly *a: pointer to input polynomial
**************************************************/
void polybeta_pack(uint8_t *r, const poly *a) {
  unsigned int i;
  uint8_t t[8];
  DBENCH_START();


  for(i = 0; i < N/8; ++i) {
    t[0] = BETA - a->coeffs[8*i+0];
    t[1] = BETA - a->coeffs[8*i+1];
    t[2] = BETA - a->coeffs[8*i+2];
    t[3] = BETA - a->coeffs[8*i+3];
    t[4] = BETA - a->coeffs[8*i+4];
    t[5] = BETA - a->coeffs[8*i+5];
    t[6] = BETA - a->coeffs[8*i+6];
    t[7] = BETA - a->coeffs[8*i+7];

    r[3*i+0]  = (t[0] >> 0) | (t[1] << 3) | (t[2] << 6);
    r[3*i+1]  = (t[2] >> 2) | (t[3] << 1) | (t[4] << 4) | (t[5] << 7);
    r[3*i+2]  = (t[5] >> 1) | (t[6] << 2) | (t[7] << 5);
  }



  DBENCH_STOP(*tpack);
}

/*************************************************
* Name:        polybeta_unpack
*
* Description: Unpack polynomial with coefficients in [-BETA,BETA].
*
* Arguments:   - poly *r: pointer to output polynomial
*              - const uint8_t *a: byte array with bit-packed polynomial
**************************************************/
void polybeta_unpack(poly *r, const uint8_t *a) {
  unsigned int i;
  DBENCH_START();


  for(i = 0; i < N/8; ++i) {
    r->coeffs[8*i+0] =  (a[3*i+0] >> 0) & 7;
    r->coeffs[8*i+1] =  (a[3*i+0] >> 3) & 7;
    r->coeffs[8*i+2] = ((a[3*i+0] >> 6) | (a[3*i+1] << 2)) & 7;
    r->coeffs[8*i+3] =  (a[3*i+1] >> 1) & 7;
    r->coeffs[8*i+4] =  (a[3*i+1] >> 4) & 7;
    r->coeffs[8*i+5] = ((a[3*i+1] >> 7) | (a[3*i+2] << 1)) & 7;
    r->coeffs[8*i+6] =  (a[3*i+2] >> 2) & 7;
    r->coeffs[8*i+7] =  (a[3*i+2] >> 5) & 7;

    r->coeffs[8*i+0] = BETA - r->coeffs[8*i+0];
    r->coeffs[8*i+1] = BETA - r->coeffs[8*i+1];
    r->coeffs[8*i+2] = BETA - r->coeffs[8*i+2];
    r->coeffs[8*i+3] = BETA - r->coeffs[8*i+3];
    r->coeffs[8*i+4] = BETA - r->coeffs[8*i+4];
    r->coeffs[8*i+5] = BETA - r->coeffs[8*i+5];
    r->coeffs[8*i+6] = BETA - r->coeffs[8*i+6];
    r->coeffs[8*i+7] = BETA - r->coeffs[8*i+7];
  }

  DBENCH_STOP(*tpack);
}


/*************************************************
* Name:        polyeta_pack
*
* Description: Bit-pack polynomial with coefficients
*              in [-ETA,ETA].
*
* Arguments:   - uint8_t *r: pointer to output byte array with at least
*                            POLYETA_PACKEDBYTES bytes
*              - const poly *a: pointer to input polynomial
**************************************************/
void polyeta_pack(uint8_t *r, const poly *a) {
  unsigned int i;
  uint32_t t[4];
  DBENCH_START();


  for(i = 0; i < N/4; ++i) {
    t[0] = GAMMA - a->coeffs[4*i+0];
    t[1] = GAMMA - a->coeffs[4*i+1];
    t[2] = GAMMA - a->coeffs[4*i+2];
    t[3] = GAMMA - a->coeffs[4*i+3];

    r[9*i+0]  = t[0];
    r[9*i+1]  = t[0] >> 8;
    r[9*i+2]  = t[0] >> 16;
    r[9*i+2] |= t[1] << 2;
    r[9*i+3]  = t[1] >> 6;
    r[9*i+4]  = t[1] >> 14;
    r[9*i+4] |= t[2] << 4;
    r[9*i+5]  = t[2] >> 4;
    r[9*i+6]  = t[2] >> 12;
    r[9*i+6] |= t[3] << 6;
    r[9*i+7]  = t[3] >> 2;
    r[9*i+8]  = t[3] >> 10;
  }


  DBENCH_STOP(*tpack);
}

/*************************************************
* Name:        polyeta_unpack
*
* Description: Unpack polynomial z with coefficients
*              in [-ETA,ETA].
*
* Arguments:   - poly *r: pointer to output polynomial
*              - const uint8_t *a: byte array with bit-packed polynomial
**************************************************/
void polyeta_unpack(poly *r, const uint8_t *a) {
  unsigned int i;
  DBENCH_START();


  for(i = 0; i < N/4; ++i) {
    r->coeffs[4*i+0]  = a[9*i+0];
    r->coeffs[4*i+0] |= (uint32_t)a[9*i+1] << 8;
    r->coeffs[4*i+0] |= (uint32_t)a[9*i+2] << 16;
    r->coeffs[4*i+0] &= 0x3FFFF;

    r->coeffs[4*i+1]  = a[9*i+2] >> 2;
    r->coeffs[4*i+1] |= (uint32_t)a[9*i+3] << 6;
    r->coeffs[4*i+1] |= (uint32_t)a[9*i+4] << 14;
    r->coeffs[4*i+1] &= 0x3FFFF;

    r->coeffs[4*i+2]  = a[9*i+4] >> 4;
    r->coeffs[4*i+2] |= (uint32_t)a[9*i+5] << 4;
    r->coeffs[4*i+2] |= (uint32_t)a[9*i+6] << 12;
    r->coeffs[4*i+2] &= 0x3FFFF;

    r->coeffs[4*i+3]  = a[9*i+6] >> 6;
    r->coeffs[4*i+3] |= (uint32_t)a[9*i+7] << 2;
    r->coeffs[4*i+3] |= (uint32_t)a[9*i+8] << 10;
    r->coeffs[4*i+3] &= 0x3FFFF;

    r->coeffs[4*i+0] = GAMMA - r->coeffs[4*i+0];
    r->coeffs[4*i+1] = GAMMA - r->coeffs[4*i+1];
    r->coeffs[4*i+2] = GAMMA - r->coeffs[4*i+2];
    r->coeffs[4*i+3] = GAMMA - r->coeffs[4*i+3];
  }

  DBENCH_STOP(*tpack);
}


/*************************************************
* Name:        polyq_pack
*
* Description: Bit-pack polynomial with coefficients
*              in [-q/2,q/2].
*
* Arguments:   - uint8_t *r: pointer to output byte array with at least
*                            POLYQ_PACKEDBYTES bytes
*              - const poly *a: pointer to input polynomial
**************************************************/

void polyq_pack(uint8_t *r, const poly *a){
    unsigned int i;
    uint32_t t[8];

    for(i = 0; i < N/8; ++i){
        t[0] = a->coeffs[8 * i + 0] + ((a->coeffs[8 * i + 0]>>31) & Q);
        t[1] = a->coeffs[8 * i + 1] + ((a->coeffs[8 * i + 1]>>31) & Q);
        t[2] = a->coeffs[8 * i + 2] + ((a->coeffs[8 * i + 2]>>31) & Q);
        t[3] = a->coeffs[8 * i + 3] + ((a->coeffs[8 * i + 3]>>31) & Q);
        t[4] = a->coeffs[8 * i + 4] + ((a->coeffs[8 * i + 4]>>31) & Q);
        t[5] = a->coeffs[8 * i + 5] + ((a->coeffs[8 * i + 5]>>31) & Q);
        t[6] = a->coeffs[8 * i + 6] + ((a->coeffs[8 * i + 6]>>31) & Q);
        t[7] = a->coeffs[8 * i + 7] + ((a->coeffs[8 * i + 7]>>31) & Q);

        r[i * 23 + 0] = t[0];
        r[i * 23 + 1] = (t[0] >> 8);  //0: 7 8 8
        r[i * 23 + 2] = ((t[0] >> 16) | (t[1] << 7));

        r[i * 23 + 3] = (t[1] >> 1);
        r[i * 23 + 4] = (t[1] >> 9);  //1: 6 8 8 1
        r[i * 23 + 5] = ((t[1] >> 17) | (t[2] << 6));

        r[i * 23 + 6] = (t[2] >> 2);
        r[i * 23 + 7] = (t[2] >> 10);  //2: 5 8 8 2
        r[i * 23 + 8] = ((t[2] >> 18) | (t[3] << 5));

        r[i * 23 + 9] = (t[3] >> 3);
        r[i * 23 + 10] = (t[3] >> 11);  //3: 4 8 8 3
        r[i * 23 + 11] = ((t[3] >> 19) | (t[4] << 4));

        r[i * 23 + 12] = (t[4] >> 4);
        r[i * 23 + 13] = (t[4] >> 12); //4: 3 8 8 4
        r[i * 23 + 14] = ((t[4] >> 20) | (t[5] << 3));

        r[i * 23 + 15] = (t[5] >> 5);
        r[i * 23 + 16] = (t[5] >> 13);  //5: 2 8 8 5
        r[i * 23 + 17] = ((t[5] >> 21) | (t[6] << 2));

        r[i * 23 + 18] = (t[6] >> 6);
        r[i * 23 + 19] = (t[6] >> 14); //6: 1 8 8 6
        r[i * 23 + 20] = ((t[6] >> 22) | (t[7] << 1));

        r[i * 23 + 21] = (t[7] >> 7);
        r[i * 23 + 22] = (t[7] >> 15);
    }
}

/*************************************************
* Name:        polyq_unpack
*
* Description: Unpack polynomial z with coefficients
*              in R_q.
*
* Arguments:   - poly *r: pointer to output polynomial
*              - const uint8_t *a: byte array with bit-packed polynomial
**************************************************/
void polyq_unpack(poly *r, const uint8_t *a){
    unsigned int i;

    for(i = 0; i < N/8; ++i){
        r->coeffs[8 * i + 0] = a[23 * i + 0];
        r->coeffs[8 * i + 0] |= (uint32_t)a[23 * i + 1] << 8;
        r->coeffs[8 * i + 0] |= (uint32_t)a[23 * i + 2] << 16;
        r->coeffs[8 * i + 0] &= 0x7FFFFF;

        r->coeffs[8 * i + 1] = a[23 * i + 2] >> 7;
        r->coeffs[8 * i + 1] |= (uint32_t)a[23 * i + 3] << 1;
        r->coeffs[8 * i + 1] |= (uint32_t)a[23 * i + 4] << 9;
        r->coeffs[8 * i + 1] |= (uint32_t)a[23 * i + 5] << 17;
        r->coeffs[8 * i + 1] &= 0x7FFFFF;

        r->coeffs[8 * i + 2] = a[23 * i + 5] >> 6;
        r->coeffs[8 * i + 2] |= (uint32_t)a[23 * i + 6] << 2;
        r->coeffs[8 * i + 2] |= (uint32_t)a[23 * i + 7] << 10;
        r->coeffs[8 * i + 2] |= (uint32_t)a[23 * i + 8] << 18;
        r->coeffs[8 * i + 2] &= 0x7FFFFF;

        r->coeffs[8 * i + 3] = a[23 * i + 8] >> 5;
        r->coeffs[8 * i + 3] |= (uint32_t)a[23 * i + 9] << 3;
        r->coeffs[8 * i + 3] |= (uint32_t)a[23 * i + 10] << 11;
        r->coeffs[8 * i + 3] |= (uint32_t)a[23 * i + 11] << 19;
        r->coeffs[8 * i + 3] &= 0x7FFFFF;

        r->coeffs[8 * i + 4] = a[23 * i + 11] >> 4;
        r->coeffs[8 * i + 4] |= (uint32_t)a[23 * i + 12] << 4;
        r->coeffs[8 * i + 4] |= (uint32_t)a[23 * i + 13] << 12;
        r->coeffs[8 * i + 4] |= (uint32_t)a[23 * i + 14] << 20;
        r->coeffs[8 * i + 4] &= 0x7FFFFF;

        r->coeffs[8 * i + 5] = a[23 * i + 14] >> 3;
        r->coeffs[8 * i + 5] |= (uint32_t)a[23 * i + 15] << 5;
        r->coeffs[8 * i + 5] |= (uint32_t)a[23 * i + 16] << 13;
        r->coeffs[8 * i + 5] |= (uint32_t)a[23 * i + 17] << 21;
        r->coeffs[8 * i + 5] &= 0x7FFFFF;

        r->coeffs[8 * i + 6] = a[23 * i + 17] >> 2;
        r->coeffs[8 * i + 6] |= (uint32_t)a[23 * i + 18] << 6;
        r->coeffs[8 * i + 6] |= (uint32_t)a[23 * i + 19] << 14;
        r->coeffs[8 * i + 6] |= (uint32_t)a[23 * i + 20] << 22;
        r->coeffs[8 * i + 6] &= 0x7FFFFF;

        r->coeffs[8 * i + 7] = a[23 * i + 20] >> 1;
        r->coeffs[8 * i + 7] |= (uint32_t)a[23 * i + 21] << 7;
        r->coeffs[8 * i + 7] |= (uint32_t)a[23 * i + 22] << 15;
        r->coeffs[8 * i + 7] &= 0x7FFFFF;

        r->coeffs[8 * i + 0] = r->coeffs[8 * i + 0] - (((Q_HALF - r->coeffs[8 * i + 0])>>31)&Q);
        r->coeffs[8 * i + 1] = r->coeffs[8 * i + 1] - (((Q_HALF - r->coeffs[8 * i + 1])>>31)&Q);
        r->coeffs[8 * i + 2] = r->coeffs[8 * i + 2] - (((Q_HALF - r->coeffs[8 * i + 2])>>31)&Q);
        r->coeffs[8 * i + 3] = r->coeffs[8 * i + 3] - (((Q_HALF - r->coeffs[8 * i + 3])>>31)&Q);
        r->coeffs[8 * i + 4] = r->coeffs[8 * i + 4] - (((Q_HALF - r->coeffs[8 * i + 4])>>31)&Q);
        r->coeffs[8 * i + 5] = r->coeffs[8 * i + 5] - (((Q_HALF - r->coeffs[8 * i + 5])>>31)&Q);
        r->coeffs[8 * i + 6] = r->coeffs[8 * i + 6] - (((Q_HALF - r->coeffs[8 * i + 6])>>31)&Q);
        r->coeffs[8 * i + 7] = r->coeffs[8 * i + 7] - (((Q_HALF - r->coeffs[8 * i + 7])>>31)&Q);

    }

}

/*************************************************
* Name:        polyrho_pack
*
* Description: Bit-pack polynomial with coefficients
*              in S_{\beta * \kappa}.
*
* Arguments:   - uint8_t *r: pointer to output byte array with at least
*                            POLYRHO_PACKEDBYTES bytes
*              - const poly *a: pointer to input polynomial
**************************************************/
void polyrho_pack(uint8_t *r, const poly *a){
    unsigned int i;
    for(i = 0; i < N; ++i){
        r[i] = TAU - a->coeffs[i];
    }
}

/*************************************************
* Name:        polyrho_unpack
*
* Description: Unpack polynomial z with coefficients
*              in S_{\beta * \kappa}.
*
* Arguments:   - poly *r: pointer to output polynomial
*              - const uint8_t *a: byte array with bit-packed polynomial
**************************************************/
void polyrho_unpack(poly *r, const uint8_t *a){
    unsigned int i ;
    for(i = 0; i < N; ++i){
        r->coeffs[i] = TAU - a[i];
    }
}

