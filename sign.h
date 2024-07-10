#ifndef SIGN_H
#define SIGN_H

#include <stddef.h>
#include <stdint.h>
#include <time.h>
#include "params.h"
#include "polyvec.h"
#include "poly.h"

typedef struct SignInParams{
    uint8_t *y_list;
    uint8_t *id_list;
    uint8_t *sk;
    uint8_t *y_TA;
} SIGN_IN_PARAMS;

typedef struct VerifyInParmas{
    uint8_t *y_list;
    uint8_t *id_list;
    uint8_t *y_TA;
} VERIFY_IN_PARAMS;

typedef struct RevokeInParams{
    uint8_t *ctag;
    uint8_t *id_list;
    uint8_t *sk_TA;
} REVOKE_IN_PARMAS;

#define crypto_setup LARRS_NAMESPACE(_crypto_setup)
void crypto_setup(void) ;

#define crypto_sign_keypair LARRS_NAMESPACE(_keypair)
int crypto_sign_keypair(uint8_t *pk_y, uint8_t *pk_id, uint8_t *sk);

#define crypto_sign_keypair1 LARRS_NAMESPACE(_keypair1)
int crypto_sign_keypair1(uint8_t *pk, uint8_t *var,int TA);


#define crypto_sign_signature LARRS_NAMESPACE(_signature)
int crypto_sign_signature(uint8_t *sig, size_t *siglen,
                          const uint8_t *m, size_t mlen,
                          const SIGN_IN_PARAMS *inParams);

#define crypto_sign LARRS_NAMESPACE(_sign)
int crypto_sign(uint8_t *sm, size_t *smlen,
                const uint8_t *m, size_t mlen,
                const SIGN_IN_PARAMS *inParams);

#define crypto_sign_verify LARRS_NAMESPACE(_verify)
int crypto_sign_verify(const uint8_t *sig, size_t siglen,
                       const uint8_t *m, size_t mlen,
                       const VERIFY_IN_PARAMS *inParams);

#define crypto_sign_open LARRS_NAMESPACE(_open)
int crypto_sign_open(size_t mlen,
                     const uint8_t *sm, size_t smlen,
                     const VERIFY_IN_PARAMS *inParams);

#define crypto_revoke LARRS_NAMESPACE(_revoke)
int crypto_revoke(REVOKE_IN_PARMAS *inParmas);
#endif
