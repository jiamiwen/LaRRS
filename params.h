#ifndef PARAMS_H
#define PARAMS_H

#include "config.h"

#define SEEDBYTES 32
#define CRHBYTES 48
#define N 256
#define Q 8380417
#define Q_HALF  ((Q + 1) / 2)
#define ROOT_OF_UNITY 1753
#define USER_NUM 4
#define PI 0
#define R2 2365951

#if LARRS_MODE == 1
#define H 3
#define V 3
#define BETA 3
#define KAPPA 30

#define GAMMA (1 << 17)

#elif LARRS_MODE == 2
#define H 4
#define V 4
#define BETA 2
#define KAPPA 39

#define GAMMA (1 << 17)

#endif

#define ETA (GAMMA - KAPPA * BETA)
#define TAU (BETA * KAPPA)

//gamma-kappa*beta, ri,ti
#define POLYETA_PACKEDBYTES   576
//si
#define POLYBETA_PACKEDBYTES  96
//Zq
#define POLYQ_PACKEDBYTES  736
// rho_i
#define POLYRHO_PACKEDBYTES 256

//id,pk
#define Y_BYTES (H * POLYQ_PACKEDBYTES)  //2944
#define CRYPTO_PUBLICKEYBYTES (SEEDBYTES + H * POLYQ_PACKEDBYTES)

// s1 s2
#define CRYPTO_SECRETKEYBYTES ((V + H) * POLYBETA_PACKEDBYTES)


//e1  ri ti  rhoi Cpi T
#define CRYPTO_BYTES (SEEDBYTES + USER_NUM * (V + H + H) * POLYETA_PACKEDBYTES \
                        + USER_NUM * POLYRHO_PACKEDBYTES                       \
                        + (V + V + 1) * POLYQ_PACKEDBYTES + 4)


#define C1_BYTES (POLYQ_PACKEDBYTES * V)
#define C2_BYTES POLYQ_PACKEDBYTES
#define CTAG_BYTES (C1_BYTES + C2_BYTES)

#define EPSILON_BYTES (V * POLYQ_PACKEDBYTES)

#define R_BYTES ((V + H) * POLYETA_PACKEDBYTES)
#define T_BYTES (H * POLYETA_PACKEDBYTES)
#define RESPONSE_BYTES (R_BYTES + T_BYTES + POLYRHO_PACKEDBYTES)

//sig = epsilon||c1||c2||T||e1||r1||t1||rho1||...||rn||tn||rhon
#define CTAG_POS EPSILON_BYTES
#define TIMESTAMP_POS (CTAG_BYTES +  CTAG_POS)
#define TIMESTAMP_BYTES 4
#define E1_POS (TIMESTAMP_POS + TIMESTAMP_BYTES)
#define RESPONSE_POS (E1_POS + SEEDBYTES)
#define T_POS R_BYTES
#define RHO_POS (T_POS + T_BYTES)

#endif
