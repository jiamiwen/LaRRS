#ifndef NTT_H
#define NTT_H

#include <stdint.h>
#include "params.h"

#define ntt LARRS_NAMESPACE(_ntt)
void ntt(int32_t a[N]);

#define invntt_tomont LARRS_NAMESPACE(_invntt_tomont)
void invntt_tomont(int32_t a[N], int32_t f);

#endif
