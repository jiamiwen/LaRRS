#ifndef REDUCE_H
#define REDUCE_H

#include <stdint.h>
#include "params.h"

#define MONT -4186625 // 2^32 % Q
#define QINV 58728449 // q^(-1) mod 2^32
#define QNINV 4236238847 //-q^(-1) mod 2^32
#define montgomery_reduce LARRS_NAMESPACE(_montgomery_reduce)
int32_t montgomery_reduce(int64_t a);

#define reduce32 LARRS_NAMESPACE(_reduce32)
int32_t reduce32(int32_t a);

#define caddq LARRS_NAMESPACE(_caddq)
int32_t caddq(int32_t a);

#define freeze LARRS_NAMESPACE(_freeze)
int32_t freeze(int32_t a);


#endif
