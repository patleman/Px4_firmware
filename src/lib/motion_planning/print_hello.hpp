#pragma once
#include<stdio.h>
#include<stdlib.h>
#include<inttypes.h>
typedef uint32_t             mp_digit;
void printing_hello_world();
void s_mp_zero_digs(mp_digit *d, int digits);


 const int MP_MUL_KARATSUBA_CUTOFF = 80;
 const int MP_SQR_KARATSUBA_CUTOFF = 120;
 const int MP_MUL_TOOM_CUTOFF = 350;
 const int MP_SQR_TOOM_CUTOFF = 400;


