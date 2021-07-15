#pragma once

#include<limits.h>
#include <stddef.h>
#include <stdint.h>
#include <stdbool.h>
#include<stdlib.h>
#include<stdio.h>
#include <fcntl.h>
#include <errno.h>
#include <unistd.h>
#include<time.h>

////////////// headers from file_operations.h
#include<inttypes.h>
#include <uORB/topics/vehicle_gps_position.h>
#include <uORB/Publication.hpp>
#include<stdarg.h>
#include<string.h>
#include <lib/parameters/param.h>
#include <uORB/topics/sensor_gps.h>

/////////


#include <drivers/drv_hrt.h>
//#include <HealthFlags.h>
#include <lib/parameters/param.h>
#include <systemlib/mavlink_log.h>
#include <uORB/Subscription.hpp>
#include <uORB/topics/vehicle_gps_position.h>
#include <uORB/topics/pa_data.h>


//////////////



#   define MP_2EXPT_C
#   define MP_ABS_C
#   define MP_ADD_C
#   define MP_ADD_D_C
#   define MP_ADDMOD_C
#   define MP_AND_C
#   define MP_CLAMP_C
#   define MP_CLEAR_C
#   define MP_CLEAR_MULTI_C
#   define MP_CMP_C
#   define MP_CMP_D_C
#   define MP_CMP_MAG_C
#   define MP_CNT_LSB_C
#   define MP_COMPLEMENT_C
#   define MP_COPY_C
#   define MP_COUNT_BITS_C
#   define MP_CUTOFFS_C
#   define MP_DIV_C
#   define MP_DIV_2_C
#   define MP_DIV_2D_C
#   define MP_DIV_D_C
#   define MP_DR_IS_MODULUS_C
#   define MP_DR_REDUCE_C
#   define MP_DR_SETUP_C
#   define MP_ERROR_TO_STRING_C
#   define MP_EXCH_C
#   define MP_EXPT_N_C
#   define MP_EXPTMOD_C
#   define MP_EXTEUCLID_C
#   define MP_FREAD_C
#   define MP_FROM_SBIN_C
#   define MP_FROM_UBIN_C
#   define MP_FWRITE_C
#   define MP_GCD_C
#   define MP_GET_DOUBLE_C
#   define MP_GET_I32_C
#   define MP_GET_I64_C
#   define MP_GET_L_C
#   define MP_GET_MAG_U32_C
#   define MP_GET_MAG_U64_C
#   define MP_GET_MAG_UL_C
#   define MP_GROW_C
#   define MP_INIT_C
#   define MP_INIT_COPY_C
#   define MP_INIT_I32_C
#   define MP_INIT_I64_C
#   define MP_INIT_L_C
#   define MP_INIT_MULTI_C
#   define MP_INIT_SET_C
#   define MP_INIT_SIZE_C
#   define MP_INIT_U32_C
#   define MP_INIT_U64_C
#   define MP_INIT_UL_C
#   define MP_INVMOD_C
#   define MP_IS_SQUARE_C
#   define MP_KRONECKER_C
#   define MP_LCM_C
#   define MP_LOG_N_C
#   define MP_LSHD_C
#   define MP_MOD_C
#   define MP_MOD_2D_C
#   define MP_MONTGOMERY_CALC_NORMALIZATION_C
#   define MP_MONTGOMERY_REDUCE_C
#   define MP_MONTGOMERY_SETUP_C
#   define MP_MUL_C
#   define MP_MUL_2_C
#   define MP_MUL_2D_C
#   define MP_MUL_D_C
#   define MP_MULMOD_C
#   define MP_NEG_C
#   define MP_OR_C
#   define MP_PACK_C
#   define MP_PACK_COUNT_C
#   define MP_PRIME_FERMAT_C
#   define MP_PRIME_FROBENIUS_UNDERWOOD_C
#   define MP_PRIME_IS_PRIME_C
#   define MP_PRIME_MILLER_RABIN_C
#   define MP_PRIME_NEXT_PRIME_C
#   define MP_PRIME_RABIN_MILLER_TRIALS_C
#   define MP_PRIME_RAND_C
#   define MP_PRIME_STRONG_LUCAS_SELFRIDGE_C
#   define MP_RADIX_SIZE_C
#   define MP_RADIX_SIZE_OVERESTIMATE_C
#   define MP_RAND_C
#   define MP_RAND_SOURCE_C
#   define MP_READ_RADIX_C
#   define MP_REDUCE_C
#   define MP_REDUCE_2K_C
#   define MP_REDUCE_2K_L_C
#   define MP_REDUCE_2K_SETUP_C
#   define MP_REDUCE_2K_SETUP_L_C
#   define MP_REDUCE_IS_2K_C
#   define MP_REDUCE_IS_2K_L_C
#   define MP_REDUCE_SETUP_C
#   define MP_ROOT_N_C
#   define MP_RSHD_C
#   define MP_SBIN_SIZE_C
#   define MP_SET_C
#   define MP_SET_DOUBLE_C
#   define MP_SET_I32_C
#   define MP_SET_I64_C
#   define MP_SET_L_C
#   define MP_SET_U32_C
#   define MP_SET_U64_C
#   define MP_SET_UL_C
#   define MP_SHRINK_C
#   define MP_SIGNED_RSH_C
#   define MP_SQRMOD_C
#   define MP_SQRT_C
#   define MP_SQRTMOD_PRIME_C
#   define MP_SUB_C
#   define MP_SUB_D_C
#   define MP_SUBMOD_C
#   define MP_TO_RADIX_C
#   define MP_TO_SBIN_C
#   define MP_TO_UBIN_C
#   define MP_UBIN_SIZE_C
#   define MP_UNPACK_C
#   define MP_XOR_C
#   define MP_ZERO_C
#   define S_MP_ADD_C
#   define S_MP_COPY_DIGS_C
#   define S_MP_DIV_3_C
#   define S_MP_DIV_RECURSIVE_C
#   define S_MP_DIV_SCHOOL_C
#   define S_MP_DIV_SMALL_C
#   define S_MP_EXPTMOD_C
#   define S_MP_EXPTMOD_FAST_C
#   define S_MP_GET_BIT_C
#   define S_MP_INVMOD_C
#   define S_MP_INVMOD_ODD_C
#   define S_MP_LOG_C
#   define S_MP_LOG_2EXPT_C
#   define S_MP_LOG_D_C
#   define S_MP_MONTGOMERY_REDUCE_COMBA_C
#   define S_MP_MUL_C
#   define S_MP_MUL_BALANCE_C
#   define S_MP_MUL_COMBA_C
#   define S_MP_MUL_HIGH_C
#   define S_MP_MUL_HIGH_COMBA_C
#   define S_MP_MUL_KARATSUBA_C
#   define S_MP_MUL_TOOM_C
#   define S_MP_PRIME_IS_DIVISIBLE_C
#   define S_MP_PRIME_TAB_C
#   define S_MP_RADIX_MAP_C
#   define S_MP_RADIX_SIZE_OVERESTIMATE_C
#   define S_MP_RAND_PLATFORM_C
#   define S_MP_SQR_C
#   define S_MP_SQR_COMBA_C
#   define S_MP_SQR_KARATSUBA_C
#   define S_MP_SQR_TOOM_C
#   define S_MP_SUB_C
#   define S_MP_ZERO_BUF_C
#   define S_MP_ZERO_DIGS_C


#include "print_hello.hpp"
///////////
using namespace std;
//typedef uint32_t             mp_digit;
typedef uint64_t mp_word;

typedef struct{
   int used,alloc,sign;
   mp_digit *dp;
}mp_int;

#define MP_DIGIT_BIT 28

#define MP_MASK          ((((mp_digit)1)<<((mp_digit)MP_DIGIT_BIT))-((mp_digit)1))
#define MP_DIGIT_MAX     MP_MASK
/*
static const int MP_MUL_KARATSUBA_CUTOFF = 80;
static const int MP_SQR_KARATSUBA_CUTOFF = 120;
static const int MP_MUL_TOOM_CUTOFF = 350;
static const int MP_SQR_TOOM_CUTOFF = 400;
*/



/* Primality generation flags */
#define MP_PRIME_BBS      0x0001 /* BBS style prime */
#define MP_PRIME_SAFE     0x0002 /* Safe prime (p-1)/2 == prime */
#define MP_PRIME_2MSB_ON  0x0008 /* force 2nd MSB to 1 */
#define MP_MALLOC(size)                   malloc(size)
#define MP_REALLOC(mem, oldsize, newsize) realloc((mem), (newsize))
#define MP_CALLOC(nmemb, size)            calloc((nmemb), (size))
#define MP_FREE(mem, size)                free(mem)
#define MP_DEV_URANDOM "/dev/urandom"
#define mp_iszero(a) ((a)->used == 0)
#define mp_isneg(a)  ((a)->sign == MP_NEG)
#define mp_iseven(a) (((a)->used == 0) || (((a)->dp[0] & 1u) == 0u))
#define mp_isodd(a)  (!mp_iseven(a))


//#  define MP_FREE_BUF(mem, size)   MP_FREE((mem), (size))
#define MP_STRINGIZE(x)  MP__STRINGIZE(x)
#define MP__STRINGIZE(x) ""#x""
#define MP_HAS(x)        (sizeof(MP_STRINGIZE(x##_C)) == 1u)

#define MP_MIN(x, y) (((x) < (y)) ? (x) : (y))
#define MP_MAX(x, y) (((x) > (y)) ? (x) : (y))
#define MP_EXCH(t, a, b) do { t _c = a; a = b; b = _c; } while (0)
//#define CHAR_BIT 8;
#define MP_SIZEOF_BITS(type)    ((size_t)CHAR_BIT * sizeof(type))

#define MP_MAX_COMBA            (int)(1uL << (MP_SIZEOF_BITS(mp_word) - (2u * (size_t)MP_DIGIT_BIT)))
#define MP_WARRAY               (int)(1uL << ((MP_SIZEOF_BITS(mp_word) - (2u * (size_t)MP_DIGIT_BIT)) + 1u))

#define MP_TOUPPER(c) ((((c) >= 'a') && ((c) <= 'z')) ? (((c) + 'A') - 'a') : (c))



/*
#ifndef passing
#define passing
void s_mp_zero_digs(mp_digit *d, int digits);
void s_mp_zero_digs(mp_digit *d, int digits)
{
   while (digits-- > 0) {
      *d++ = 0;
   }
}
#endif

*/
#  define MP_FREE_BUF(mem, size)                        \
do {                                                    \
   size_t fs_ = (size);                                 \
   void* fm_ = (mem);                                   \
   if (fm_ != NULL) {                                   \
      s_mp_zero_buf(fm_, fs_);                          \
      MP_FREE(fm_, fs_);                                \
   }                                                    \
} while (0)
#  define MP_FREE_DIGS(mem, digits)                     \
do {                                                    \
   int fd_ = (digits);                                  \
   mp_digit* fm_ = (mem);                               \
   if (fm_ != NULL) {                                   \
      s_mp_zero_digs(fm_, fd_);                         \
      MP_FREE(fm_, sizeof (mp_digit) * (size_t)fd_);    \
   }                                                    \
} while (0)



//#define MP_EXCH(t, a, b) do { t _c = a; a = b; b = _c; } while (0)

#define MP_IS_2EXPT(x) (((x) != 0u) && (((x) & ((x) - 1u)) == 0u))
/* b = a*a  */
#define mp_sqr(a, b) mp_mul((a), (a), (b))
#define MP_MAX_DIGIT_COUNT ((INT_MAX - 2) / MP_DIGIT_BIT)
#define MP_MIN_DIGIT_COUNT MP_MAX(3, (((int)MP_SIZEOF_BITS(uint64_t) + MP_DIGIT_BIT) - 1) / MP_DIGIT_BIT)
#   define TAB_SIZE 256
#   define MAX_WINSIZE 0

#define MP_DEFAULT_DIGIT_COUNT 32

#define MP_RADIX_MAP_REVERSE_SIZE 80u

#define MP_SET_SIGNED(name, uname, type, utype)          \
    void name(mp_int * a, type b)                        \
    {                                                    \
        uname(a, (b < 0) ? -(utype)b : (utype)b);        \
        if (b < 0) { a->sign = MP_NEG; }                 \
    }



/* code-generating macros */
#define MP_SET_UNSIGNED(name, type)                                                    \
    void name(mp_int * a, type b)                                                      \
    {                                                                                  \
        int i = 0;                                                                     \
        while (b != 0u) {                                                              \
            a->dp[i++] = ((mp_digit)b & MP_MASK);                                      \
            if (MP_SIZEOF_BITS(type) <= MP_DIGIT_BIT) { break; }                       \
            b >>= ((MP_SIZEOF_BITS(type) <= MP_DIGIT_BIT) ? 0 : MP_DIGIT_BIT);         \
        }                                                                              \
        a->used = i;                                                                   \
        a->sign = MP_ZPOS;                                                             \
        s_mp_zero_digs(a->dp + a->used, a->alloc - a->used);                         \
    }

typedef enum {
   MP_ZPOS = 0,   /* positive */
   MP_NEG = 1     /* negative */
} mp_sign;

typedef enum {
   MP_LT = -1,    /* less than */
   MP_EQ = 0,     /* equal */
   MP_GT = 1      /* greater than */
} mp_ord;

typedef enum {
   MP_OKAY  = 0,   /* no error */
   MP_ERR   = -1,  /* unknown error */
   MP_MEM   = -2,  /* out of mem */
   MP_VAL   = -3,  /* invalid input */
   MP_ITER  = -4,  /* maximum iterations reached */
   MP_BUF   = -5,  /* buffer overflow, supplied buffer too small */
   MP_OVF   = -6   /* mp_int overflow, too many digits */
} mp_err;

//////////////////////////////////////////



/* ---> init and deinit bignum functions <--- */
/* init a bignum */
mp_err mp_init(mp_int *a) ;

/* free a bignum */
void mp_clear(mp_int *a);

/* init a null terminated series of arguments */
mp_err mp_init_multi(mp_int *mp, ...) ;

/* clear a null terminated series of arguments */
void mp_clear_multi(mp_int *mp, ...) ;

/* exchange two ints */
void mp_exch(mp_int *a, mp_int *b);

/* shrink ram required for a bignum */
mp_err mp_shrink(mp_int *a) ;

/* grow an int to a given size */
mp_err mp_grow(mp_int *a, int size) ;

/* init to a given number of digits */
mp_err mp_init_size(mp_int *a, int size);

/* ---> Basic Manipulations <--- */
#define mp_iszero(a) ((a)->used == 0)
#define mp_isneg(a)  ((a)->sign == MP_NEG)
#define mp_iseven(a) (((a)->used == 0) || (((a)->dp[0] & 1u) == 0u))
#define mp_isodd(a)  (!mp_iseven(a))

/* set to zero */
void mp_zero(mp_int *a);



/* get integer, set integer and init with integer (int32_t) */
int32_t mp_get_i32(const mp_int *a) ;
void mp_set_i32(mp_int *a, int32_t b);
mp_err mp_init_i32(mp_int *a, int32_t b) ;

/* get integer, set integer and init with integer, behaves like two complement for negative numbers (uint32_t) */
#define mp_get_u32(a) ((uint32_t)mp_get_i32(a))
void mp_set_u32(mp_int *a, uint32_t b);
mp_err mp_init_u32(mp_int *a, uint32_t b);

/* get integer, set integer and init with integer (int64_t) */
int64_t mp_get_i64(const mp_int *a) ;
void mp_set_i64(mp_int *a, int64_t b);
mp_err mp_init_i64(mp_int *a, int64_t b) ;

/* get integer, set integer and init with integer, behaves like two complement for negative numbers (uint64_t) */
#define mp_get_u64(a) ((uint64_t)mp_get_i64(a))
void mp_set_u64(mp_int *a, uint64_t b);
mp_err mp_init_u64(mp_int *a, uint64_t b) ;

/* get magnitude */
uint32_t mp_get_mag_u32(const mp_int *a) ;
uint64_t mp_get_mag_u64(const mp_int *a) ;
unsigned long mp_get_mag_ul(const mp_int *a) ;

/* get integer, set integer (long) */
long mp_get_l(const mp_int *a) ;
void mp_set_l(mp_int *a, long b);
mp_err mp_init_l(mp_int *a, long b) ;

/* get integer, set integer (unsigned long) */
#define mp_get_ul(a) ((unsigned long)mp_get_l(a))
void mp_set_ul(mp_int *a, unsigned long b);
mp_err mp_init_ul(mp_int *a, unsigned long b) ;

/* set to single unsigned digit, up to MP_DIGIT_MAX */
void mp_set(mp_int *a, mp_digit b);
mp_err mp_init_set(mp_int *a, mp_digit b) ;

/* copy, b = a */
mp_err mp_copy(const mp_int *a, mp_int *b) ;

/* inits and copies, a = b */
mp_err mp_init_copy(mp_int *a, const mp_int *b) ;

/* trim unused digits */
void mp_clamp(mp_int *a);



/* pack binary data */
size_t mp_pack_count(const mp_int *a, size_t nails, size_t size) ;


/* ---> digit manipulation <--- */

/* right shift by "b" digits */
void mp_rshd(mp_int *a, int b);

/* left shift by "b" digits */
mp_err mp_lshd(mp_int *a, int b) ;

/* c = a / 2**b, implemented as c = a >> b */
mp_err mp_div_2d(const mp_int *a, int b, mp_int *c, mp_int *d) ;

/* b = a/2 */
mp_err mp_div_2(const mp_int *a, mp_int *b) ;

/* c = a * 2**b, implemented as c = a << b */
mp_err mp_mul_2d(const mp_int *a, int b, mp_int *c) ;

/* b = a*2 */
mp_err mp_mul_2(const mp_int *a, mp_int *b) ;

/* c = a mod 2**b */
mp_err mp_mod_2d(const mp_int *a, int b, mp_int *c) ;

/* computes a = 2**b */
mp_err mp_2expt(mp_int *a, int b) ;

/* Counts the number of lsbs which are zero before the first zero bit */
int mp_cnt_lsb(const mp_int *a) ;

/* I Love Earth! */

/* makes a pseudo-random mp_int of a given size */
mp_err mp_rand(mp_int *a, int digits) ;
/* use custom random data source instead of source provided the platform */
void mp_rand_source(mp_err(*source)(void *out, size_t size));

/* ---> binary operations <--- */

/* c = a XOR b (two complement) */
mp_err mp_xor(const mp_int *a, const mp_int *b, mp_int *c) ;

/* c = a OR b (two complement) */
mp_err mp_or(const mp_int *a, const mp_int *b, mp_int *c) ;

/* c = a AND b (two complement) */
mp_err mp_and(const mp_int *a, const mp_int *b, mp_int *c) ;

/* b = ~a (bitwise not, two complement) */
mp_err mp_complement(const mp_int *a, mp_int *b) ;

/* right shift with sign extension */
mp_err mp_signed_rsh(const mp_int *a, int b, mp_int *c) ;

/* ---> Basic arithmetic <--- */

/* b = -a */
mp_err mp_neg(const mp_int *a, mp_int *b) ;

/* b = |a| */
mp_err mp_abs(const mp_int *a, mp_int *b) ;

/* compare a to b */
mp_ord mp_cmp(const mp_int *a, const mp_int *b) ;

/* compare |a| to |b| */
mp_ord mp_cmp_mag(const mp_int *a, const mp_int *b) ;

/* c = a + b */
mp_err mp_add(const mp_int *a, const mp_int *b, mp_int *c) ;

/* c = a - b */
mp_err mp_sub(const mp_int *a, const mp_int *b, mp_int *c) ;
/* c = a * b */
mp_err mp_mul(const mp_int *a, const mp_int *b, mp_int *c) ;

/* b = a*a  */
#define mp_sqr(a, b) mp_mul((a), (a), (b))

/* a/b => cb + d == a */
mp_err mp_div(const mp_int *a, const mp_int *b, mp_int *c, mp_int *d) ;

/* c = a mod b, 0 <= c < b  */
mp_err mp_mod(const mp_int *a, const mp_int *b, mp_int *c) ;

/* Increment "a" by one like "a++". Changes input! */
#define mp_incr(a) mp_add_d((a), 1u, (a))

/* Decrement "a" by one like "a--". Changes input! */
#define mp_decr(a) mp_sub_d((a), 1u, (a))

/* ---> single digit functions <--- */

/* compare against a single digit */
mp_ord mp_cmp_d(const mp_int *a, mp_digit b);

/* c = a + b */
mp_err mp_add_d(const mp_int *a, mp_digit b, mp_int *c) ;

/* c = a - b */
mp_err mp_sub_d(const mp_int *a, mp_digit b, mp_int *c) ;

/* c = a * b */
mp_err mp_mul_d(const mp_int *a, mp_digit b, mp_int *c) ;

/* a/b => cb + d == a */
mp_err mp_div_d(const mp_int *a, mp_digit b, mp_int *c, mp_digit *d) ;

/* c = a mod b, 0 <= c < b  */
#define mp_mod_d(a, b, c) mp_div_d((a), (b), NULL, (c))

/* ---> number theory <--- */

/* d = a + b (mod c) */
mp_err mp_addmod(const mp_int *a, const mp_int *b, const mp_int *c, mp_int *d) ;

/* d = a - b (mod c) */
mp_err mp_submod(const mp_int *a, const mp_int *b, const mp_int *c, mp_int *d) ;

/* d = a * b (mod c) */
mp_err mp_mulmod(const mp_int *a, const mp_int *b, const mp_int *c, mp_int *d) ;

/* c = a * a (mod b) */
mp_err mp_sqrmod(const mp_int *a, const mp_int *b, mp_int *c) ;

/* c = 1/a (mod b) */
mp_err mp_invmod(const mp_int *a, const mp_int *b, mp_int *c) ;

/* c = (a, b) */
mp_err mp_gcd(const mp_int *a, const mp_int *b, mp_int *c) ;

/* produces value such that U1*a + U2*b = U3 */
mp_err mp_exteuclid(const mp_int *a, const mp_int *b, mp_int *U1, mp_int *U2, mp_int *U3) ;

/* c = [a, b] or (a*b)/(a, b) */
mp_err mp_lcm(const mp_int *a, const mp_int *b, mp_int *c) ;

/* Integer logarithm to integer base */
mp_err mp_log_n(const mp_int *a, int base, int *c) ;

/* c = a**b */
mp_err mp_expt_n(const mp_int *a, int b, mp_int *c) ;

/* finds one of the b'th root of a, such that |c|**b <= |a|
 *
 * returns error if a < 0 and b is even
 */
mp_err mp_root_n(const mp_int *a, int b, mp_int *c) ;

/* special sqrt algo */
mp_err mp_sqrt(const mp_int *arg, mp_int *ret) ;

/* special sqrt (mod prime) */
mp_err mp_sqrtmod_prime(const mp_int *n, const mp_int *prime, mp_int *ret) ;

/* is number a square? */
mp_err mp_is_square(const mp_int *arg, bool *ret) ;

/* computes the Kronecker symbol c = (a | p) (like jacobi() but with {a,p} in Z */
mp_err mp_kronecker(const mp_int *a, const mp_int *p, int *c) ;

/* used to setup the Barrett reduction for a given modulus b */
mp_err mp_reduce_setup(mp_int *a, const mp_int *b) ;

/* Barrett Reduction, computes a (mod b) with a precomputed value c
 *
 * Assumes that 0 < x <= m*m, note if 0 > x > -(m*m) then you can merely
 * compute the reduction as -1 * mp_reduce(mp_abs(x)) [pseudo code].
 */
mp_err mp_reduce(mp_int *x, const mp_int *m, const mp_int *mu) ;

/* setups the montgomery reduction */
mp_err mp_montgomery_setup(const mp_int *n, mp_digit *rho) ;

/* computes a = B**n mod b without division or multiplication useful for
 * normalizing numbers in a Montgomery system.
 */
mp_err mp_montgomery_calc_normalization(mp_int *a, const mp_int *b) ;

/* computes x/R == x (mod N) via Montgomery Reduction */
mp_err mp_montgomery_reduce(mp_int *x, const mp_int *n, mp_digit rho) ;

/* returns 1 if a is a valid DR modulus */
bool mp_dr_is_modulus(const mp_int *a) ;

/* sets the value of "d" required for mp_dr_reduce */
void mp_dr_setup(const mp_int *a, mp_digit *d);

/* reduces a modulo n using the Diminished Radix method */
mp_err mp_dr_reduce(mp_int *x, const mp_int *n, mp_digit k) ;

/* returns true if a can be reduced with mp_reduce_2k */
bool mp_reduce_is_2k(const mp_int *a) ;

/* determines k value for 2k reduction */
mp_err mp_reduce_2k_setup(const mp_int *a, mp_digit *d) ;

/* reduces a modulo b where b is of the form 2**p - k [0 <= a] */
mp_err mp_reduce_2k(mp_int *a, const mp_int *n, mp_digit d) ;

/* returns true if a can be reduced with mp_reduce_2k_l */
bool mp_reduce_is_2k_l(const mp_int *a) ;

/* determines k value for 2k reduction */
mp_err mp_reduce_2k_setup_l(const mp_int *a, mp_int *d) ;

/* reduces a modulo b where b is of the form 2**p - k [0 <= a] */
mp_err mp_reduce_2k_l(mp_int *a, const mp_int *n, const mp_int *d) ;

/* Y = G**X (mod P) */
mp_err mp_exptmod(const mp_int *G, const mp_int *X, const mp_int *P, mp_int *Y) ;

/* ---> Primes <--- */

/* performs one Fermat test of "a" using base "b".
 * Sets result to 0 if composite or 1 if probable prime
 */
mp_err mp_prime_fermat(const mp_int *a, const mp_int *b, bool *result) ;

/* performs one Miller-Rabin test of "a" using base "b".
 * Sets result to 0 if composite or 1 if probable prime
 */
mp_err mp_prime_miller_rabin(const mp_int *a, const mp_int *b, bool *result) ;

/* This gives [for a given bit size] the number of trials required
 * such that Miller-Rabin gives a prob of failure lower than 2^-96
 */
int mp_prime_rabin_miller_trials(int size);

/* performs one strong Lucas-Selfridge test of "a".
 * Sets result to 0 if composite or 1 if probable prime
 */
mp_err mp_prime_strong_lucas_selfridge(const mp_int *a, bool *result) ;

/* performs one Frobenius test of "a" as described by Paul Underwood.
 * Sets result to 0 if composite or 1 if probable prime
 */
mp_err mp_prime_frobenius_underwood(const mp_int *N, bool *result) ;

/* performs t random rounds of Miller-Rabin on "a" additional to
 * bases 2 and 3.  Also performs an initial sieve of trial
 * division.  Determines if "a" is prime with probability
 * of error no more than (1/4)**t.
 * Both a strong Lucas-Selfridge to complete the BPSW test
 * and a separate Frobenius test are available at compile time.
 * With t<0 a deterministic test is run for primes up to
 * 318665857834031151167461. With t<13 (abs(t)-13) additional
 * tests with sequential small primes are run starting at 43.
 * Is Fips 186.4 compliant if called with t as computed by
 * mp_prime_rabin_miller_trials();
 *
 * Sets result to 1 if probably prime, 0 otherwise
 */
mp_err mp_prime_is_prime(const mp_int *a, int t, bool *result) ;

/* finds the next prime after the number "a" using "t" trials
 * of Miller-Rabin.
 *
 * bbs_style = true means the prime must be congruent to 3 mod 4
 */
mp_err mp_prime_next_prime(mp_int *a, int t, bool bbs_style) ;

/* makes a truly random prime of a given size (bits),
 *
 * Flags are as follows:
 *
 *   MP_PRIME_BBS      - make prime congruent to 3 mod 4
 *   MP_PRIME_SAFE     - make sure (p-1)/2 is prime as well (implies MP_PRIME_BBS)
 *   MP_PRIME_2MSB_ON  - make the 2nd highest bit one
 *
 * You have to supply a callback which fills in a buffer with random bytes.  "dat" is a parameter you can
 * have passed to the callback (e.g. a state or something).  This function doesn't use "dat" itself
 * so it can be NULL
 *
 */
mp_err mp_prime_rand(mp_int *a, int t, int size, int flags) ;

/* ---> radix conversion <--- */
int mp_count_bits(const mp_int *a) ;

size_t mp_ubin_size(const mp_int *a) ;
mp_err mp_from_ubin(mp_int *a, const uint8_t *buf, size_t size) ;
mp_err mp_to_ubin(const mp_int *a, uint8_t *buf, size_t maxlen, size_t *written) ;

size_t mp_sbin_size(const mp_int *a) ;
mp_err mp_from_sbin(mp_int *a, const uint8_t *buf, size_t size) ;
mp_err mp_to_sbin(const mp_int *a, uint8_t *buf, size_t maxlen, size_t *written) ;

mp_err mp_read_radix(mp_int *a, const char *str, int radix) ;
mp_err mp_to_radix(const mp_int *a, char *str, size_t maxlen, size_t *written, int radix) ;

mp_err mp_radix_size(const mp_int *a, int radix, size_t *size) ;
mp_err mp_radix_size_overestimate(const mp_int *a, const int radix, size_t *size) ;

#ifndef MP_NO_FILE
mp_err mp_fread(mp_int *a, int radix, FILE *stream) ;
mp_err mp_fwrite(const mp_int *a, int radix, FILE *stream) ;
#endif

#define mp_to_binary(M, S, N)  mp_to_radix((M), (S), (N), NULL, 2)
#define mp_to_octal(M, S, N)   mp_to_radix((M), (S), (N), NULL, 8)
#define mp_to_decimal(M, S, N) mp_to_radix((M), (S), (N), NULL, 10)
#define mp_to_hex(M, S, N)     mp_to_radix((M), (S), (N), NULL, 16)


/* lowlevel functions, do not call! */
 bool s_mp_get_bit(const mp_int *a, int b) ;
 int s_mp_log_2expt(const mp_int *a, mp_digit base) ;
 int s_mp_log_d(mp_digit base, mp_digit n) ;
 mp_err s_mp_add(const mp_int *a, const mp_int *b, mp_int *c) ;
 mp_err s_mp_div_3(const mp_int *a, mp_int *c, mp_digit *d) ;
 mp_err s_mp_div_recursive(const mp_int *a, const mp_int *b, mp_int *q, mp_int *r) ;
mp_err s_mp_div_school(const mp_int *a, const mp_int *b, mp_int *c, mp_int *d) ;
 mp_err s_mp_div_small(const mp_int *a, const mp_int *b, mp_int *c, mp_int *d) ;
 mp_err s_mp_exptmod(const mp_int *G, const mp_int *X, const mp_int *P, mp_int *Y, int redmode) ;
 mp_err s_mp_exptmod_fast(const mp_int *G, const mp_int *X, const mp_int *P, mp_int *Y, int redmode) ;
 mp_err s_mp_invmod(const mp_int *a, const mp_int *b, mp_int *c) ;
 mp_err s_mp_invmod_odd(const mp_int *a, const mp_int *b, mp_int *c) ;
 mp_err s_mp_log(const mp_int *a, mp_digit base, int *c) ;
 mp_err s_mp_montgomery_reduce_comba(mp_int *x, const mp_int *n, mp_digit rho) ;
 mp_err s_mp_mul(const mp_int *a, const mp_int *b, mp_int *c, int digs) ;
 mp_err s_mp_mul_balance(const mp_int *a, const mp_int *b, mp_int *c) ;
 mp_err s_mp_mul_comba(const mp_int *a, const mp_int *b, mp_int *c, int digs) ;
 mp_err s_mp_mul_high(const mp_int *a, const mp_int *b, mp_int *c, int digs) ;
 mp_err s_mp_mul_high_comba(const mp_int *a, const mp_int *b, mp_int *c, int digs) ;
 mp_err s_mp_mul_karatsuba(const mp_int *a, const mp_int *b, mp_int *c) ;
 mp_err s_mp_mul_toom(const mp_int *a, const mp_int *b, mp_int *c) ;
 mp_err s_mp_prime_is_divisible(const mp_int *a, bool *result) ;
 mp_err s_mp_rand_platform(void *p, size_t n) ;
 mp_err s_mp_sqr(const mp_int *a, mp_int *b) ;
 mp_err s_mp_sqr_comba(const mp_int *a, mp_int *b) ;
 mp_err s_mp_sqr_karatsuba(const mp_int *a, mp_int *b) ;
 mp_err s_mp_sqr_toom(const mp_int *a, mp_int *b) ;
 mp_err s_mp_sub(const mp_int *a, const mp_int *b, mp_int *c) ;
 void s_mp_copy_digs(mp_digit *d, const mp_digit *s, int digits);
 void s_mp_zero_buf(void *mem, size_t size);

 mp_err s_mp_radix_size_overestimate(const mp_int *a, const int radix, size_t *size);


///////////////////////////////

static const char *const BASE64_DIGITS ="ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";

static const char *const HEX_DIGITS = "0123456789abcdef";

class SHA256 {

public:
	SHA256();
	void update(uint8_t * data, size_t length);
	uint8_t * digest();


private:
	uint8_t  m_data[64];
	uint32_t m_blocklen;
	uint64_t m_bitlen;
	uint32_t m_state[8]; //A, B, C, D, E, F, G, H

    uint32_t K[64] = {
		0x428a2f98,0x71374491,0xb5c0fbcf,0xe9b5dba5,
		0x3956c25b,0x59f111f1,0x923f82a4,0xab1c5ed5,
		0xd807aa98,0x12835b01,0x243185be,0x550c7dc3,
		0x72be5d74,0x80deb1fe,0x9bdc06a7,0xc19bf174,
		0xe49b69c1,0xefbe4786,0x0fc19dc6,0x240ca1cc,
		0x2de92c6f,0x4a7484aa,0x5cb0a9dc,0x76f988da,
		0x983e5152,0xa831c66d,0xb00327c8,0xbf597fc7,
		0xc6e00bf3,0xd5a79147,0x06ca6351,0x14292967,
		0x27b70a85,0x2e1b2138,0x4d2c6dfc,0x53380d13,
		0x650a7354,0x766a0abb,0x81c2c92e,0x92722c85,
		0xa2bfe8a1,0xa81a664b,0xc24b8b70,0xc76c51a3,
		0xd192e819,0xd6990624,0xf40e3585,0x106aa070,
		0x19a4c116,0x1e376c08,0x2748774c,0x34b0bcb5,
		0x391c0cb3,0x4ed8aa4a,0x5b9cca4f,0x682e6ff3,
		0x748f82ee,0x78a5636f,0x84c87814,0x8cc70208,
		0x90befffa,0xa4506ceb,0xbef9a3f7,0xc67178f2
	};

	static uint32_t rotr(uint32_t x, uint32_t n);
	static uint32_t choose(uint32_t e, uint32_t f, uint32_t g);
	static uint32_t majority(uint32_t a, uint32_t b, uint32_t c);
	static uint32_t sig0(uint32_t x);
	static uint32_t sig1(uint32_t x);
	void transform();
	void pad();
	void revert(uint8_t * hash);

};

////////////////////////////////
struct Date_time{
       int Hours;
       int Minutes;
       int Seconds;
       int date;
       int Month;
       int Year;
};
struct pair_set{
    char tag[30];
    char value[1200];
};

struct key{
char modulus[1000];
char private_exponent[1000];
};


static const char HEX_ele[17]="0123456789abcdef";

void fetch_tag(char *content, char *target, char *result);

int isSubstring(char *s1, char *s2);///s1 is the sub string ; s2 is the larger string

int find_int_hex(char ptr);

void base64Encoder(char input_str[], int len_str, char *result);

// signs the content with signKey(a key) and puts the result inside result
void signing_support(key signKey,char *content,char*result);

void pair_file_write(pair_set *ptr,int ptr_quant ,char *file_name , key Skey);


//this function is for validating the txt file wrt  the given key
int inBase64(char *d);

void base64decoder(char *base64_ptr,char *hex);

char small_letter(char a);

int Validating_File(char *file , key key);

// this function is for encrypting files that has to be remained inside the RFM
void encrypting_File(char *content, key key, char *fname);


// this function is for decrypting files that has to be remained inside the RFM
void decrypting_File(char *file , key key, char *result);

// this function is for fetching value of a particular tag in the given string
void fetch_tag(char *content, char *target, char *result);


//function to generate DroneID.txt
void DroneIDcreation( );

// this function  creates .txt file when an amendment takes
// place in the hardware(currently only gps)
void HardwareInuseCreation(int gps);


void get_RFM_Key(key *key1);


/// this function combines validation and fetching of tag values of files
// coming from MC/MS  or staying inside RFM
int file_read(char *fname,char *tag,char *result,int file_type);


//amendment in KeyLog.txt with fileID of a file
void KeyLog_Regen(char *fileID);




// this function will tell if the current time lies inside the start and end
bool In_Time(Date_time current, Date_time start, Date_time end);


// this function will convert string format of time to Date_time struct
void conversion_to_dateTime(Date_time *dt,char *str);

//this function will tell if current time lies inside the start and end time
int check_time(char *start,char *end);

int check_recentPA(char *paID);


//////////////////////////////////////// headers from PA_EXTRACT.hpp




struct stack{
    char list_att[40][40];
    int counter;
};

struct tag_value{
    char tag[100];
    char value[100];
};



struct GEO_DATE_TIME_XML{// for storing xml data
    Date_time start;
    Date_time end;
};

struct ParsedData{
   char start_time[50];
   char end_time[50];
   double long_lat_coords[20][10];
   int num_coords;

};


void updateStack_remove(char *ptr, stack *ss);

void updateStack_add(char *ptr,int set, stack *ss);

void replaceTag(char *ptr, stack *ss,char *result,int *ptr_count);

int isSignature(char *ptr);

void Reference_canon(char *file_name, char *res);

void cleanerXML(char *input,char *output);

void SignedInfo_canon(char *file_name, char *res);

void getTagvalue(char *Tag, char *certificate,char *file_nam_perm);

char* strip(char *str);

void formPairstrings(char *tag,char *value,int st1,int st2,int end1,int end2,char *sample);

int min(int a,int b);

void sortAndmake(tag_value *aux,int number_of_tags,char *result);

void xmlExeclusiveCanon(char *sample,char *result);


int Is_PA_VAlid();

ParsedData parse_artifact();

bool In_Place(ParsedData geo,int latti, int longi);

bool In_Time(Date_time current, GEO_DATE_TIME_XML Xml);

void data_fetch_delete();


int date_time_extract_and_check();

int DroneIDverification(char *paID);

int DroneIDverification();
////////////////////////////////////

//function for maintaining hardware integrity
//return types
//0:not genuine harware used
//1 : genuine hardware used
int RPAS_identifier();

void   call_DroneIDcreation();


int call_file_read(char *fname,char *tag,char *result,int file_type);

int CHECK_REUSAGE(char *file_id);

//for calling KeyLog_Regen
void call_KeyLog_Regen(char *fileID);



// Key rotation function
//two files are generated
//1)PublicKeyNew.txt: later sent to MS
//2)PublicPrivateNew.txt:kept inside RFM encrypted
void Key_rotation_start(char *file_id);



/*In this function contents of PublicPrivateNew.txt are fetched. The corresponding tags and values
inside the PublicPrivateInuse.txt are replaced by the ones that are fetched earlier.
The process invloves 1)decrypting of PublicPrivateNew.txt
                     2)fetching Modulus, Privatekey, PublicKey
                     3)making new PublicPrivateInuse.txt
*/
void KEY_CHANGE_INITIATION();

/*This function comes in handy when ParamChangePerm.txt is present and one
has to modify already existig ParamInuse.txt file
read ParamChangePerm.txt
if any similar then change the value
if any new then add the value
*/
void ParamInuseModify();


/*
function to start setting the parameters as specified in ParamInuse.txt
*/
void ParamSetfile();

//this function will make the amendments in recentPA.txt file
/*
fetch_required
freq_done
previous_log_hash
*/
void update_recentPA(int type,char *value);
