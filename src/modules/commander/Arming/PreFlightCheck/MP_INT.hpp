/* steps for key generation

1) generation of two prime numbers  p,q (024)

2) multiplication of two primes, modulus,  n=p*q

3) phi(n) = (p-1)*(q-1)

4)  e=65537 public exponent

5) d=(1/e) mod phi(n)

6) Private key = (d,n)
   Public Key = (e,n)

*/
#pragma once

#include<limits.h>
#include <stddef.h>
#include <stdint.h>
#include <stdbool.h>
#include<stdlib.h>
#include<stdio.h>
#include<inttypes.h>
#include <fcntl.h>
#include <errno.h>
#include <unistd.h>
#include<time.h>
#include<stdarg.h>

/////////



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



///////////
using namespace std;
typedef uint32_t             mp_digit;
typedef uint64_t mp_word;

typedef struct{
   int used,alloc,sign;
   mp_digit *dp;
}mp_int;

#define MP_DIGIT_BIT 28

#define MP_MASK          ((((mp_digit)1)<<((mp_digit)MP_DIGIT_BIT))-((mp_digit)1))
#define MP_DIGIT_MAX     MP_MASK

int MP_MUL_KARATSUBA_CUTOFF = 80;
int MP_SQR_KARATSUBA_CUTOFF = 120;
int MP_MUL_TOOM_CUTOFF = 350;
int MP_SQR_TOOM_CUTOFF = 400;




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


void s_mp_zero_digs(mp_digit *d, int digits)
{
   while (digits-- > 0) {
      *d++ = 0;
   }
}

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
 void s_mp_zero_digs(mp_digit *d, int digits);
 mp_err s_mp_radix_size_overestimate(const mp_int *a, const int radix, size_t *size);


//////////////////////////////////////////


/* reverse an array, used for radix code */
static void s_reverse(char *s, size_t len)
{
   size_t ix = 0, iy = len - 1u;
   while (ix < iy) {
      MP_EXCH(char, s[ix], s[iy]);
      ++ix;
      --iy;
   }
}
const char s_mp_radix_map[] = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz+/";
/* stores a bignum as a ASCII string in a given radix (2..64)
 *
 * Stores upto "size - 1" chars and always a NULL byte, puts the number of characters
 * written, including the '\0', in "written".
 */
mp_err mp_to_radix(const mp_int *a, char *str, size_t maxlen, size_t *written, int radix)
{
   size_t  digs;
   mp_err  err;
   mp_int  t;
   mp_digit d;
   char   *_s = str;

   /* check range of radix and size*/
   if (maxlen < 2u) {
      return MP_BUF;
   }
   if ((radix < 2) || (radix > 64)) {
      return MP_VAL;
   }

   /* quick out if its zero */
   if (mp_iszero(a)) {
      *str++ = '0';
      *str = '\0';
      if (written != NULL) {
         *written = 2u;
      }
      return MP_OKAY;
   }

   if ((err = mp_init_copy(&t, a)) != MP_OKAY) {
      return err;
   }

   /* if it is negative output a - */
   if (mp_isneg(&t)) {
      /* we have to reverse our digits later... but not the - sign!! */
      ++_s;

      /* store the flag and mark the number as positive */
      *str++ = '-';
      t.sign = MP_ZPOS;

      /* subtract a char */
      --maxlen;
   }
   digs = 0u;
   while (!mp_iszero(&t)) {
      if (--maxlen < 1u) {
         /* no more room */
         err = MP_BUF;
         goto LBL_ERR;
      }
      if ((err = mp_div_d(&t, (mp_digit)radix, &t, &d)) != MP_OKAY) {
         goto LBL_ERR;
      }
      *str++ = s_mp_radix_map[d];
      ++digs;
   }
   /* reverse the digits of the string.  In this case _s points
    * to the first digit [excluding the sign] of the number
    */
   s_reverse(_s, digs);

   /* append a NULL so the string is properly terminated */
   *str = '\0';
   digs++;

   if (written != NULL) {
      *written = mp_isneg(a) ? (digs + 1u): digs;
   }

LBL_ERR:
   mp_clear(&t);
   return err;
}


/* this is a modified version of s_mp_mul_comba that only produces
 * output digits *above* digs.  See the comments for s_mp_mul_comba
 * to see how it works.
 *
 * This is used in the Barrett reduction since for one of the multiplications
 * only the higher digits were needed.  This essentially halves the work.
 *
 * Based on Algorithm 14.12 on pp.595 of HAC.
 */
mp_err s_mp_mul_high_comba(const mp_int *a, const mp_int *b, mp_int *c, int digs)
{
   int     oldused, pa, ix;
   mp_err   err;
   mp_digit W[MP_WARRAY];
   mp_word  _W;

   /* grow the destination as required */
   pa = a->used + b->used;
   if ((err = mp_grow(c, pa)) != MP_OKAY) {
      return err;
   }

   /* number of output digits to produce */
   pa = a->used + b->used;
   _W = 0;
   for (ix = digs; ix < pa; ix++) {
      int      tx, ty, iy, iz;

      /* get offsets into the two bignums */
      ty = MP_MIN(b->used-1, ix);
      tx = ix - ty;

      /* this is the number of times the loop will iterrate, essentially its
         while (tx++ < a->used && ty-- >= 0) { ... }
       */
      iy = MP_MIN(a->used-tx, ty+1);

      /* execute loop */
      for (iz = 0; iz < iy; iz++) {
         _W += (mp_word)a->dp[tx + iz] * (mp_word)b->dp[ty - iz];
      }

      /* store term */
      W[ix] = (mp_digit)_W & MP_MASK;

      /* make next carry */
      _W = _W >> (mp_word)MP_DIGIT_BIT;
   }

   /* setup dest */
   oldused  = c->used;
   c->used = pa;

   for (ix = digs; ix < pa; ix++) {
      /* now extract the previous digit [below the carry] */
      c->dp[ix] = W[ix];
   }

   /* clear unused digits [that existed in the old copy of c] */
   s_mp_zero_digs(c->dp + c->used, oldused - c->used);

   mp_clamp(c);
   return MP_OKAY;
}

/* multiplies |a| * |b| and does not compute the lower digs digits
 * [meant to get the higher part of the product]
 */
mp_err s_mp_mul_high(const mp_int *a, const mp_int *b, mp_int *c, int digs)
{
   mp_int   t;
   int      pa, pb, ix;
   mp_err   err;

   /* can we use the fast multiplier? */
   if ( MP_HAS(S_MP_MUL_HIGH_COMBA) && ((a->used + b->used + 1) < MP_WARRAY) && (MP_MIN(a->used, b->used) < MP_MAX_COMBA) )
      {
      return s_mp_mul_high_comba(a, b, c, digs);
      }

   if ((err = mp_init_size(&t, a->used + b->used + 1)) != MP_OKAY) {
      return err;
   }
   t.used = a->used + b->used + 1;

   pa = a->used;
   pb = b->used;
   for (ix = 0; ix < pa; ix++) {
      int iy;
      mp_digit u = 0;

      for (iy = digs - ix; iy < pb; iy++) {
         /* calculate the double precision result */
         mp_word r = (mp_word)t.dp[ix + iy] +
                     ((mp_word)a->dp[ix] * (mp_word)b->dp[iy]) +
                     (mp_word)u;

         /* get the lower part */
         t.dp[ix + iy] = (mp_digit)(r & (mp_word)MP_MASK);

         /* carry the carry */
         u       = (mp_digit)(r >> (mp_word)MP_DIGIT_BIT);
      }
      t.dp[ix + pb] = u;
   }
   mp_clamp(&t);
   mp_exch(&t, c);
   mp_clear(&t);
   return MP_OKAY;
}






/* reduce "x" in place modulo "n" using the Diminished Radix algorithm.
 *
 * Based on algorithm from the paper
 *
 * "Generating Efficient Primes for Discrete Log Cryptosystems"
 *                 Chae Hoon Lim, Pil Joong Lee,
 *          POSTECH Information Research Laboratories
 *
 * The modulus must be of a special format [see manual]
 *
 * Has been modified to use algorithm 7.10 from the LTM book instead
 *
 * Input x must be in the range 0 <= x <= (n-1)**2
 */
mp_err mp_dr_reduce(mp_int *x, const mp_int *n, mp_digit k)
{
   mp_err err;

   /* m = digits in modulus */
   int m = n->used;

   /* ensure that "x" has at least 2m digits */
   if ((err = mp_grow(x, m + m)) != MP_OKAY) {
      return err;
   }

   /* top of loop, this is where the code resumes if
    * another reduction pass is required.
    */
   for (;;) {
      int i;
      mp_digit mu = 0;

      /* compute (x mod B**m) + k * [x/B**m] inline and inplace */
      for (i = 0; i < m; i++) {
         mp_word r         = ((mp_word)x->dp[i + m] * (mp_word)k) + x->dp[i] + mu;
         x->dp[i]  = (mp_digit)(r & MP_MASK);
         mu        = (mp_digit)(r >> ((mp_word)MP_DIGIT_BIT));
      }

      /* set final carry */
      x->dp[i] = mu;

      /* zero words above m */
      s_mp_zero_digs(x->dp + m + 1, (x->used - m) - 1);

      /* clamp, sub and return */
      mp_clamp(x);

      /* if x >= n then subtract and reduce again
       * Each successive "recursion" makes the input smaller and smaller.
       */
      if (mp_cmp_mag(x, n) == MP_LT) {
         break;
      }

      if ((err = s_mp_sub(x, n, x)) != MP_OKAY) {
         return err;
      }
   }
   return MP_OKAY;
}
/* reduces a modulo n where n is of the form 2**p - d */
mp_err mp_reduce_2k(mp_int *a, const mp_int *n, mp_digit d)
{
   mp_int q;
   mp_err err;
   int p;

   if ((err = mp_init(&q)) != MP_OKAY) {
      return err;
   }

   p = mp_count_bits(n);
   for (;;) {
      /* q = a/2**p, a = a mod 2**p */
      if ((err = mp_div_2d(a, p, &q, a)) != MP_OKAY) {
         goto LBL_ERR;
      }

      if (d != 1u) {
         /* q = q * d */
         if ((err = mp_mul_d(&q, d, &q)) != MP_OKAY) {
            goto LBL_ERR;
         }
      }

      /* a = a + q */
      if ((err = s_mp_add(a, &q, a)) != MP_OKAY) {
         goto LBL_ERR;
      }

      if (mp_cmp_mag(a, n) == MP_LT) {
         break;
      }
      if ((err = s_mp_sub(a, n, a)) != MP_OKAY) {
         goto LBL_ERR;
      }
   }

LBL_ERR:
   mp_clear(&q);
   return err;
}


/*
 * shifts with subtractions when the result is greater than b.
 *
 * The method is slightly modified to shift B unconditionally upto just under
 * the leading bit of b.  This saves alot of multiple precision shifting.
 */
mp_err mp_montgomery_calc_normalization(mp_int *a, const mp_int *b)
{
   int    x, bits;
   mp_err err;

   /* how many bits of last digit does b use */
   bits = mp_count_bits(b) % MP_DIGIT_BIT;

   if (b->used > 1) {
      if ((err = mp_2expt(a, ((b->used - 1) * MP_DIGIT_BIT) + bits - 1)) != MP_OKAY) {
         return err;
      }
   } else {
      mp_set(a, 1uL);
      bits = 1;
   }

   /* now compute C = A * B mod b */
   for (x = bits - 1; x < (int)MP_DIGIT_BIT; x++) {
      if ((err = mp_mul_2(a, a)) != MP_OKAY) {
         return err;
      }
      if (mp_cmp_mag(a, b) != MP_LT) {
         if ((err = s_mp_sub(a, b, a)) != MP_OKAY) {
            return err;
         }
      }
   }

   return MP_OKAY;
}


/* reduces a modulo n where n is of the form 2**p - d
   This differs from reduce_2k since "d" can be larger
   than a single digit.
*/
mp_err mp_reduce_2k_l(mp_int *a, const mp_int *n, const mp_int *d)
{
   mp_int q;
   mp_err err;
   int    p;

   if ((err = mp_init(&q)) != MP_OKAY) {
      return err;
   }

   p = mp_count_bits(n);

   for (;;) {
      /* q = a/2**p, a = a mod 2**p */
      if ((err = mp_div_2d(a, p, &q, a)) != MP_OKAY) {
         goto LBL_ERR;
      }

      /* q = q * d */
      if ((err = mp_mul(&q, d, &q)) != MP_OKAY) {
         goto LBL_ERR;
      }

      /* a = a + q */
      if ((err = s_mp_add(a, &q, a)) != MP_OKAY) {
         goto LBL_ERR;
      }

      if (mp_cmp_mag(a, n) == MP_LT) {
         break;
      }
      if ((err = s_mp_sub(a, n, a)) != MP_OKAY) {
         goto LBL_ERR;
      }

   }

LBL_ERR:
   mp_clear(&q);
   return err;
}

/* determines the setup value */
mp_err mp_reduce_2k_setup_l(const mp_int *a, mp_int *d)
{
   mp_err err;
   mp_int tmp;

   if ((err = mp_init(&tmp)) != MP_OKAY) {
      return err;
   }

   if ((err = mp_2expt(&tmp, mp_count_bits(a))) != MP_OKAY) {
      goto LBL_ERR;
   }

   if ((err = s_mp_sub(&tmp, a, d)) != MP_OKAY) {
      goto LBL_ERR;
   }

LBL_ERR:
   mp_clear(&tmp);
   return err;
}





/* reduces x mod m, assumes 0 < x < m**2, mu is
 * precomputed via mp_reduce_setup.
 * From HAC pp.604 Algorithm 14.42
 */
mp_err mp_reduce(mp_int *x, const mp_int *m, const mp_int *mu)
{
   mp_int  q;
   mp_err  err;
   int     um = m->used;

   /* q = x */
   if ((err = mp_init_copy(&q, x)) != MP_OKAY) {
      return err;
   }

   /* q1 = x / b**(k-1)  */
   mp_rshd(&q, um - 1);

   /* according to HAC this optimization is ok */
   if ((mp_digit)um > ((mp_digit)1 << (MP_DIGIT_BIT - 1))) {
      if ((err = mp_mul(&q, mu, &q)) != MP_OKAY) {
         goto LBL_ERR;
      }
   } else if (MP_HAS(S_MP_MUL_HIGH)) {
      if ((err = s_mp_mul_high(&q, mu, &q, um)) != MP_OKAY) {
         goto LBL_ERR;
      }
   } else if (MP_HAS(S_MP_MUL_HIGH_COMBA)) {
      if ((err = s_mp_mul_high_comba(&q, mu, &q, um)) != MP_OKAY) {
         goto LBL_ERR;
      }
   } else {
      err = MP_VAL;
      goto LBL_ERR;
   }

   /* q3 = q2 / b**(k+1) */
   mp_rshd(&q, um + 1);

   /* x = x mod b**(k+1), quick (no division) */
   if ((err = mp_mod_2d(x, MP_DIGIT_BIT * (um + 1), x)) != MP_OKAY) {
      goto LBL_ERR;
   }

   /* q = q * m mod b**(k+1), quick (no division) */
   if ((err = s_mp_mul(&q, m, &q, um + 1)) != MP_OKAY) {
      goto LBL_ERR;
   }

   /* x = x - q */
   if ((err = mp_sub(x, &q, x)) != MP_OKAY) {
      goto LBL_ERR;
   }

   /* If x < 0, add b**(k+1) to it */
   if (mp_cmp_d(x, 0uL) == MP_LT) {
      mp_set(&q, 1uL);
      if ((err = mp_lshd(&q, um + 1)) != MP_OKAY) {
         goto LBL_ERR;
      }
      if ((err = mp_add(x, &q, x)) != MP_OKAY) {
         goto LBL_ERR;
      }
   }

   /* Back off if it's too big */
   while (mp_cmp(x, m) != MP_LT) {
      if ((err = s_mp_sub(x, m, x)) != MP_OKAY) {
         goto LBL_ERR;
      }
   }

LBL_ERR:
   mp_clear(&q);

   return err;
}


/* d = a * b (mod c) */
mp_err mp_mulmod(const mp_int *a, const mp_int *b, const mp_int *c, mp_int *d)
{
   mp_err err;
   if ((err = mp_mul(a, b, d)) != MP_OKAY) {
      return err;
   }
   return mp_mod(d, c, d);
}




/* computes xR**-1 == x (mod N) via Montgomery Reduction */
mp_err mp_montgomery_reduce(mp_int *x, const mp_int *n, mp_digit rho)
{
   mp_err err;
   int ix, digs;

   /* can the fast reduction [comba] method be used?
    *
    * Note that unlike in mul you're safely allowed *less*
    * than the available columns [255 per default] since carries
    * are fixed up in the inner loop.
    */
   digs = (n->used * 2) + 1;
   if ((digs < MP_WARRAY) &&
       (x->used <= MP_WARRAY) &&
       (n->used < MP_MAX_COMBA)) {
      return s_mp_montgomery_reduce_comba(x, n, rho);
   }

   /* grow the input as required */
   if ((err = mp_grow(x, digs)) != MP_OKAY) {
      return err;
   }
   x->used = digs;

   for (ix = 0; ix < n->used; ix++) {
      int iy;
      mp_digit u, mu;

      /* mu = ai * rho mod b
       *
       * The value of rho must be precalculated via
       * montgomery_setup() such that
       * it equals -1/n0 mod b this allows the
       * following inner loop to reduce the
       * input one digit at a time
       */
      mu = (mp_digit)(((mp_word)x->dp[ix] * (mp_word)rho) & MP_MASK);

      /* a = a + mu * m * b**i */

      /* Multiply and add in place */
      u = 0;
      for (iy = 0; iy < n->used; iy++) {
         /* compute product and sum */
         mp_word r = ((mp_word)mu * (mp_word)n->dp[iy]) +
                     (mp_word)u + (mp_word)x->dp[ix + iy];

         /* get carry */
         u       = (mp_digit)(r >> (mp_word)MP_DIGIT_BIT);

         /* fix digit */
         x->dp[ix + iy] = (mp_digit)(r & (mp_word)MP_MASK);
      }
      /* At this point the ix'th digit of x should be zero */

      /* propagate carries upwards as required*/
      while (u != 0u) {
         x->dp[ix + iy]   += u;
         u        = x->dp[ix + iy] >> MP_DIGIT_BIT;
         x->dp[ix + iy] &= MP_MASK;
         ++iy;
      }
   }

   /* at this point the n.used'th least
    * significant digits of x are all zero
    * which means we can shift x to the
    * right by n.used digits and the
    * residue is unchanged.
    */

   /* x = x/b**n.used */
   mp_clamp(x);
   mp_rshd(x, n->used);

   /* if x >= n then x = x - n */
   if (mp_cmp_mag(x, n) != MP_LT) {
      return s_mp_sub(x, n, x);
   }

   return MP_OKAY;
}






/* determines the setup value */
void mp_dr_setup(const mp_int *a, mp_digit *d)
{
   /* the casts are required if MP_DIGIT_BIT is one less than
    * the number of bits in a mp_digit [e.g. MP_DIGIT_BIT==31]
    */
   *d = (mp_digit)(((mp_word)1 << (mp_word)MP_DIGIT_BIT) - (mp_word)a->dp[0]);
}







void s_mp_zero_buf(void *mem, size_t size)
{
#ifdef MP_USE_MEMOPS
   memset(mem, 0, size);
#else
   char *m = (char *)mem;
   while (size-- > 0u) {
      *m++ = '\0';
   }
#endif
}




/* computes xR**-1 == x (mod N) via Montgomery Reduction
 *
 * This is an optimized implementation of montgomery_reduce
 * which uses the comba method to quickly calculate the columns of the
 * reduction.
 *
 * Based on Algorithm 14.32 on pp.601 of HAC.
*/
mp_err s_mp_montgomery_reduce_comba(mp_int *x, const mp_int *n, mp_digit rho)
{
   int     ix, oldused;
   mp_err  err;
   mp_word W[MP_WARRAY];

   if (x->used > MP_WARRAY) {
      return MP_VAL;
   }

   /* get old used count */
   oldused = x->used;

   /* grow a as required */
   if ((err = mp_grow(x, n->used + 1)) != MP_OKAY) {
      return err;
   }

   /* first we have to get the digits of the input into
    * an array of double precision words W[...]
    */

   /* copy the digits of a into W[0..a->used-1] */
   for (ix = 0; ix < x->used; ix++) {
      W[ix] = x->dp[ix];
   }

   /* zero the high words of W[a->used..m->used*2] */
   if (ix < ((n->used * 2) + 1)) {
      s_mp_zero_buf(W + x->used, sizeof(mp_word) * (size_t)(((n->used * 2) + 1) - ix));
   }

   /* now we proceed to zero successive digits
    * from the least significant upwards
    */
   for (ix = 0; ix < n->used; ix++) {
      int iy;
      mp_digit mu;

      /* mu = ai * m' mod b
       *
       * We avoid a double precision multiplication (which isn't required)
       * by casting the value down to a mp_digit.  Note this requires
       * that W[ix-1] have  the carry cleared (see after the inner loop)
       */
      mu = ((W[ix] & MP_MASK) * rho) & MP_MASK;

      /* a = a + mu * m * b**i
       *
       * This is computed in place and on the fly.  The multiplication
       * by b**i is handled by offseting which columns the results
       * are added to.
       *
       * Note the comba method normally doesn't handle carries in the
       * inner loop In this case we fix the carry from the previous
       * column since the Montgomery reduction requires digits of the
       * result (so far) [see above] to work.  This is
       * handled by fixing up one carry after the inner loop.  The
       * carry fixups are done in order so after these loops the
       * first m->used words of W[] have the carries fixed
       */
      for (iy = 0; iy < n->used; iy++) {
         W[ix + iy] += (mp_word)mu * (mp_word)n->dp[iy];
      }

      /* now fix carry for next digit, W[ix+1] */
      W[ix + 1] += W[ix] >> (mp_word)MP_DIGIT_BIT;
   }

   /* now we have to propagate the carries and
    * shift the words downward [all those least
    * significant digits we zeroed].
    */

   for (; ix < (n->used * 2); ix++) {
      W[ix + 1] += W[ix] >> (mp_word)MP_DIGIT_BIT;
   }

   /* copy out, A = A/b**n
    *
    * The result is A/b**n but instead of converting from an
    * array of mp_word to mp_digit than calling mp_rshd
    * we just copy them in the right order
    */

   for (ix = 0; ix < (n->used + 1); ix++) {
      x->dp[ix] = W[n->used + ix] & (mp_word)MP_MASK;
   }

   /* set the max used */
   x->used = n->used + 1;

   /* zero oldused digits, if the input a was larger than
    * m->used+1 we'll have to clear the digits
    */
   s_mp_zero_digs(x->dp + x->used, oldused - x->used);

   mp_clamp(x);

   /* if A >= m then A = A - m */
   if (mp_cmp_mag(x, n) != MP_LT) {
      return s_mp_sub(x, n, x);
   }
   return MP_OKAY;
}

/* hac 14.61, pp608 */
mp_err s_mp_invmod(const mp_int *a, const mp_int *b, mp_int *c)
{
   mp_int  x, y, u, v, A, B, C, D;
   mp_err  err;

   /* b cannot be negative */
   if ((b->sign == MP_NEG) || mp_iszero(b)) {
      return MP_VAL;
   }

   /* init temps */
   if ((err = mp_init_multi(&x, &y, &u, &v,
                            &A, &B, &C, &D, NULL)) != MP_OKAY) {
      return err;
   }

   /* x = a, y = b */
   if ((err = mp_mod(a, b, &x)) != MP_OKAY)                       goto LBL_ERR;
   if ((err = mp_copy(b, &y)) != MP_OKAY)                         goto LBL_ERR;

   /* 2. [modified] if x,y are both even then return an error! */
   if (mp_iseven(&x) && mp_iseven(&y)) {
      err = MP_VAL;
      goto LBL_ERR;
   }

   /* 3. u=x, v=y, A=1, B=0, C=0,D=1 */
   if ((err = mp_copy(&x, &u)) != MP_OKAY)                        goto LBL_ERR;
   if ((err = mp_copy(&y, &v)) != MP_OKAY)                        goto LBL_ERR;
   mp_set(&A, 1uL);
   mp_set(&D, 1uL);

   do {
      /* 4.  while u is even do */
      while (mp_iseven(&u)) {
         /* 4.1 u = u/2 */
         if ((err = mp_div_2(&u, &u)) != MP_OKAY)                    goto LBL_ERR;

         /* 4.2 if A or B is odd then */
         if (mp_isodd(&A) || mp_isodd(&B)) {
            /* A = (A+y)/2, B = (B-x)/2 */
            if ((err = mp_add(&A, &y, &A)) != MP_OKAY)               goto LBL_ERR;
            if ((err = mp_sub(&B, &x, &B)) != MP_OKAY)               goto LBL_ERR;
         }
         /* A = A/2, B = B/2 */
         if ((err = mp_div_2(&A, &A)) != MP_OKAY)                    goto LBL_ERR;
         if ((err = mp_div_2(&B, &B)) != MP_OKAY)                    goto LBL_ERR;
      }

      /* 5.  while v is even do */
      while (mp_iseven(&v)) {
         /* 5.1 v = v/2 */
         if ((err = mp_div_2(&v, &v)) != MP_OKAY)                    goto LBL_ERR;

         /* 5.2 if C or D is odd then */
         if (mp_isodd(&C) || mp_isodd(&D)) {
            /* C = (C+y)/2, D = (D-x)/2 */
            if ((err = mp_add(&C, &y, &C)) != MP_OKAY)               goto LBL_ERR;
            if ((err = mp_sub(&D, &x, &D)) != MP_OKAY)               goto LBL_ERR;
         }
         /* C = C/2, D = D/2 */
         if ((err = mp_div_2(&C, &C)) != MP_OKAY)                    goto LBL_ERR;
         if ((err = mp_div_2(&D, &D)) != MP_OKAY)                    goto LBL_ERR;
      }

      /* 6.  if u >= v then */
      if (mp_cmp(&u, &v) != MP_LT) {
         /* u = u - v, A = A - C, B = B - D */
         if ((err = mp_sub(&u, &v, &u)) != MP_OKAY)                  goto LBL_ERR;

         if ((err = mp_sub(&A, &C, &A)) != MP_OKAY)                  goto LBL_ERR;

         if ((err = mp_sub(&B, &D, &B)) != MP_OKAY)                  goto LBL_ERR;
      } else {
         /* v - v - u, C = C - A, D = D - B */
         if ((err = mp_sub(&v, &u, &v)) != MP_OKAY)                  goto LBL_ERR;

         if ((err = mp_sub(&C, &A, &C)) != MP_OKAY)                  goto LBL_ERR;

         if ((err = mp_sub(&D, &B, &D)) != MP_OKAY)                  goto LBL_ERR;
      }

      /* if not zero goto step 4 */
   } while (!mp_iszero(&u));

   /* now a = C, b = D, gcd == g*v */

   /* if v != 1 then there is no inverse */
   if (mp_cmp_d(&v, 1uL) != MP_EQ) {
      err = MP_VAL;
      goto LBL_ERR;
   }

   /* if its too low */
   while (mp_cmp_d(&C, 0uL) == MP_LT) {
      if ((err = mp_add(&C, b, &C)) != MP_OKAY)                   goto LBL_ERR;
   }

   /* too big */
   while (mp_cmp_mag(&C, b) != MP_LT) {
      if ((err = mp_sub(&C, b, &C)) != MP_OKAY)                   goto LBL_ERR;
   }

   /* C is now the inverse */
   mp_exch(&C, c);

LBL_ERR:
   mp_clear_multi(&x, &y, &u, &v, &A, &B, &C, &D, NULL);
   return err;
}


#define MP_IS_2EXPT(x) (((x) != 0u) && (((x) & ((x) - 1u)) == 0u))

/* computes a = 2**b
 *
 * Simple algorithm which zeroes the int, grows it then just sets one bit
 * as required.
 */
mp_err mp_2expt(mp_int *a, int b)
{
   mp_err    err;

   /* zero a as per default */
   mp_zero(a);

   /* grow a to accomodate the single bit */
   if ((err = mp_grow(a, (b / MP_DIGIT_BIT) + 1)) != MP_OKAY) {
      return err;
   }

   /* set the used count of where the bit will go */
   a->used = (b / MP_DIGIT_BIT) + 1;

   /* put the single bit in its place */
   a->dp[b / MP_DIGIT_BIT] = (mp_digit)1 << (mp_digit)(b % MP_DIGIT_BIT);

   return MP_OKAY;
}


/* setups the montgomery reduction stuff */
mp_err mp_montgomery_setup(const mp_int *n, mp_digit *rho)
{
   mp_digit x, b;

   /* fast inversion mod 2**k
    *
    * Based on the fact that
    *
    * XA = 1 (mod 2**n)  =>  (X(2-XA)) A = 1 (mod 2**2n)
    *                    =>  2*X*A - X*X*A*A = 1
    *                    =>  2*(1) - (1)     = 1
    */
   b = n->dp[0];

   if ((b & 1u) == 0u) {
      return MP_VAL;
   }

   x = (((b + 2u) & 4u) << 1) + b; /* here x*a==1 mod 2**4 */
   x *= 2u - (b * x);              /* here x*a==1 mod 2**8 */
   x *= 2u - (b * x);              /* here x*a==1 mod 2**16 */
#if defined(MP_64BIT) || !(defined(MP_16BIT))
   x *= 2u - (b * x);              /* here x*a==1 mod 2**32 */
#endif
#ifdef MP_64BIT
   x *= 2u - (b * x);              /* here x*a==1 mod 2**64 */
#endif

   /* rho = -1/m mod b */
   *rho = (mp_digit)(((mp_word)1 << (mp_word)MP_DIGIT_BIT) - x) & MP_MASK;

   return MP_OKAY;
}


/* chars used in radix conversions */

const uint8_t s_mp_radix_map_reverse[] = {
   0x3e, 0xff, 0xff, 0xff, 0x3f, 0x00, 0x01, 0x02, 0x03, 0x04, /* +,-./01234 */
   0x05, 0x06, 0x07, 0x08, 0x09, 0xff, 0xff, 0xff, 0xff, 0xff, /* 56789:;<=> */
   0xff, 0xff, 0x0a, 0x0b, 0x0c, 0x0d, 0x0e, 0x0f, 0x10, 0x11, /* ?@ABCDEFGH */
   0x12, 0x13, 0x14, 0x15, 0x16, 0x17, 0x18, 0x19, 0x1a, 0x1b, /* IJKLMNOPQR */
   0x1c, 0x1d, 0x1e, 0x1f, 0x20, 0x21, 0x22, 0x23, 0xff, 0xff, /* STUVWXYZ[\ */
   0xff, 0xff, 0xff, 0xff, 0x24, 0x25, 0x26, 0x27, 0x28, 0x29, /* ]^_`abcdef */
   0x2a, 0x2b, 0x2c, 0x2d, 0x2e, 0x2f, 0x30, 0x31, 0x32, 0x33, /* ghijklmnop */
   0x34, 0x35, 0x36, 0x37, 0x38, 0x39, 0x3a, 0x3b, 0x3c, 0x3d  /* qrstuvwxyz */
};



/* divide by three (based on routine from MPI and the GMP manual) */
mp_err s_mp_div_3(const mp_int *a, mp_int *c, mp_digit *d)
{
   mp_int   q;
   mp_word  w;
   mp_digit b;
   mp_err   err;
   int      ix;

   /* b = 2**MP_DIGIT_BIT / 3 */
   b = ((mp_word)1 << (mp_word)MP_DIGIT_BIT) / (mp_word)3;

   if ((err = mp_init_size(&q, a->used)) != MP_OKAY) {
      return err;
   }

   q.used = a->used;
   q.sign = a->sign;
   w = 0;
   for (ix = a->used; ix --> 0;) {
      mp_word t;
      w = (w << (mp_word)MP_DIGIT_BIT) | (mp_word)a->dp[ix];

      if (w >= 3u) {
         /* multiply w by [1/3] */
         t = (w * (mp_word)b) >> (mp_word)MP_DIGIT_BIT;

         /* now subtract 3 * [w/3] from w, to get the remainder */
         w -= t+t+t;

         /* fixup the remainder as required since
          * the optimization is not exact.
          */
         while (w >= 3u) {
            t += 1u;
            w -= 3u;
         }
      } else {
         t = 0;
      }
      q.dp[ix] = (mp_digit)t;
   }

   /* [optional] store the remainder */
   if (d != NULL) {
      *d = (mp_digit)w;
   }

   /* [optional] store the quotient */
   if (c != NULL) {
      mp_clamp(&q);
      mp_exch(&q, c);
   }
   mp_clear(&q);

   return MP_OKAY;
}



/* computes the modular inverse via binary extended euclidean algorithm,
 * that is c = 1/a mod b
 *
 * Based on slow invmod except this is optimized for the case where b is
 * odd as per HAC Note 14.64 on pp. 610
 */
mp_err s_mp_invmod_odd(const mp_int *a, const mp_int *b, mp_int *c)
{
   mp_int  x, y, u, v, B, D;
   mp_sign sign;
   mp_err  err;

   /* 2. [modified] b must be odd   */
   if (mp_iseven(b)) {
      return MP_VAL;
   }

   /* init all our temps */
   if ((err = mp_init_multi(&x, &y, &u, &v, &B, &D, NULL)) != MP_OKAY) {
      return err;
   }

   /* x == modulus, y == value to invert */
   if ((err = mp_copy(b, &x)) != MP_OKAY)                         goto LBL_ERR;

   /* we need y = |a| */
   if ((err = mp_mod(a, b, &y)) != MP_OKAY)                       goto LBL_ERR;

   /* if one of x,y is zero return an error! */
   if (mp_iszero(&x) || mp_iszero(&y)) {
      err = MP_VAL;
      goto LBL_ERR;
   }

   /* 3. u=x, v=y, A=1, B=0, C=0,D=1 */
   if ((err = mp_copy(&x, &u)) != MP_OKAY)                        goto LBL_ERR;
   if ((err = mp_copy(&y, &v)) != MP_OKAY)                        goto LBL_ERR;
   mp_set(&D, 1uL);

   do {
      /* 4.  while u is even do */
      while (mp_iseven(&u)) {
         /* 4.1 u = u/2 */
         if ((err = mp_div_2(&u, &u)) != MP_OKAY)                    goto LBL_ERR;

         /* 4.2 if B is odd then */
         if (mp_isodd(&B)) {
            if ((err = mp_sub(&B, &x, &B)) != MP_OKAY)               goto LBL_ERR;
         }
         /* B = B/2 */
         if ((err = mp_div_2(&B, &B)) != MP_OKAY)                    goto LBL_ERR;
      }

      /* 5.  while v is even do */
      while (mp_iseven(&v)) {
         /* 5.1 v = v/2 */
         if ((err = mp_div_2(&v, &v)) != MP_OKAY)                    goto LBL_ERR;

         /* 5.2 if D is odd then */
         if (mp_isodd(&D)) {
            /* D = (D-x)/2 */
            if ((err = mp_sub(&D, &x, &D)) != MP_OKAY)               goto LBL_ERR;
         }
         /* D = D/2 */
         if ((err = mp_div_2(&D, &D)) != MP_OKAY)                    goto LBL_ERR;
      }

      /* 6.  if u >= v then */
      if (mp_cmp(&u, &v) != MP_LT) {
         /* u = u - v, B = B - D */
         if ((err = mp_sub(&u, &v, &u)) != MP_OKAY)                  goto LBL_ERR;

         if ((err = mp_sub(&B, &D, &B)) != MP_OKAY)                  goto LBL_ERR;
      } else {
         /* v - v - u, D = D - B */
         if ((err = mp_sub(&v, &u, &v)) != MP_OKAY)                  goto LBL_ERR;

         if ((err = mp_sub(&D, &B, &D)) != MP_OKAY)                  goto LBL_ERR;
      }

      /* if not zero goto step 4 */
   } while (!mp_iszero(&u));

   /* now a = C, b = D, gcd == g*v */

   /* if v != 1 then there is no inverse */
   if (mp_cmp_d(&v, 1uL) != MP_EQ) {
      err = MP_VAL;
      goto LBL_ERR;
   }

   /* b is now the inverse */
   sign = (a->sign > 0) ? MP_NEG:MP_ZPOS;
   while (mp_isneg(&D)) {
      if ((err = mp_add(&D, b, &D)) != MP_OKAY)                   goto LBL_ERR;
   }

   /* too big */
   while (mp_cmp_mag(&D, b) != MP_LT) {
      if ((err = mp_sub(&D, b, &D)) != MP_OKAY)                   goto LBL_ERR;
   }

   mp_exch(&D, c);
   c->sign = sign;
   err = MP_OKAY;

LBL_ERR:
   mp_clear_multi(&x, &y, &u, &v, &B, &D, NULL);
   return err;
}



/* pre-calculate the value required for Barrett reduction
 * For a given modulus "b" it calulates the value required in "a"
 */
mp_err mp_reduce_setup(mp_int *a, const mp_int *b)
{
   mp_err err;
   if ((err = mp_2expt(a, b->used * 2 * MP_DIGIT_BIT)) != MP_OKAY) {
      return err;
   }
   return mp_div(a, b, a, NULL);
}


/* determines the setup value */
mp_err mp_reduce_2k_setup(const mp_int *a, mp_digit *d)
{
   mp_err err;
   mp_int tmp;

   if ((err = mp_init(&tmp)) != MP_OKAY) {
      return err;
   }

   if ((err = mp_2expt(&tmp, mp_count_bits(a))) != MP_OKAY) {
      goto LBL_ERR;
   }

   if ((err = s_mp_sub(&tmp, a, &tmp)) != MP_OKAY) {
      goto LBL_ERR;
   }

   *d = tmp.dp[0];

LBL_ERR:
   mp_clear(&tmp);
   return err;
}


#define MP_GET_MAG(name, type)                                                         \
    type name(const mp_int* a)                                                         \
    {                                                                                  \
        int i = MP_MIN(a->used, (int)((MP_SIZEOF_BITS(type) + MP_DIGIT_BIT - 1) / MP_DIGIT_BIT)); \
        type res = 0u;                                                                 \
        while (i --> 0) {                                                              \
            res <<= ((MP_SIZEOF_BITS(type) <= MP_DIGIT_BIT) ? 0 : MP_DIGIT_BIT);       \
            res |= (type)a->dp[i];                                                     \
            if (MP_SIZEOF_BITS(type) <= MP_DIGIT_BIT) { break; }                       \
        }                                                                              \
        return res;                                                                    \
    }

MP_GET_MAG(mp_get_mag_u32, uint32_t)

#define MP_GET_SIGNED(name, mag, type, utype)                 \
    type name(const mp_int* a)                                \
    {                                                         \
        utype res = mag(a);                                   \
        return mp_isneg(a) ? (type)-res : (type)res;          \
    }

MP_GET_SIGNED(mp_get_i32, mp_get_mag_u32, int32_t, uint32_t)
#define mp_get_u32(a) ((uint32_t)mp_get_i32(a))

MP_SET_UNSIGNED(mp_set_u32, uint32_t)


MP_SET_SIGNED(mp_set_i32, mp_set_u32, int32_t, uint32_t)

/* multiplies |a| * |b| and only computes upto digs digits of result
 * HAC pp. 595, Algorithm 14.12  Modified so you can control how
 * many digits of output are created.
 */
mp_err s_mp_mul(const mp_int *a, const mp_int *b, mp_int *c, int digs)
{
   mp_int  t;
   mp_err  err;
   int     pa, ix;

   /* can we use the fast multiplier? */
   if ((digs < MP_WARRAY) &&
       (MP_MIN(a->used, b->used) < MP_MAX_COMBA)) {
      return s_mp_mul_comba(a, b, c, digs);
   }

   if ((err = mp_init_size(&t, digs)) != MP_OKAY) {
      return err;
   }
   t.used = digs;

   /* compute the digits of the product directly */
   pa = a->used;
   for (ix = 0; ix < pa; ix++) {
      int iy, pb;
      mp_digit u = 0;

      /* limit ourselves to making digs digits of output */
      pb = MP_MIN(b->used, digs - ix);

      /* compute the columns of the output and propagate the carry */
      for (iy = 0; iy < pb; iy++) {
         /* compute the column as a mp_word */
         mp_word r = (mp_word)t.dp[ix + iy] +
                     ((mp_word)a->dp[ix] * (mp_word)b->dp[iy]) +
                     (mp_word)u;

         /* the new column is the lower part of the result */
         t.dp[ix + iy] = (mp_digit)(r & (mp_word)MP_MASK);

         /* get the carry word from the result */
         u       = (mp_digit)(r >> (mp_word)MP_DIGIT_BIT);
      }
      /* set carry if it is placed below digs */
      if ((ix + iy) < digs) {
         t.dp[ix + pb] = u;
      }
   }

   mp_clamp(&t);
   mp_exch(&t, c);

   mp_clear(&t);
   return MP_OKAY;
}




/* Fast (comba) multiplier
 *
 * This is the fast column-array [comba] multiplier.  It is
 * designed to compute the columns of the product first
 * then handle the carries afterwards.  This has the effect
 * of making the nested loops that compute the columns very
 * simple and schedulable on super-scalar processors.
 *
 * This has been modified to produce a variable number of
 * digits of output so if say only a half-product is required
 * you don't have to compute the upper half (a feature
 * required for fast Barrett reduction).
 *
 * Based on Algorithm 14.12 on pp.595 of HAC.
 *
 */
mp_err s_mp_mul_comba(const mp_int *a, const mp_int *b, mp_int *c, int digs)
{
   int      oldused, pa, ix;
   mp_err   err;
   mp_digit W[MP_WARRAY];
   mp_word  _W;

   /* grow the destination as required */
   if ((err = mp_grow(c, digs)) != MP_OKAY) {
      return err;
   }

   /* number of output digits to produce */
   pa = MP_MIN(digs, a->used + b->used);

   /* clear the carry */
   _W = 0;
   for (ix = 0; ix < pa; ix++) {
      int tx, ty, iy, iz;

      /* get offsets into the two bignums */
      ty = MP_MIN(b->used-1, ix);
      tx = ix - ty;

      /* this is the number of times the loop will iterrate, essentially
         while (tx++ < a->used && ty-- >= 0) { ... }
       */
      iy = MP_MIN(a->used-tx, ty+1);

      /* execute loop */
      for (iz = 0; iz < iy; ++iz) {
         _W += (mp_word)a->dp[tx + iz] * (mp_word)b->dp[ty - iz];
      }

      /* store term */
      W[ix] = (mp_digit)_W & MP_MASK;

      /* make next carry */
      _W = _W >> (mp_word)MP_DIGIT_BIT;
   }

   /* setup dest */
   oldused  = c->used;
   c->used = pa;

   for (ix = 0; ix < pa; ix++) {
      /* now extract the previous digit [below the carry] */
      c->dp[ix] = W[ix];
   }

   /* clear unused digits [that existed in the old copy of c] */
   s_mp_zero_digs(c->dp + c->used, oldused - c->used);

   mp_clamp(c);
   return MP_OKAY;
}






/* c = |a| * |b| using Karatsuba Multiplication using
 * three half size multiplications
 *
 * Let B represent the radix [e.g. 2**MP_DIGIT_BIT] and
 * let n represent half of the number of digits in
 * the min(a,b)
 *
 * a = a1 * B**n + a0
 * b = b1 * B**n + b0
 *
 * Then, a * b =>
   a1b1 * B**2n + ((a1 + a0)(b1 + b0) - (a0b0 + a1b1)) * B + a0b0
 *
 * Note that a1b1 and a0b0 are used twice and only need to be
 * computed once.  So in total three half size (half # of
 * digit) multiplications are performed, a0b0, a1b1 and
 * (a1+b1)(a0+b0)
 *
 * Note that a multiplication of half the digits requires
 * 1/4th the number of single precision multiplications so in
 * total after one call 25% of the single precision multiplications
 * are saved.  Note also that the call to mp_mul can end up back
 * in this function if the a0, a1, b0, or b1 are above the threshold.
 * This is known as divide-and-conquer and leads to the famous
 * O(N**lg(3)) or O(N**1.584) work which is asymptopically lower than
 * the standard O(N**2) that the baseline/comba methods use.
 * Generally though the overhead of this method doesn't pay off
 * until a certain size (N ~ 80) is reached.
 */
mp_err s_mp_mul_karatsuba(const mp_int *a, const mp_int *b, mp_int *c)
{
   mp_int  x0, x1, y0, y1, t1, x0y0, x1y1;
   int  B;
   mp_err  err;

   /* min # of digits */
   B = MP_MIN(a->used, b->used);

   /* now divide in two */
   B = B >> 1;

   /* init copy all the temps */
   if ((err = mp_init_size(&x0, B)) != MP_OKAY) {
      goto LBL_ERR;
   }
   if ((err = mp_init_size(&x1, a->used - B)) != MP_OKAY) {
      goto X0;
   }
   if ((err = mp_init_size(&y0, B)) != MP_OKAY) {
      goto X1;
   }
   if ((err = mp_init_size(&y1, b->used - B)) != MP_OKAY) {
      goto Y0;
   }

   /* init temps */
   if ((err = mp_init_size(&t1, B * 2)) != MP_OKAY) {
      goto Y1;
   }
   if ((err = mp_init_size(&x0y0, B * 2)) != MP_OKAY) {
      goto T1;
   }
   if ((err = mp_init_size(&x1y1, B * 2)) != MP_OKAY) {
      goto X0Y0;
   }

   /* now shift the digits */
   x0.used = y0.used = B;
   x1.used = a->used - B;
   y1.used = b->used - B;

   /* we copy the digits directly instead of using higher level functions
    * since we also need to shift the digits
    */
   s_mp_copy_digs(x0.dp, a->dp, x0.used);
   s_mp_copy_digs(y0.dp, b->dp, y0.used);
   s_mp_copy_digs(x1.dp, a->dp + B, x1.used);
   s_mp_copy_digs(y1.dp, b->dp + B, y1.used);

   /* only need to clamp the lower words since by definition the
    * upper words x1/y1 must have a known number of digits
    */
   mp_clamp(&x0);
   mp_clamp(&y0);

   /* now calc the products x0y0 and x1y1 */
   /* after this x0 is no longer required, free temp [x0==t2]! */
   if ((err = mp_mul(&x0, &y0, &x0y0)) != MP_OKAY) {
      goto X1Y1;          /* x0y0 = x0*y0 */
   }
   if ((err = mp_mul(&x1, &y1, &x1y1)) != MP_OKAY) {
      goto X1Y1;          /* x1y1 = x1*y1 */
   }

   /* now calc x1+x0 and y1+y0 */
   if ((err = s_mp_add(&x1, &x0, &t1)) != MP_OKAY) {
      goto X1Y1;          /* t1 = x1 - x0 */
   }
   if ((err = s_mp_add(&y1, &y0, &x0)) != MP_OKAY) {
      goto X1Y1;          /* t2 = y1 - y0 */
   }
   if ((err = mp_mul(&t1, &x0, &t1)) != MP_OKAY) {
      goto X1Y1;          /* t1 = (x1 + x0) * (y1 + y0) */
   }

   /* add x0y0 */
   if ((err = mp_add(&x0y0, &x1y1, &x0)) != MP_OKAY) {
      goto X1Y1;          /* t2 = x0y0 + x1y1 */
   }
   if ((err = s_mp_sub(&t1, &x0, &t1)) != MP_OKAY) {
      goto X1Y1;          /* t1 = (x1+x0)*(y1+y0) - (x1y1 + x0y0) */
   }

   /* shift by B */
   if ((err = mp_lshd(&t1, B)) != MP_OKAY) {
      goto X1Y1;          /* t1 = (x0y0 + x1y1 - (x1-x0)*(y1-y0))<<B */
   }
   if ((err = mp_lshd(&x1y1, B * 2)) != MP_OKAY) {
      goto X1Y1;          /* x1y1 = x1y1 << 2*B */
   }

   if ((err = mp_add(&x0y0, &t1, &t1)) != MP_OKAY) {
      goto X1Y1;          /* t1 = x0y0 + t1 */
   }
   if ((err = mp_add(&t1, &x1y1, c)) != MP_OKAY) {
      goto X1Y1;          /* t1 = x0y0 + t1 + x1y1 */
   }

X1Y1:
   mp_clear(&x1y1);
X0Y0:
   mp_clear(&x0y0);
T1:
   mp_clear(&t1);
Y1:
   mp_clear(&y1);
Y0:
   mp_clear(&y0);
X1:
   mp_clear(&x1);
X0:
   mp_clear(&x0);
LBL_ERR:
   return err;
}





/*
   Setup from

     Chung, Jaewook, and M. Anwar Hasan. "Asymmetric squaring formulae."
     18th IEEE Symposium on Computer Arithmetic (ARITH'07). IEEE, 2007.

   The interpolation from above needed one temporary variable more
   than the interpolation here:

     Bodrato, Marco, and Alberto Zanoni. "What about Toom-Cook matrices optimality."
     Centro Vito Volterra Universita di Roma Tor Vergata (2006)
*/

mp_err s_mp_mul_toom(const mp_int *a, const mp_int *b, mp_int *c)
{
   mp_int S1, S2, T1, a0, a1, a2, b0, b1, b2;
   int B;
   mp_err err;

   /* init temps */
   if ((err = mp_init_multi(&S1, &S2, &T1, NULL)) != MP_OKAY) {
      return err;
   }

   /* B */
   B = MP_MIN(a->used, b->used) / 3;

   /** a = a2 * x^2 + a1 * x + a0; */
   if ((err = mp_init_size(&a0, B)) != MP_OKAY)                   goto LBL_ERRa0;
   if ((err = mp_init_size(&a1, B)) != MP_OKAY)                   goto LBL_ERRa1;
   if ((err = mp_init_size(&a2, a->used - 2 * B)) != MP_OKAY)     goto LBL_ERRa2;

   a0.used = a1.used = B;
   a2.used = a->used - 2 * B;
   s_mp_copy_digs(a0.dp, a->dp, a0.used);
   s_mp_copy_digs(a1.dp, a->dp + B, a1.used);
   s_mp_copy_digs(a2.dp, a->dp + 2 * B, a2.used);
   mp_clamp(&a0);
   mp_clamp(&a1);
   mp_clamp(&a2);

   /** b = b2 * x^2 + b1 * x + b0; */
   if ((err = mp_init_size(&b0, B)) != MP_OKAY)                   goto LBL_ERRb0;
   if ((err = mp_init_size(&b1, B)) != MP_OKAY)                   goto LBL_ERRb1;
   if ((err = mp_init_size(&b2, b->used - 2 * B)) != MP_OKAY)     goto LBL_ERRb2;

   b0.used = b1.used = B;
   b2.used = b->used - 2 * B;
   s_mp_copy_digs(b0.dp, b->dp, b0.used);
   s_mp_copy_digs(b1.dp, b->dp + B, b1.used);
   s_mp_copy_digs(b2.dp, b->dp + 2 * B, b2.used);
   mp_clamp(&b0);
   mp_clamp(&b1);
   mp_clamp(&b2);

   /** \\ S1 = (a2+a1+a0) * (b2+b1+b0); */
   /** T1 = a2 + a1; */
   if ((err = mp_add(&a2, &a1, &T1)) != MP_OKAY)                  goto LBL_ERR;

   /** S2 = T1 + a0; */
   if ((err = mp_add(&T1, &a0, &S2)) != MP_OKAY)                  goto LBL_ERR;

   /** c = b2 + b1; */
   if ((err = mp_add(&b2, &b1, c)) != MP_OKAY)                    goto LBL_ERR;

   /** S1 = c + b0; */
   if ((err = mp_add(c, &b0, &S1)) != MP_OKAY)                    goto LBL_ERR;

   /** S1 = S1 * S2; */
   if ((err = mp_mul(&S1, &S2, &S1)) != MP_OKAY)                  goto LBL_ERR;

   /** \\S2 = (4*a2+2*a1+a0) * (4*b2+2*b1+b0); */
   /** T1 = T1 + a2; */
   if ((err = mp_add(&T1, &a2, &T1)) != MP_OKAY)                  goto LBL_ERR;

   /** T1 = T1 << 1; */
   if ((err = mp_mul_2(&T1, &T1)) != MP_OKAY)                     goto LBL_ERR;

   /** T1 = T1 + a0; */
   if ((err = mp_add(&T1, &a0, &T1)) != MP_OKAY)                  goto LBL_ERR;

   /** c = c + b2; */
   if ((err = mp_add(c, &b2, c)) != MP_OKAY)                      goto LBL_ERR;

   /** c = c << 1; */
   if ((err = mp_mul_2(c, c)) != MP_OKAY)                         goto LBL_ERR;

   /** c = c + b0; */
   if ((err = mp_add(c, &b0, c)) != MP_OKAY)                      goto LBL_ERR;

   /** S2 = T1 * c; */
   if ((err = mp_mul(&T1, c, &S2)) != MP_OKAY)                    goto LBL_ERR;

   /** \\S3 = (a2-a1+a0) * (b2-b1+b0); */
   /** a1 = a2 - a1; */
   if ((err = mp_sub(&a2, &a1, &a1)) != MP_OKAY)                  goto LBL_ERR;

   /** a1 = a1 + a0; */
   if ((err = mp_add(&a1, &a0, &a1)) != MP_OKAY)                  goto LBL_ERR;

   /** b1 = b2 - b1; */
   if ((err = mp_sub(&b2, &b1, &b1)) != MP_OKAY)                  goto LBL_ERR;

   /** b1 = b1 + b0; */
   if ((err = mp_add(&b1, &b0, &b1)) != MP_OKAY)                  goto LBL_ERR;

   /** a1 = a1 * b1; */
   if ((err = mp_mul(&a1, &b1, &a1)) != MP_OKAY)                  goto LBL_ERR;

   /** b1 = a2 * b2; */
   if ((err = mp_mul(&a2, &b2, &b1)) != MP_OKAY)                  goto LBL_ERR;

   /** \\S2 = (S2 - S3)/3; */
   /** S2 = S2 - a1; */
   if ((err = mp_sub(&S2, &a1, &S2)) != MP_OKAY)                  goto LBL_ERR;

   /** S2 = S2 / 3; \\ this is an exact division  */
   if ((err = s_mp_div_3(&S2, &S2, NULL)) != MP_OKAY)             goto LBL_ERR;

   /** a1 = S1 - a1; */
   if ((err = mp_sub(&S1, &a1, &a1)) != MP_OKAY)                  goto LBL_ERR;

   /** a1 = a1 >> 1; */
   if ((err = mp_div_2(&a1, &a1)) != MP_OKAY)                     goto LBL_ERR;

   /** a0 = a0 * b0; */
   if ((err = mp_mul(&a0, &b0, &a0)) != MP_OKAY)                  goto LBL_ERR;

   /** S1 = S1 - a0; */
   if ((err = mp_sub(&S1, &a0, &S1)) != MP_OKAY)                  goto LBL_ERR;

   /** S2 = S2 - S1; */
   if ((err = mp_sub(&S2, &S1, &S2)) != MP_OKAY)                  goto LBL_ERR;

   /** S2 = S2 >> 1; */
   if ((err = mp_div_2(&S2, &S2)) != MP_OKAY)                     goto LBL_ERR;

   /** S1 = S1 - a1; */
   if ((err = mp_sub(&S1, &a1, &S1)) != MP_OKAY)                  goto LBL_ERR;

   /** S1 = S1 - b1; */
   if ((err = mp_sub(&S1, &b1, &S1)) != MP_OKAY)                  goto LBL_ERR;

   /** T1 = b1 << 1; */
   if ((err = mp_mul_2(&b1, &T1)) != MP_OKAY)                     goto LBL_ERR;

   /** S2 = S2 - T1; */
   if ((err = mp_sub(&S2, &T1, &S2)) != MP_OKAY)                  goto LBL_ERR;

   /** a1 = a1 - S2; */
   if ((err = mp_sub(&a1, &S2, &a1)) != MP_OKAY)                  goto LBL_ERR;


   /** P = b1*x^4+ S2*x^3+ S1*x^2+ a1*x + a0; */
   if ((err = mp_lshd(&b1, 4 * B)) != MP_OKAY)                    goto LBL_ERR;
   if ((err = mp_lshd(&S2, 3 * B)) != MP_OKAY)                    goto LBL_ERR;
   if ((err = mp_add(&b1, &S2, &b1)) != MP_OKAY)                  goto LBL_ERR;
   if ((err = mp_lshd(&S1, 2 * B)) != MP_OKAY)                    goto LBL_ERR;
   if ((err = mp_add(&b1, &S1, &b1)) != MP_OKAY)                  goto LBL_ERR;
   if ((err = mp_lshd(&a1, 1 * B)) != MP_OKAY)                    goto LBL_ERR;
   if ((err = mp_add(&b1, &a1, &b1)) != MP_OKAY)                  goto LBL_ERR;
   if ((err = mp_add(&b1, &a0, c)) != MP_OKAY)                    goto LBL_ERR;

   /** a * b - P */


LBL_ERR:
   mp_clear(&b2);
LBL_ERRb2:
   mp_clear(&b1);
LBL_ERRb1:
   mp_clear(&b0);
LBL_ERRb0:
   mp_clear(&a2);
LBL_ERRa2:
   mp_clear(&a1);
LBL_ERRa1:
   mp_clear(&a0);
LBL_ERRa0:
   mp_clear_multi(&S1, &S2, &T1, NULL);
   return err;
}



/* single-digit multiplication with the smaller number as the single-digit */
mp_err s_mp_mul_balance(const mp_int *a, const mp_int *b, mp_int *c)
{
   mp_int a0, tmp, r;
   mp_err err;
   int i, j,
       nblocks = MP_MAX(a->used, b->used) / MP_MIN(a->used, b->used),
       bsize = MP_MIN(a->used, b->used);

   if ((err = mp_init_size(&a0, bsize + 2)) != MP_OKAY) {
      return err;
   }
   if ((err = mp_init_multi(&tmp, &r, NULL)) != MP_OKAY) {
      mp_clear(&a0);
      return err;
   }

   /* Make sure that A is the larger one*/
   if (a->used < b->used) {
      MP_EXCH(const mp_int *, a, b);
   }

   for (i = 0, j=0; i < nblocks; i++) {
      /* Cut a slice off of a */
      a0.used = bsize;
      s_mp_copy_digs(a0.dp, a->dp + j, a0.used);
      j += a0.used;
      mp_clamp(&a0);

      /* Multiply with b */
      if ((err = mp_mul(&a0, b, &tmp)) != MP_OKAY) {
         goto LBL_ERR;
      }
      /* Shift tmp to the correct position */
      if ((err = mp_lshd(&tmp, bsize * i)) != MP_OKAY) {
         goto LBL_ERR;
      }
      /* Add to output. No carry needed */
      if ((err = mp_add(&r, &tmp, &r)) != MP_OKAY) {
         goto LBL_ERR;
      }
   }
   /* The left-overs; there are always left-overs */
   if (j < a->used) {
      a0.used = a->used - j;
      s_mp_copy_digs(a0.dp, a->dp + j, a0.used);
      j += a0.used;
      mp_clamp(&a0);

      if ((err = mp_mul(&a0, b, &tmp)) != MP_OKAY) {
         goto LBL_ERR;
      }
      if ((err = mp_lshd(&tmp, bsize * i)) != MP_OKAY) {
         goto LBL_ERR;
      }
      if ((err = mp_add(&r, &tmp, &r)) != MP_OKAY) {
         goto LBL_ERR;
      }
   }

   mp_exch(&r,c);
LBL_ERR:
   mp_clear_multi(&a0, &tmp, &r,NULL);
   return err;
}



/* low level squaring, b = a*a, HAC pp.596-597, Algorithm 14.16 */
mp_err s_mp_sqr(const mp_int *a, mp_int *b)
{
   mp_int   t;
   int      ix, pa;
   mp_err   err;

   pa = a->used;
   if ((err = mp_init_size(&t, (2 * pa) + 1)) != MP_OKAY) {
      return err;
   }

   /* default used is maximum possible size */
   t.used = (2 * pa) + 1;

   for (ix = 0; ix < pa; ix++) {
      mp_digit u;
      int iy;

      /* first calculate the digit at 2*ix */
      /* calculate double precision result */
      mp_word r = (mp_word)t.dp[2*ix] +
                  ((mp_word)a->dp[ix] * (mp_word)a->dp[ix]);

      /* store lower part in result */
      t.dp[ix+ix] = (mp_digit)(r & (mp_word)MP_MASK);

      /* get the carry */
      u           = (mp_digit)(r >> (mp_word)MP_DIGIT_BIT);

      for (iy = ix + 1; iy < pa; iy++) {
         /* first calculate the product */
         r       = (mp_word)a->dp[ix] * (mp_word)a->dp[iy];

         /* now calculate the double precision result, note we use
          * addition instead of *2 since it's easier to optimize
          */
         r       = (mp_word)t.dp[ix + iy] + r + r + (mp_word)u;

         /* store lower part */
         t.dp[ix + iy] = (mp_digit)(r & (mp_word)MP_MASK);

         /* get carry */
         u       = (mp_digit)(r >> (mp_word)MP_DIGIT_BIT);
      }
      /* propagate upwards */
      while (u != 0uL) {
         r       = (mp_word)t.dp[ix + iy] + (mp_word)u;
         t.dp[ix + iy] = (mp_digit)(r & (mp_word)MP_MASK);
         u       = (mp_digit)(r >> (mp_word)MP_DIGIT_BIT);
         ++iy;
      }
   }

   mp_clamp(&t);
   mp_exch(&t, b);
   mp_clear(&t);
   return MP_OKAY;
}



/* the jist of squaring...
 * you do like mult except the offset of the tmpx [one that
 * starts closer to zero] can't equal the offset of tmpy.
 * So basically you set up iy like before then you min it with
 * (ty-tx) so that it never happens.  You double all those
 * you add in the inner loop

After that loop you do the squares and add them in.
*/

mp_err s_mp_sqr_comba(const mp_int *a, mp_int *b)
{
   int       oldused, pa, ix;
   mp_digit  W[MP_WARRAY];
   mp_word   W1;
   mp_err err;

   /* grow the destination as required */
   pa = a->used + a->used;
   if ((err = mp_grow(b, pa)) != MP_OKAY) {
      return err;
   }

   /* number of output digits to produce */
   W1 = 0;
   for (ix = 0; ix < pa; ix++) {
      int      tx, ty, iy, iz;
      mp_word  _W;

      /* clear counter */
      _W = 0;

      /* get offsets into the two bignums */
      ty = MP_MIN(a->used-1, ix);
      tx = ix - ty;

      /* this is the number of times the loop will iterrate, essentially
         while (tx++ < a->used && ty-- >= 0) { ... }
       */
      iy = MP_MIN(a->used-tx, ty+1);

      /* now for squaring tx can never equal ty
       * we halve the distance since they approach at a rate of 2x
       * and we have to round because odd cases need to be executed
       */
      iy = MP_MIN(iy, ((ty-tx)+1)>>1);

      /* execute loop */
      for (iz = 0; iz < iy; iz++) {
         _W += (mp_word)a->dp[tx + iz] * (mp_word)a->dp[ty - iz];
      }

      /* double the inner product and add carry */
      _W = _W + _W + W1;

      /* even columns have the square term in them */
      if (((unsigned)ix & 1u) == 0u) {
         _W += (mp_word)a->dp[ix>>1] * (mp_word)a->dp[ix>>1];
      }

      /* store it */
      W[ix] = (mp_digit)_W & MP_MASK;

      /* make next carry */
      W1 = _W >> (mp_word)MP_DIGIT_BIT;
   }

   /* setup dest */
   oldused  = b->used;
   b->used = a->used+a->used;

   for (ix = 0; ix < pa; ix++) {
      b->dp[ix] = W[ix] & MP_MASK;
   }

   /* clear unused digits [that existed in the old copy of c] */
   s_mp_zero_digs(b->dp + b->used, oldused - b->used);

   mp_clamp(b);
   return MP_OKAY;
}



/* Karatsuba squaring, computes b = a*a using three
 * half size squarings
 *
 * See comments of mul_karatsuba for details.  It
 * is essentially the same algorithm but merely
 * tuned to perform recursive squarings.
 */
mp_err s_mp_sqr_karatsuba(const mp_int *a, mp_int *b)
{
   mp_int  x0, x1, t1, t2, x0x0, x1x1;
   int B;
   mp_err  err;

   /* min # of digits */
   B = a->used;

   /* now divide in two */
   B = B >> 1;

   /* init copy all the temps */
   if ((err = mp_init_size(&x0, B)) != MP_OKAY)
      goto LBL_ERR;
   if ((err = mp_init_size(&x1, a->used - B)) != MP_OKAY)
      goto X0;

   /* init temps */
   if ((err = mp_init_size(&t1, a->used * 2)) != MP_OKAY)
      goto X1;
   if ((err = mp_init_size(&t2, a->used * 2)) != MP_OKAY)
      goto T1;
   if ((err = mp_init_size(&x0x0, B * 2)) != MP_OKAY)
      goto T2;
   if ((err = mp_init_size(&x1x1, (a->used - B) * 2)) != MP_OKAY)
      goto X0X0;

   /* now shift the digits */
   x0.used = B;
   x1.used = a->used - B;
   s_mp_copy_digs(x0.dp, a->dp, x0.used);
   s_mp_copy_digs(x1.dp, a->dp + B, x1.used);
   mp_clamp(&x0);

   /* now calc the products x0*x0 and x1*x1 */
   if ((err = mp_sqr(&x0, &x0x0)) != MP_OKAY)
      goto X1X1;           /* x0x0 = x0*x0 */
   if ((err = mp_sqr(&x1, &x1x1)) != MP_OKAY)
      goto X1X1;           /* x1x1 = x1*x1 */

   /* now calc (x1+x0)**2 */
   if ((err = s_mp_add(&x1, &x0, &t1)) != MP_OKAY)
      goto X1X1;           /* t1 = x1 - x0 */
   if ((err = mp_sqr(&t1, &t1)) != MP_OKAY)
      goto X1X1;           /* t1 = (x1 - x0) * (x1 - x0) */

   /* add x0y0 */
   if ((err = s_mp_add(&x0x0, &x1x1, &t2)) != MP_OKAY)
      goto X1X1;           /* t2 = x0x0 + x1x1 */
   if ((err = s_mp_sub(&t1, &t2, &t1)) != MP_OKAY)
      goto X1X1;           /* t1 = (x1+x0)**2 - (x0x0 + x1x1) */

   /* shift by B */
   if ((err = mp_lshd(&t1, B)) != MP_OKAY)
      goto X1X1;           /* t1 = (x0x0 + x1x1 - (x1-x0)*(x1-x0))<<B */
   if ((err = mp_lshd(&x1x1, B * 2)) != MP_OKAY)
      goto X1X1;           /* x1x1 = x1x1 << 2*B */

   if ((err = mp_add(&x0x0, &t1, &t1)) != MP_OKAY)
      goto X1X1;           /* t1 = x0x0 + t1 */
   if ((err = mp_add(&t1, &x1x1, b)) != MP_OKAY)
      goto X1X1;           /* t1 = x0x0 + t1 + x1x1 */

X1X1:
   mp_clear(&x1x1);
X0X0:
   mp_clear(&x0x0);
T2:
   mp_clear(&t2);
T1:
   mp_clear(&t1);
X1:
   mp_clear(&x1);
X0:
   mp_clear(&x0);
LBL_ERR:
   return err;
}



/* squaring using Toom-Cook 3-way algorithm */

/*
   This file contains code from J. Arndt's book  "Matters Computational"
   and the accompanying FXT-library with permission of the author.
*/

/* squaring using Toom-Cook 3-way algorithm */
/*
   Setup and interpolation from algorithm SQR_3 in

     Chung, Jaewook, and M. Anwar Hasan. "Asymmetric squaring formulae."
     18th IEEE Symposium on Computer Arithmetic (ARITH'07). IEEE, 2007.

*/
mp_err s_mp_sqr_toom(const mp_int *a, mp_int *b)
{
   mp_int S0, a0, a1, a2;
   int B;
   mp_err err;

   /* init temps */
   if ((err = mp_init(&S0)) != MP_OKAY) {
      return err;
   }

   /* B */
   B = a->used / 3;

   /** a = a2 * x^2 + a1 * x + a0; */
   if ((err = mp_init_size(&a0, B)) != MP_OKAY)                   goto LBL_ERRa0;
   if ((err = mp_init_size(&a1, B)) != MP_OKAY)                   goto LBL_ERRa1;
   if ((err = mp_init_size(&a2, a->used - (2 * B))) != MP_OKAY)   goto LBL_ERRa2;

   a0.used = a1.used = B;
   a2.used = a->used - 2 * B;
   s_mp_copy_digs(a0.dp, a->dp, a0.used);
   s_mp_copy_digs(a1.dp, a->dp + B, a1.used);
   s_mp_copy_digs(a2.dp, a->dp + 2 * B, a2.used);
   mp_clamp(&a0);
   mp_clamp(&a1);
   mp_clamp(&a2);

   /** S0 = a0^2;  */
   if ((err = mp_sqr(&a0, &S0)) != MP_OKAY)                       goto LBL_ERR;

   /** \\S1 = (a2 + a1 + a0)^2 */
   /** \\S2 = (a2 - a1 + a0)^2  */
   /** \\S1 = a0 + a2; */
   /** a0 = a0 + a2; */
   if ((err = mp_add(&a0, &a2, &a0)) != MP_OKAY)                  goto LBL_ERR;
   /** \\S2 = S1 - a1; */
   /** b = a0 - a1; */
   if ((err = mp_sub(&a0, &a1, b)) != MP_OKAY)                    goto LBL_ERR;
   /** \\S1 = S1 + a1; */
   /** a0 = a0 + a1; */
   if ((err = mp_add(&a0, &a1, &a0)) != MP_OKAY)                  goto LBL_ERR;
   /** \\S1 = S1^2;  */
   /** a0 = a0^2; */
   if ((err = mp_sqr(&a0, &a0)) != MP_OKAY)                       goto LBL_ERR;
   /** \\S2 = S2^2;  */
   /** b = b^2; */
   if ((err = mp_sqr(b, b)) != MP_OKAY)                           goto LBL_ERR;

   /** \\ S3 = 2 * a1 * a2  */
   /** \\S3 = a1 * a2;  */
   /** a1 = a1 * a2; */
   if ((err = mp_mul(&a1, &a2, &a1)) != MP_OKAY)                  goto LBL_ERR;
   /** \\S3 = S3 << 1;  */
   /** a1 = a1 << 1; */
   if ((err = mp_mul_2(&a1, &a1)) != MP_OKAY)                     goto LBL_ERR;

   /** \\S4 = a2^2;  */
   /** a2 = a2^2; */
   if ((err = mp_sqr(&a2, &a2)) != MP_OKAY)                       goto LBL_ERR;

   /** \\ tmp = (S1 + S2)/2  */
   /** \\tmp = S1 + S2; */
   /** b = a0 + b; */
   if ((err = mp_add(&a0, b, b)) != MP_OKAY)                      goto LBL_ERR;
   /** \\tmp = tmp >> 1; */
   /** b = b >> 1; */
   if ((err = mp_div_2(b, b)) != MP_OKAY)                         goto LBL_ERR;

   /** \\ S1 = S1 - tmp - S3  */
   /** \\S1 = S1 - tmp; */
   /** a0 = a0 - b; */
   if ((err = mp_sub(&a0, b, &a0)) != MP_OKAY)                    goto LBL_ERR;
   /** \\S1 = S1 - S3;  */
   /** a0 = a0 - a1; */
   if ((err = mp_sub(&a0, &a1, &a0)) != MP_OKAY)                  goto LBL_ERR;

   /** \\S2 = tmp - S4 -S0  */
   /** \\S2 = tmp - S4;  */
   /** b = b - a2; */
   if ((err = mp_sub(b, &a2, b)) != MP_OKAY)                      goto LBL_ERR;
   /** \\S2 = S2 - S0;  */
   /** b = b - S0; */
   if ((err = mp_sub(b, &S0, b)) != MP_OKAY)                      goto LBL_ERR;


   /** \\P = S4*x^4 + S3*x^3 + S2*x^2 + S1*x + S0; */
   /** P = a2*x^4 + a1*x^3 + b*x^2 + a0*x + S0; */

   if ((err = mp_lshd(&a2, 4 * B)) != MP_OKAY)                    goto LBL_ERR;
   if ((err = mp_lshd(&a1, 3 * B)) != MP_OKAY)                    goto LBL_ERR;
   if ((err = mp_lshd(b, 2 * B)) != MP_OKAY)                      goto LBL_ERR;
   if ((err = mp_lshd(&a0, 1 * B)) != MP_OKAY)                    goto LBL_ERR;
   if ((err = mp_add(&a2, &a1, &a2)) != MP_OKAY)                  goto LBL_ERR;
   if ((err = mp_add(&a2, b, b)) != MP_OKAY)                      goto LBL_ERR;
   if ((err = mp_add(b, &a0, b)) != MP_OKAY)                      goto LBL_ERR;
   if ((err = mp_add(b, &S0, b)) != MP_OKAY)                      goto LBL_ERR;
   /** a^2 - P  */


LBL_ERR:
   mp_clear(&a2);
LBL_ERRa2:
   mp_clear(&a1);
LBL_ERRa1:
   mp_clear(&a0);
LBL_ERRa0:
   mp_clear(&S0);

   return err;
}


/* low level subtraction (assumes |a| > |b|), HAC pp.595 Algorithm 14.9 */
mp_err s_mp_sub(const mp_int *a, const mp_int *b, mp_int *c)
{
   int oldused = c->used, min = b->used, max = a->used, i;
   mp_digit u;
   mp_err err;

   /* init result */
   if ((err = mp_grow(c, max)) != MP_OKAY) {
      return err;
   }

   c->used = max;

   /* set carry to zero */
   u = 0;
   for (i = 0; i < min; i++) {
      /* T[i] = A[i] - B[i] - U */
      c->dp[i] = (a->dp[i] - b->dp[i]) - u;

      /* U = carry bit of T[i]
       * Note this saves performing an AND operation since
       * if a carry does occur it will propagate all the way to the
       * MSB.  As a result a single shift is enough to get the carry
       */
      u = c->dp[i] >> (MP_SIZEOF_BITS(mp_digit) - 1u);

      /* Clear carry from T[i] */
      c->dp[i] &= MP_MASK;
   }

   /* now copy higher words if any, e.g. if A has more digits than B  */
   for (; i < max; i++) {
      /* T[i] = A[i] - U */
      c->dp[i] = a->dp[i] - u;

      /* U = carry bit of T[i] */
      u = c->dp[i] >> (MP_SIZEOF_BITS(mp_digit) - 1u);

      /* Clear carry from T[i] */
      c->dp[i] &= MP_MASK;
   }

   /* clear digits above used (since we may not have grown result above) */
   s_mp_zero_digs(c->dp + c->used, oldused - c->used);

   mp_clamp(c);
   return MP_OKAY;
}

/* low level addition, based on HAC pp.594, Algorithm 14.7 */
mp_err s_mp_add(const mp_int *a, const mp_int *b, mp_int *c)
{
   int oldused, min, max, i;
   mp_digit u;
   mp_err err;

   /* find sizes, we let |a| <= |b| which means we have to sort
    * them.  "x" will point to the input with the most digits
    */
   if (a->used < b->used) {
      MP_EXCH(const mp_int *, a, b);
   }

   min = b->used;
   max = a->used;

   /* init result */
   if ((err = mp_grow(c, max + 1)) != MP_OKAY) {
      return err;
   }

   /* get old used digit count and set new one */
   oldused = c->used;
   c->used = max + 1;

   /* zero the carry */
   u = 0;
   for (i = 0; i < min; i++) {
      /* Compute the sum at one digit, T[i] = A[i] + B[i] + U */
      c->dp[i] = a->dp[i] + b->dp[i] + u;

      /* U = carry bit of T[i] */
      u = c->dp[i] >> (mp_digit)MP_DIGIT_BIT;

      /* take away carry bit from T[i] */
      c->dp[i] &= MP_MASK;
   }

   /* now copy higher words if any, that is in A+B
    * if A or B has more digits add those in
    */
   if (min != max) {
      for (; i < max; i++) {
         /* T[i] = A[i] + U */
         c->dp[i] = a->dp[i] + u;

         /* U = carry bit of T[i] */
         u = c->dp[i] >> (mp_digit)MP_DIGIT_BIT;

         /* take away carry bit from T[i] */
         c->dp[i] &= MP_MASK;
      }
   }

   /* add carry */
   c->dp[i] = u;

   /* clear digits above oldused */
   s_mp_zero_digs(c->dp + c->used, oldused - c->used);

   mp_clamp(c);
   return MP_OKAY;
}


/* slower bit-bang division... also smaller */
mp_err s_mp_div_small(const mp_int *a, const mp_int *b, mp_int *c, mp_int *d)
{
   mp_int ta, tb, tq, q;
   int n;
   bool neg;
   mp_err err;

   /* init our temps */
   if ((err = mp_init_multi(&ta, &tb, &tq, &q, NULL)) != MP_OKAY) {
      return err;
   }

   mp_set(&tq, 1uL);
   n = mp_count_bits(a) - mp_count_bits(b);
   if ((err = mp_abs(a, &ta)) != MP_OKAY)                         goto LBL_ERR;
   if ((err = mp_abs(b, &tb)) != MP_OKAY)                         goto LBL_ERR;
   if ((err = mp_mul_2d(&tb, n, &tb)) != MP_OKAY)                 goto LBL_ERR;
   if ((err = mp_mul_2d(&tq, n, &tq)) != MP_OKAY)                 goto LBL_ERR;

   while (n-- >= 0) {
      if (mp_cmp(&tb, &ta) != MP_GT) {
         if ((err = mp_sub(&ta, &tb, &ta)) != MP_OKAY)            goto LBL_ERR;
         if ((err = mp_add(&q, &tq, &q)) != MP_OKAY)              goto LBL_ERR;
      }
      if ((err = mp_div_2d(&tb, 1, &tb, NULL)) != MP_OKAY)        goto LBL_ERR;
      if ((err = mp_div_2d(&tq, 1, &tq, NULL)) != MP_OKAY)        goto LBL_ERR;
   }

   /* now q == quotient and ta == remainder */

   neg = (a->sign != b->sign);
   if (c != NULL) {
      mp_exch(c, &q);
      c->sign = ((neg && !mp_iszero(c)) ? MP_NEG : MP_ZPOS);
   }
   if (d != NULL) {
      mp_exch(d, &ta);
      d->sign = (mp_iszero(d) ? MP_ZPOS : a->sign);
   }
LBL_ERR:
   mp_clear_multi(&ta, &tb, &tq, &q, NULL);
   return err;
}



/* integer signed division.
 * c*b + d == a [e.g. a/b, c=quotient, d=remainder]
 * HAC pp.598 Algorithm 14.20
 *
 * Note that the description in HAC is horribly
 * incomplete.  For example, it doesn't consider
 * the case where digits are removed from 'x' in
 * the inner loop.  It also doesn't consider the
 * case that y has fewer than three digits, etc..
 *
 * The overall algorithm is as described as
 * 14.20 from HAC but fixed to treat these cases.
*/
mp_err s_mp_div_school(const mp_int *a, const mp_int *b, mp_int *c, mp_int *d)
{
   mp_int q, x, y, t1, t2;
   int n, t, i, norm;
   bool neg;
   mp_err err;

   if ((err = mp_init_size(&q, a->used + 2)) != MP_OKAY) {
      return err;
   }
   q.used = a->used + 2;

   if ((err = mp_init(&t1)) != MP_OKAY)                           goto LBL_Q;
   if ((err = mp_init(&t2)) != MP_OKAY)                           goto LBL_T1;
   if ((err = mp_init_copy(&x, a)) != MP_OKAY)                    goto LBL_T2;
   if ((err = mp_init_copy(&y, b)) != MP_OKAY)                    goto LBL_X;

   /* fix the sign */
   neg = (a->sign != b->sign);
   x.sign = y.sign = MP_ZPOS;

   /* normalize both x and y, ensure that y >= b/2, [b == 2**MP_DIGIT_BIT] */
   norm = mp_count_bits(&y) % MP_DIGIT_BIT;
   if (norm < (MP_DIGIT_BIT - 1)) {
      norm = (MP_DIGIT_BIT - 1) - norm;
      if ((err = mp_mul_2d(&x, norm, &x)) != MP_OKAY)             goto LBL_Y;
      if ((err = mp_mul_2d(&y, norm, &y)) != MP_OKAY)             goto LBL_Y;
   } else {
      norm = 0;
   }

   /* note hac does 0 based, so if used==5 then its 0,1,2,3,4, e.g. use 4 */
   n = x.used - 1;
   t = y.used - 1;

   /* while (x >= y*b**n-t) do { q[n-t] += 1; x -= y*b**{n-t} } */
   /* y = y*b**{n-t} */
   if ((err = mp_lshd(&y, n - t)) != MP_OKAY)                     goto LBL_Y;

   while (mp_cmp(&x, &y) != MP_LT) {
      ++(q.dp[n - t]);
      if ((err = mp_sub(&x, &y, &x)) != MP_OKAY)                  goto LBL_Y;
   }

   /* reset y by shifting it back down */
   mp_rshd(&y, n - t);

   /* step 3. for i from n down to (t + 1) */
   for (i = n; i >= (t + 1); i--) {
      if (i > x.used) {
         continue;
      }

      /* step 3.1 if xi == yt then set q{i-t-1} to b-1,
       * otherwise set q{i-t-1} to (xi*b + x{i-1})/yt */
      if (x.dp[i] == y.dp[t]) {
         q.dp[(i - t) - 1] = ((mp_digit)1 << (mp_digit)MP_DIGIT_BIT) - (mp_digit)1;
      } else {
         mp_word tmp;
         tmp = (mp_word)x.dp[i] << (mp_word)MP_DIGIT_BIT;
         tmp |= (mp_word)x.dp[i - 1];
         tmp /= (mp_word)y.dp[t];
         if (tmp > (mp_word)MP_MASK) {
            tmp = MP_MASK;
         }
         q.dp[(i - t) - 1] = (mp_digit)(tmp & (mp_word)MP_MASK);
      }

      /* while (q{i-t-1} * (yt * b + y{t-1})) >
               xi * b**2 + xi-1 * b + xi-2

         do q{i-t-1} -= 1;
      */
      q.dp[(i - t) - 1] = (q.dp[(i - t) - 1] + 1uL) & (mp_digit)MP_MASK;
      do {
         q.dp[(i - t) - 1] = (q.dp[(i - t) - 1] - 1uL) & (mp_digit)MP_MASK;

         /* find left hand */
         mp_zero(&t1);
         t1.dp[0] = ((t - 1) < 0) ? 0u : y.dp[t - 1];
         t1.dp[1] = y.dp[t];
         t1.used = 2;
         if ((err = mp_mul_d(&t1, q.dp[(i - t) - 1], &t1)) != MP_OKAY)   goto LBL_Y;

         /* find right hand */
         t2.dp[0] = ((i - 2) < 0) ? 0u : x.dp[i - 2];
         t2.dp[1] = x.dp[i - 1]; /* i >= 1 always holds */
         t2.dp[2] = x.dp[i];
         t2.used = 3;
      } while (mp_cmp_mag(&t1, &t2) == MP_GT);

      /* step 3.3 x = x - q{i-t-1} * y * b**{i-t-1} */
      if ((err = mp_mul_d(&y, q.dp[(i - t) - 1], &t1)) != MP_OKAY)       goto LBL_Y;
      if ((err = mp_lshd(&t1, (i - t) - 1)) != MP_OKAY)                  goto LBL_Y;
      if ((err = mp_sub(&x, &t1, &x)) != MP_OKAY)                        goto LBL_Y;

      /* if x < 0 then { x = x + y*b**{i-t-1}; q{i-t-1} -= 1; } */
      if (mp_isneg(&x)) {
         if ((err = mp_copy(&y, &t1)) != MP_OKAY)                        goto LBL_Y;
         if ((err = mp_lshd(&t1, (i - t) - 1)) != MP_OKAY)               goto LBL_Y;
         if ((err = mp_add(&x, &t1, &x)) != MP_OKAY)                     goto LBL_Y;

         q.dp[(i - t) - 1] = (q.dp[(i - t) - 1] - 1uL) & MP_MASK;
      }
   }

   /* now q is the quotient and x is the remainder
    * [which we have to normalize]
    */

   /* get sign before writing to c */
   x.sign = mp_iszero(&x) ? MP_ZPOS : a->sign;

   if (c != NULL) {
      mp_clamp(&q);
      mp_exch(&q, c);
      c->sign = (neg ? MP_NEG : MP_ZPOS);
   }

   if (d != NULL) {
      if ((err = mp_div_2d(&x, norm, &x, NULL)) != MP_OKAY)       goto LBL_Y;
      mp_exch(&x, d);
   }

LBL_Y:
   mp_clear(&y);
LBL_X:
   mp_clear(&x);
LBL_T2:
   mp_clear(&t2);
LBL_T1:
   mp_clear(&t1);
LBL_Q:
   mp_clear(&q);
   return err;
}



#define mp_decr(a) mp_sub_d((a), 1u, (a))

/*
   Direct implementation of algorithms 1.8 "RecursiveDivRem" and 1.9 "UnbalancedDivision"
   from:

      Brent, Richard P., and Paul Zimmermann. "Modern computer arithmetic"
      Vol. 18. Cambridge University Press, 2010
      Available online at https://arxiv.org/pdf/1004.4710

   pages 19ff. in the above online document.
*/

static mp_err s_recursion(const mp_int *a, const mp_int *b, mp_int *q, mp_int *r)
{
   mp_err err;
   mp_int A1, A2, B1, B0, Q1, Q0, R1, R0, t;
   int m = a->used - b->used, k = m/2;

   if (m < (MP_MUL_KARATSUBA_CUTOFF)) {
      return s_mp_div_school(a, b, q, r);
   }

   if ((err = mp_init_multi(&A1, &A2, &B1, &B0, &Q1, &Q0, &R1, &R0, &t, NULL)) != MP_OKAY) {
      goto LBL_ERR;
   }

   /* B1 = b / beta^k, B0 = b % beta^k*/
   if ((err = mp_div_2d(b, k * MP_DIGIT_BIT, &B1, &B0)) != MP_OKAY)        goto LBL_ERR;

   /* (Q1, R1) =  RecursiveDivRem(A / beta^(2k), B1) */
   if ((err = mp_div_2d(a, 2*k * MP_DIGIT_BIT, &A1, &t)) != MP_OKAY)       goto LBL_ERR;
   if ((err = s_recursion(&A1, &B1, &Q1, &R1)) != MP_OKAY)                 goto LBL_ERR;

   /* A1 = (R1 * beta^(2k)) + (A % beta^(2k)) - (Q1 * B0 * beta^k) */
   if ((err = mp_lshd(&R1, 2*k)) != MP_OKAY)                               goto LBL_ERR;
   if ((err = mp_add(&R1, &t, &A1)) != MP_OKAY)                            goto LBL_ERR;
   if ((err = mp_mul(&Q1, &B0, &t)) != MP_OKAY)                            goto LBL_ERR;
   if ((err = mp_lshd(&t, k)) != MP_OKAY)                                  goto LBL_ERR;
   if ((err = mp_sub(&A1, &t, &A1)) != MP_OKAY)                            goto LBL_ERR;

   /* while A1 < 0 do Q1 = Q1 - 1, A1 = A1 + (beta^k * B) */
   if (mp_cmp_d(&A1, 0uL) == MP_LT) {
      if ((err = mp_mul_2d(b, k * MP_DIGIT_BIT, &t)) != MP_OKAY)           goto LBL_ERR;
      do {
         if ((err = mp_decr(&Q1)) != MP_OKAY)                              goto LBL_ERR;
         if ((err = mp_add(&A1, &t, &A1)) != MP_OKAY)                      goto LBL_ERR;
      } while (mp_cmp_d(&A1, 0uL) == MP_LT);
   }
   /* (Q0, R0) =  RecursiveDivRem(A1 / beta^(k), B1) */
   if ((err = mp_div_2d(&A1, k * MP_DIGIT_BIT, &A1, &t)) != MP_OKAY)       goto LBL_ERR;
   if ((err = s_recursion(&A1, &B1, &Q0, &R0)) != MP_OKAY)                 goto LBL_ERR;

   /* A2 = (R0*beta^k) +  (A1 % beta^k) - (Q0*B0) */
   if ((err = mp_lshd(&R0, k)) != MP_OKAY)                                 goto LBL_ERR;
   if ((err = mp_add(&R0, &t, &A2)) != MP_OKAY)                            goto LBL_ERR;
   if ((err = mp_mul(&Q0, &B0, &t)) != MP_OKAY)                            goto LBL_ERR;
   if ((err = mp_sub(&A2, &t, &A2)) != MP_OKAY)                            goto LBL_ERR;

   /* while A2 < 0 do Q0 = Q0 - 1, A2 = A2 + B */
   while (mp_cmp_d(&A2, 0uL) == MP_LT) {
      if ((err = mp_decr(&Q0)) != MP_OKAY)                                 goto LBL_ERR;
      if ((err = mp_add(&A2, b, &A2)) != MP_OKAY)                          goto LBL_ERR;
   }
   /* return q = (Q1*beta^k) + Q0, r = A2 */
   if ((err = mp_lshd(&Q1, k)) != MP_OKAY)                                 goto LBL_ERR;
   if ((err = mp_add(&Q1, &Q0, q)) != MP_OKAY)                             goto LBL_ERR;

   if ((err = mp_copy(&A2, r)) != MP_OKAY)                                 goto LBL_ERR;

LBL_ERR:
   mp_clear_multi(&A1, &A2, &B1, &B0, &Q1, &Q0, &R1, &R0, &t, NULL);
   return err;
}


mp_err s_mp_div_recursive(const mp_int *a, const mp_int *b, mp_int *q, mp_int *r)
{
   int j, m, n, sigma;
   mp_err err;
   bool neg;
   mp_digit msb_b, msb;
   mp_int A, B, Q, Q1, R, A_div, A_mod;

   if ((err = mp_init_multi(&A, &B, &Q, &Q1, &R, &A_div, &A_mod, NULL)) != MP_OKAY) {
      goto LBL_ERR;
   }

   /* most significant bit of a limb */
   /* assumes  MP_DIGIT_MAX < (sizeof(mp_digit) * CHAR_BIT) */
   msb = (MP_DIGIT_MAX + (mp_digit)(1)) >> 1;
   sigma = 0;
   msb_b = b->dp[b->used - 1];
   while (msb_b < msb) {
      sigma++;
      msb_b <<= 1;
   }
   /* Use that sigma to normalize B */
   if ((err = mp_mul_2d(b, sigma, &B)) != MP_OKAY) {
      goto LBL_ERR;
   }
   if ((err = mp_mul_2d(a, sigma, &A)) != MP_OKAY) {
      goto LBL_ERR;
   }

   /* fix the sign */
   neg = (a->sign != b->sign);
   A.sign = B.sign = MP_ZPOS;

   /*
      If the magnitude of "A" is not more more than twice that of "B" we can work
      on them directly, otherwise we need to work at "A" in chunks
    */
   n = B.used;
   m = A.used - B.used;

   /* Q = 0 */
   mp_zero(&Q);
   while (m > n) {
      /* (q, r) = RecursiveDivRem(A / (beta^(m-n)), B) */
      j = (m - n) * MP_DIGIT_BIT;
      if ((err = mp_div_2d(&A, j, &A_div, &A_mod)) != MP_OKAY)                   goto LBL_ERR;
      if ((err = s_recursion(&A_div, &B, &Q1, &R)) != MP_OKAY)                goto LBL_ERR;
      /* Q = (Q*beta!(n)) + q */
      if ((err = mp_mul_2d(&Q, n * MP_DIGIT_BIT, &Q)) != MP_OKAY)                goto LBL_ERR;
      if ((err = mp_add(&Q, &Q1, &Q)) != MP_OKAY)                                goto LBL_ERR;
      /* A = (r * beta^(m-n)) + (A % beta^(m-n))*/
      if ((err = mp_mul_2d(&R, (m - n) * MP_DIGIT_BIT, &R)) != MP_OKAY)          goto LBL_ERR;
      if ((err = mp_add(&R, &A_mod, &A)) != MP_OKAY)                             goto LBL_ERR;
      /* m = m - n */
      m = m - n;
   }
   /* (q, r) = RecursiveDivRem(A, B) */
   if ((err = s_recursion(&A, &B, &Q1, &R)) != MP_OKAY)                       goto LBL_ERR;
   /* Q = (Q * beta^m) + q, R = r */
   if ((err = mp_mul_2d(&Q, m * MP_DIGIT_BIT, &Q)) != MP_OKAY)                   goto LBL_ERR;
   if ((err = mp_add(&Q, &Q1, &Q)) != MP_OKAY)                                   goto LBL_ERR;

   /* get sign before writing to c */
   R.sign = (mp_iszero(&Q) ? MP_ZPOS : a->sign);

   if (q != NULL) {
      mp_exch(&Q, q);
      q->sign = (neg ? MP_NEG : MP_ZPOS);
   }
   if (r != NULL) {
      /* de-normalize the remainder */
      if ((err = mp_div_2d(&R, sigma, &R, NULL)) != MP_OKAY)                      goto LBL_ERR;
      mp_exch(&R, r);
   }
LBL_ERR:
   mp_clear_multi(&A, &B, &Q, &Q1, &R, &A_div, &A_mod, NULL);
   return err;
}


mp_err mp_div(const mp_int *a, const mp_int *b, mp_int *c, mp_int *d)
{
   mp_err err;

   /* is divisor zero ? */
   if (mp_iszero(b)) {
      return MP_VAL;
   }

   /* if a < b then q = 0, r = a */
   if (mp_cmp_mag(a, b) == MP_LT) {
      if (d != NULL) {
         if ((err = mp_copy(a, d)) != MP_OKAY) {
            return err;
         }
      }
      if (c != NULL) {
         mp_zero(c);
      }
      return MP_OKAY;
   }

   if (MP_HAS(S_MP_DIV_RECURSIVE)
       && (b->used > (2 * MP_MUL_KARATSUBA_CUTOFF))
       && (b->used <= ((a->used/3)*2))) {
      err = s_mp_div_recursive(a, b, c, d);
   } else if (MP_HAS(S_MP_DIV_SCHOOL)) {
      err = s_mp_div_school(a, b, c, d);
   } else if (MP_HAS(S_MP_DIV_SMALL)) {
      err = s_mp_div_small(a, b, c, d);
   } else {
      err = MP_VAL;
   }

   return err;
}


/* compare maginitude of two ints (unsigned) */
mp_ord mp_cmp_mag(const mp_int *a, const mp_int *b)
{
   int n;

   /* compare based on # of non-zero digits */
   if (a->used != b->used) {
      return a->used > b->used ? MP_GT : MP_LT;
   }

   /* compare based on digits  */
   for (n = a->used; n --> 0;) {
      if (a->dp[n] != b->dp[n]) {
         return a->dp[n] > b->dp[n] ? MP_GT : MP_LT;
      }
   }

   return MP_EQ;
}

/* this function is less generic than mp_n_root, simpler and faster */
mp_err mp_sqrt(const mp_int *arg, mp_int *ret)
{
   mp_err err;
   mp_int t1, t2;

   /* must be positive */
   if (mp_isneg(arg)) {
      return MP_VAL;
   }

   /* easy out */
   if (mp_iszero(arg)) {
      mp_zero(ret);
      return MP_OKAY;
   }

   if ((err = mp_init_copy(&t1, arg)) != MP_OKAY) {
      return err;
   }

   if ((err = mp_init(&t2)) != MP_OKAY) {
      goto LBL_ERR2;
   }

   /* First approx. (not very bad for large arg) */
   mp_rshd(&t1, t1.used/2);

   /* t1 > 0  */
   if ((err = mp_div(arg, &t1, &t2, NULL)) != MP_OKAY) {
      goto LBL_ERR1;
   }
   if ((err = mp_add(&t1, &t2, &t1)) != MP_OKAY) {
      goto LBL_ERR1;
   }
   if ((err = mp_div_2(&t1, &t1)) != MP_OKAY) {
      goto LBL_ERR1;
   }
   /* And now t1 > sqrt(arg) */
   do {
      if ((err = mp_div(arg, &t1, &t2, NULL)) != MP_OKAY) {
         goto LBL_ERR1;
      }
      if ((err = mp_add(&t1, &t2, &t1)) != MP_OKAY) {
         goto LBL_ERR1;
      }
      if ((err = mp_div_2(&t1, &t1)) != MP_OKAY) {
         goto LBL_ERR1;
      }
      /* t1 >= sqrt(arg) >= t2 at this point */
   } while (mp_cmp_mag(&t1, &t2) == MP_GT);

   mp_exch(&t1, ret);

LBL_ERR1:
   mp_clear(&t2);
LBL_ERR2:
   mp_clear(&t1);
   return err;
}


#define MP_INIT_INT(name , set, type)                    \
    mp_err name(mp_int * a, type b)                      \
    {                                                    \
        mp_err err;                                      \
        if ((err = mp_init(a)) != MP_OKAY) {             \
            return err;                                  \
        }                                                \
        set(a, b);                                       \
        return MP_OKAY;                                  \
    }
MP_INIT_INT(mp_init_u32, mp_set_u32, uint32_t)

void mp_clear_multi(mp_int *mp, ...)
{
   va_list args;
   va_start(args, mp);
   while (mp != NULL) {
      mp_clear(mp);
      mp = va_arg(args, mp_int *);
   }
   va_end(args);
}


/* high level addition (handles signs) */
mp_err mp_add(const mp_int *a, const mp_int *b, mp_int *c)
{
   /* handle two cases, not four */
   if (a->sign == b->sign) {
      /* both positive or both negative */
      /* add their magnitudes, copy the sign */
      c->sign = a->sign;
      return s_mp_add(a, b, c);
   }

   /* one positive, the other negative */
   /* subtract the one with the greater magnitude from */
   /* the one of the lesser magnitude. The result gets */
   /* the sign of the one with the greater magnitude. */
   if (mp_cmp_mag(a, b) == MP_LT) {
      MP_EXCH(const mp_int *, a, b);
   }

   c->sign = a->sign;
   return s_mp_sub(a, b, c);
}


/* Get bit at position b and return true if the bit is 1, false if it is 0 */
bool s_mp_get_bit(const mp_int *a, int b)
{
   mp_digit bit;
   int limb = b / MP_DIGIT_BIT;

   if (limb < 0 || limb >= a->used) {
      return false;
   }

   bit = (mp_digit)1 << (b % MP_DIGIT_BIT);
   return ((a->dp[limb] & bit) != 0u);
}


/* high level subtraction (handles signs) */
mp_err mp_sub(const mp_int *a, const mp_int *b, mp_int *c)
{
   if (a->sign != b->sign) {
      /* subtract a negative from a positive, OR */
      /* subtract a positive from a negative. */
      /* In either case, ADD their magnitudes, */
      /* and use the sign of the first number. */
      c->sign = a->sign;
      return s_mp_add(a, b, c);
   }

   /* subtract a positive from a positive, OR */
   /* subtract a negative from a negative. */
   /* First, take the difference between their */
   /* magnitudes, then... */
   if (mp_cmp_mag(a, b) == MP_LT) {
      /* The second has a larger magnitude */
      /* The result has the *opposite* sign from */
      /* the first number. */
      c->sign = (!mp_isneg(a) ? MP_NEG : MP_ZPOS);
      MP_EXCH(const mp_int *, a, b);
   } else {
      /* The first has a larger or equal magnitude */
      /* Copy the sign from the first */
      c->sign = a->sign;
   }
   return s_mp_sub(a, b, c);
}


/*
   Kronecker symbol (a|p)
   Straightforward implementation of algorithm 1.4.10 in
   Henri Cohen: "A Course in Computational Algebraic Number Theory"

   @book{cohen2013course,
     title={A course in computational algebraic number theory},
     author={Cohen, Henri},
     volume={138},
     year={2013},
     publisher={Springer Science \& Business Media}
    }
 */
mp_err mp_kronecker(const mp_int *a, const mp_int *p, int *c)
{
   mp_int a1, p1, r;
   mp_err err;
   int v, k;

   static const char table[] = {0, 1, 0, -1, 0, -1, 0, 1};

   if (mp_iszero(p)) {
      if ((a->used == 1) && (a->dp[0] == 1u)) {
         *c = 1;
      } else {
         *c = 0;
      }
      return MP_OKAY;
   }

   if (mp_iseven(a) && mp_iseven(p)) {
      *c = 0;
      return MP_OKAY;
   }

   if ((err = mp_init_copy(&a1, a)) != MP_OKAY) {
      return err;
   }
   if ((err = mp_init_copy(&p1, p)) != MP_OKAY) {
      goto LBL_KRON_0;
   }

   v = mp_cnt_lsb(&p1);
   if ((err = mp_div_2d(&p1, v, &p1, NULL)) != MP_OKAY) {
      goto LBL_KRON_1;
   }

   if ((v & 1) == 0) {
      k = 1;
   } else {
      k = table[a->dp[0] & 7u];
   }

   if (mp_isneg(&p1)) {
      p1.sign = MP_ZPOS;
      if (mp_isneg(&a1)) {
         k = -k;
      }
   }

   if ((err = mp_init(&r)) != MP_OKAY) {
      goto LBL_KRON_1;
   }

   for (;;) {
      if (mp_iszero(&a1)) {
         if (mp_cmp_d(&p1, 1uL) == MP_EQ) {
            *c = k;
            goto LBL_KRON;
         } else {
            *c = 0;
            goto LBL_KRON;
         }
      }

      v = mp_cnt_lsb(&a1);
      if ((err = mp_div_2d(&a1, v, &a1, NULL)) != MP_OKAY) {
         goto LBL_KRON;
      }

      if ((v & 1) == 1) {
         k = k * table[p1.dp[0] & 7u];
      }

      if (mp_isneg(&a1)) {
         /*
          * Compute k = (-1)^((a1)*(p1-1)/4) * k
          * a1.dp[0] + 1 cannot overflow because the MSB
          * of the type mp_digit is not set by definition
          */
         if (((a1.dp[0] + 1u) & p1.dp[0] & 2u) != 0u) {
            k = -k;
         }
      } else {
         /* compute k = (-1)^((a1-1)*(p1-1)/4) * k */
         if ((a1.dp[0] & p1.dp[0] & 2u) != 0u) {
            k = -k;
         }
      }

      if ((err = mp_copy(&a1, &r)) != MP_OKAY) {
         goto LBL_KRON;
      }
      r.sign = MP_ZPOS;
      if ((err = mp_mod(&p1, &r, &a1)) != MP_OKAY) {
         goto LBL_KRON;
      }
      if ((err = mp_copy(&r, &p1)) != MP_OKAY) {
         goto LBL_KRON;
      }
   }

LBL_KRON:
   mp_clear(&r);
LBL_KRON_1:
   mp_clear(&p1);
LBL_KRON_0:
   mp_clear(&a1);

   return err;
}

/* Greatest Common Divisor using the binary method */
mp_err mp_gcd(const mp_int *a, const mp_int *b, mp_int *c)
{
   mp_int  u, v;
   int     k, u_lsb, v_lsb;
   mp_err err;

   /* either zero than gcd is the largest */
   if (mp_iszero(a)) {
      return mp_abs(b, c);
   }
   if (mp_iszero(b)) {
      return mp_abs(a, c);
   }

   /* get copies of a and b we can modify */
   if ((err = mp_init_copy(&u, a)) != MP_OKAY) {
      return err;
   }

   if ((err = mp_init_copy(&v, b)) != MP_OKAY) {
      goto LBL_U;
   }

   /* must be positive for the remainder of the algorithm */
   u.sign = v.sign = MP_ZPOS;

   /* B1.  Find the common power of two for u and v */
   u_lsb = mp_cnt_lsb(&u);
   v_lsb = mp_cnt_lsb(&v);
   k     = MP_MIN(u_lsb, v_lsb);

   if (k > 0) {
      /* divide the power of two out */
      if ((err = mp_div_2d(&u, k, &u, NULL)) != MP_OKAY) {
         goto LBL_V;
      }

      if ((err = mp_div_2d(&v, k, &v, NULL)) != MP_OKAY) {
         goto LBL_V;
      }
   }

   /* divide any remaining factors of two out */
   if (u_lsb != k) {
      if ((err = mp_div_2d(&u, u_lsb - k, &u, NULL)) != MP_OKAY) {
         goto LBL_V;
      }
   }

   if (v_lsb != k) {
      if ((err = mp_div_2d(&v, v_lsb - k, &v, NULL)) != MP_OKAY) {
         goto LBL_V;
      }
   }

   while (!mp_iszero(&v)) {
      /* make sure v is the largest */
      if (mp_cmp_mag(&u, &v) == MP_GT) {
         /* swap u and v to make sure v is >= u */
         mp_exch(&u, &v);
      }

      /* subtract smallest from largest */
      if ((err = s_mp_sub(&v, &u, &v)) != MP_OKAY) {
         goto LBL_V;
      }

      /* Divide out all factors of two */
      if ((err = mp_div_2d(&v, mp_cnt_lsb(&v), &v, NULL)) != MP_OKAY) {
         goto LBL_V;
      }
   }

   /* multiply by 2**k which we divided out at the beginning */
   if ((err = mp_mul_2d(&u, k, c)) != MP_OKAY) {
      goto LBL_V;
   }
   c->sign = MP_ZPOS;
   err = MP_OKAY;
LBL_V:
   mp_clear(&u);
LBL_U:
   mp_clear(&v);
   return err;
}


mp_err mp_init_multi(mp_int *mp, ...)
{
   mp_err err = MP_OKAY;
   int n = 0;                 /* Number of ok inits */
   mp_int *cur_arg = mp;
   va_list args;

   va_start(args, mp);        /* init args to next argument from caller */
   while (cur_arg != NULL) {
      err = mp_init(cur_arg);
      if (err != MP_OKAY) {
         /* Oops - error! Back-track and mp_clear what we already
            succeeded in init-ing, then return error.
         */
         va_list clean_args;

         /* now start cleaning up */
         cur_arg = mp;
         va_start(clean_args, mp);
         while (n-- != 0) {
            mp_clear(cur_arg);
            cur_arg = va_arg(clean_args, mp_int *);
         }
         va_end(clean_args);
         break;
      }
      n++;
      cur_arg = va_arg(args, mp_int *);
   }
   va_end(args);
   return err;
}


mp_err mp_mul(const mp_int *a, const mp_int *b, mp_int *c)
{
   mp_err err;
   int min = MP_MIN(a->used, b->used),
       max = MP_MAX(a->used, b->used),
       digs = a->used + b->used + 1;
   bool neg = (a->sign != b->sign);

   if ((a == b) &&
       MP_HAS(S_MP_SQR_TOOM) && /* use Toom-Cook? */
       (a->used >= MP_SQR_TOOM_CUTOFF)) {
      err = s_mp_sqr_toom(a, c);
   } else if ((a == b) &&
              MP_HAS(S_MP_SQR_KARATSUBA) &&  /* Karatsuba? */
              (a->used >= MP_SQR_KARATSUBA_CUTOFF)) {
      err = s_mp_sqr_karatsuba(a, c);
   } else if ((a == b) &&
              MP_HAS(S_MP_SQR_COMBA) && /* can we use the fast comba multiplier? */
              (((a->used * 2) + 1) < MP_WARRAY) &&
              (a->used < (MP_MAX_COMBA / 2))) {
      err = s_mp_sqr_comba(a, c);
   } else if ((a == b) &&
              MP_HAS(S_MP_SQR)) {
      err = s_mp_sqr(a, c);
   } else if (MP_HAS(S_MP_MUL_BALANCE) &&
              /* Check sizes. The smaller one needs to be larger than the Karatsuba cut-off.
               * The bigger one needs to be at least about one MP_MUL_KARATSUBA_CUTOFF bigger
               * to make some sense, but it depends on architecture, OS, position of the
               * stars... so YMMV.
               * Using it to cut the input into slices small enough for s_mp_mul_comba
               * was actually slower on the author's machine, but YMMV.
               */
              (min >= MP_MUL_KARATSUBA_CUTOFF) &&
              ((max / 2) >= MP_MUL_KARATSUBA_CUTOFF) &&
              /* Not much effect was observed below a ratio of 1:2, but again: YMMV. */
              (max >= (2 * min))) {
      err = s_mp_mul_balance(a,b,c);
   } else if (MP_HAS(S_MP_MUL_TOOM) &&
              (min >= MP_MUL_TOOM_CUTOFF)) {
      err = s_mp_mul_toom(a, b, c);
   } else if (MP_HAS(S_MP_MUL_KARATSUBA) &&
              (min >= MP_MUL_KARATSUBA_CUTOFF)) {
      err = s_mp_mul_karatsuba(a, b, c);
   } else if (MP_HAS(S_MP_MUL_COMBA) &&
              /* can we use the fast multiplier?
               *
               * The fast multiplier can be used if the output will
               * have less than MP_WARRAY digits and the number of
               * digits won't affect carry propagation
               */
              (digs < MP_WARRAY) &&
              (min <= MP_MAX_COMBA)) {
      err = s_mp_mul_comba(a, b, c, digs);
   } else if (MP_HAS(S_MP_MUL)) {
      err = s_mp_mul(a, b, c, digs);
   } else {
      err = MP_VAL;
   }
   c->sign = ((c->used > 0) && neg) ? MP_NEG : MP_ZPOS;
   return err;
}


/* multiply by a digit */
mp_err mp_mul_d(const mp_int *a, mp_digit b, mp_int *c)
{
   mp_digit u;
   mp_err   err;
   int   ix, oldused;

   if (b == 1u) {
      return mp_copy(a, c);
   }

   /* power of two ? */
   if (MP_HAS(MP_MUL_2) && (b == 2u)) {
      return mp_mul_2(a, c);
   }
   if (MP_HAS(MP_MUL_2D) && MP_IS_2EXPT(b)) {
      ix = 1;
      while ((ix < MP_DIGIT_BIT) && (b != (((mp_digit)1)<<ix))) {
         ix++;
      }
      return mp_mul_2d(a, ix, c);
   }

   /* make sure c is big enough to hold a*b */
   if ((err = mp_grow(c, a->used + 1)) != MP_OKAY) {
      return err;
   }

   /* get the original destinations used count */
   oldused = c->used;

   /* set the sign */
   c->sign = a->sign;

   /* zero carry */
   u = 0;

   /* compute columns */
   for (ix = 0; ix < a->used; ix++) {
      /* compute product and carry sum for this term */
      mp_word r       = (mp_word)u + ((mp_word)a->dp[ix] * (mp_word)b);

      /* mask off higher bits to get a single digit */
      c->dp[ix] = (mp_digit)(r & (mp_word)MP_MASK);

      /* send carry into next iteration */
      u       = (mp_digit)(r >> (mp_word)MP_DIGIT_BIT);
   }

   /* store final carry [if any] and increment ix offset  */
   c->dp[ix] = u;

   /* set used count */
   c->used = a->used + 1;

   /* now zero digits above the top */
   s_mp_zero_digs(c->dp + c->used, oldused - c->used);

   mp_clamp(c);

   return MP_OKAY;
}


/* shift right a certain amount of digits */
void mp_rshd(mp_int *a, int b)
{
   int x;

   /* if b <= 0 then ignore it */
   if (b <= 0) {
      return;
   }

   /* if b > used then simply zero it and return */
   if (a->used <= b) {
      mp_zero(a);
      return;
   }

   /* shift the digits down.
    * this is implemented as a sliding window where
    * the window is b-digits long and digits from
    * the top of the window are copied to the bottom
    *
    * e.g.

    b-2 | b-1 | b0 | b1 | b2 | ... | bb |   ---->
                /\                   |      ---->
                 \-------------------/      ---->
    */
   for (x = 0; x < (a->used - b); x++) {
      a->dp[x] = a->dp[x + b];
   }

   /* zero the top digits */
   s_mp_zero_digs(a->dp + a->used - b, b);

   /* remove excess digits */
   a->used -= b;
}


/* calc a value mod 2**b */
mp_err mp_mod_2d(const mp_int *a, int b, mp_int *c)
{
   int x;
   mp_err err;

   if (b < 0) {
      return MP_VAL;
   }

   if (b == 0) {
      mp_zero(c);
      return MP_OKAY;
   }

   /* if the modulus is larger than the value than return */
   if (b >= (a->used * MP_DIGIT_BIT)) {
      return mp_copy(a, c);
   }

   if ((err = mp_copy(a, c)) != MP_OKAY) {
      return err;
   }

   /* zero digits above the last digit of the modulus */
   x = (b / MP_DIGIT_BIT) + (((b % MP_DIGIT_BIT) == 0) ? 0 : 1);
   s_mp_zero_digs(c->dp + x, c->used - x);

   /* clear the digit that is not completely outside/inside the modulus */
   c->dp[b / MP_DIGIT_BIT] &=
      ((mp_digit)1 << (mp_digit)(b % MP_DIGIT_BIT)) - (mp_digit)1;
   mp_clamp(c);
   return MP_OKAY;
}


/* swap the elements of two integers, for cases where you can't simply swap the
 * mp_int pointers around
 */
void mp_exch(mp_int *a, mp_int *b)
{
   MP_EXCH(mp_int, *a, *b);
}


/* init an mp_init for a given size */
mp_err mp_init_size(mp_int *a, int size)
{
   size = MP_MAX(MP_MIN_DIGIT_COUNT, size);

   if (size > MP_MAX_DIGIT_COUNT) {
      return MP_OVF;
   }

   /* alloc mem */
   a->dp = (mp_digit *) MP_CALLOC((size_t)size, sizeof(mp_digit));
   if (a->dp == NULL) {
      return MP_MEM;
   }

   /* set the members */
   a->used  = 0;
   a->alloc = size;
   a->sign  = MP_ZPOS;

   return MP_OKAY;
}


/* determines if mp_reduce_2k can be used */
bool mp_reduce_is_2k(const mp_int *a)
{
   if (mp_iszero(a)) {
      return false;
   } else if (a->used == 1) {
      return true;
   } else if (a->used > 1) {
      int ix, iy = mp_count_bits(a), iw = 1;
      mp_digit iz = 1;

      /* Test every bit from the second digit up, must be 1 */
      for (ix = MP_DIGIT_BIT; ix < iy; ix++) {
         if ((a->dp[iw] & iz) == 0u) {
            return false;
         }
         iz <<= 1;
         if (iz > MP_DIGIT_MAX) {
            ++iw;
            iz = 1;
         }
      }
      return true;
   } else {
      return true;
   }
}



mp_err s_mp_exptmod_fast(const mp_int *G, const mp_int *X, const mp_int *P, mp_int *Y, int redmode)
{
   mp_int  M[TAB_SIZE], res;
   mp_digit buf, mp;
   int     bitbuf, bitcpy, bitcnt, mode, digidx, x, y, winsize;
   mp_err   err;

   /* use a pointer to the reduction algorithm.  This allows us to use
    * one of many reduction algorithms without modding the guts of
    * the code with if statements everywhere.
    */
   mp_err(*redux)(mp_int *x, const mp_int *n, mp_digit rho);

   /* find window size */
   x = mp_count_bits(X);
   if (x <= 7) {
      winsize = 2;
   } else if (x <= 36) {
      winsize = 3;
   } else if (x <= 140) {
      winsize = 4;
   } else if (x <= 450) {
      winsize = 5;
   } else if (x <= 1303) {
      winsize = 6;
   } else if (x <= 3529) {
      winsize = 7;
   } else {
      winsize = 8;
   }

   winsize = MAX_WINSIZE ? MP_MIN(MAX_WINSIZE, winsize) : winsize;

   /* init M array */
   /* init first cell */
   if ((err = mp_init_size(&M[1], P->alloc)) != MP_OKAY) {
      return err;
   }

   /* now init the second half of the array */
   for (x = 1<<(winsize-1); x < (1 << winsize); x++) {
      if ((err = mp_init_size(&M[x], P->alloc)) != MP_OKAY) {
         for (y = 1<<(winsize-1); y < x; y++) {
            mp_clear(&M[y]);
         }
         mp_clear(&M[1]);
         return err;
      }
   }

   /* determine and setup reduction code */
   if (redmode == 0) {
      if (MP_HAS(MP_MONTGOMERY_SETUP)) {
         /* now setup montgomery  */
         if ((err = mp_montgomery_setup(P, &mp)) != MP_OKAY)      goto LBL_M;
      } else {
         err = MP_VAL;
         goto LBL_M;
      }

      /* automatically pick the comba one if available (saves quite a few calls/ifs) */
      if (MP_HAS(S_MP_MONTGOMERY_REDUCE_COMBA) &&
          (((P->used * 2) + 1) < MP_WARRAY) &&
          (P->used < MP_MAX_COMBA)) {
         redux = s_mp_montgomery_reduce_comba;
      } else if (MP_HAS(MP_MONTGOMERY_REDUCE)) {
         /* use slower baseline Montgomery method */
         redux = mp_montgomery_reduce;
      } else {
         err = MP_VAL;
         goto LBL_M;
      }
   } else if (redmode == 1) {
      if (MP_HAS(MP_DR_SETUP) && MP_HAS(MP_DR_REDUCE)) {
         /* setup DR reduction for moduli of the form B**k - b */
         mp_dr_setup(P, &mp);
         redux = mp_dr_reduce;
      } else {
         err = MP_VAL;
         goto LBL_M;
      }
   } else if (MP_HAS(MP_REDUCE_2K_SETUP) && MP_HAS(MP_REDUCE_2K)) {
      /* setup DR reduction for moduli of the form 2**k - b */
      if ((err = mp_reduce_2k_setup(P, &mp)) != MP_OKAY)          goto LBL_M;
      redux = mp_reduce_2k;
   } else {
      err = MP_VAL;
      goto LBL_M;
   }

   /* setup result */
   if ((err = mp_init_size(&res, P->alloc)) != MP_OKAY)           goto LBL_M;

   /* create M table
    *

    *
    * The first half of the table is not computed though accept for M[0] and M[1]
    */

   if (redmode == 0) {
      if (MP_HAS(MP_MONTGOMERY_CALC_NORMALIZATION)) {
         /* now we need R mod m */
         if ((err = mp_montgomery_calc_normalization(&res, P)) != MP_OKAY) goto LBL_RES;

         /* now set M[1] to G * R mod m */
         if ((err = mp_mulmod(G, &res, P, &M[1])) != MP_OKAY)     goto LBL_RES;
      } else {
         err = MP_VAL;
         goto LBL_RES;
      }
   } else {
      mp_set(&res, 1uL);
      if ((err = mp_mod(G, P, &M[1])) != MP_OKAY)                 goto LBL_RES;
   }

   /* compute the value at M[1<<(winsize-1)] by squaring M[1] (winsize-1) times */
   if ((err = mp_copy(&M[1], &M[(size_t)1 << (winsize - 1)])) != MP_OKAY) goto LBL_RES;

   for (x = 0; x < (winsize - 1); x++) {
      if ((err = mp_sqr(&M[(size_t)1 << (winsize - 1)], &M[(size_t)1 << (winsize - 1)])) != MP_OKAY) goto LBL_RES;
      if ((err = redux(&M[(size_t)1 << (winsize - 1)], P, mp)) != MP_OKAY) goto LBL_RES;
   }

   /* create upper table */
   for (x = (1 << (winsize - 1)) + 1; x < (1 << winsize); x++) {
      if ((err = mp_mul(&M[x - 1], &M[1], &M[x])) != MP_OKAY)     goto LBL_RES;
      if ((err = redux(&M[x], P, mp)) != MP_OKAY)                 goto LBL_RES;
   }

   /* set initial mode and bit cnt */
   mode   = 0;
   bitcnt = 1;
   buf    = 0;
   digidx = X->used - 1;
   bitcpy = 0;
   bitbuf = 0;

   for (;;) {
      /* grab next digit as required */
      if (--bitcnt == 0) {
         /* if digidx == -1 we are out of digits so break */
         if (digidx == -1) {
            break;
         }
         /* read next digit and reset bitcnt */
         buf    = X->dp[digidx--];
         bitcnt = (int)MP_DIGIT_BIT;
      }

      /* grab the next msb from the exponent */
      y     = (mp_digit)(buf >> (MP_DIGIT_BIT - 1)) & 1uL;
      buf <<= (mp_digit)1;

      /* if the bit is zero and mode == 0 then we ignore it
       * These represent the leading zero bits before the first 1 bit
       * in the exponent.  Technically this opt is not required but it
       * does lower the # of trivial squaring/reductions used
       */
      if ((mode == 0) && (y == 0)) {
         continue;
      }

      /* if the bit is zero and mode == 1 then we square */
      if ((mode == 1) && (y == 0)) {
         if ((err = mp_sqr(&res, &res)) != MP_OKAY)               goto LBL_RES;
         if ((err = redux(&res, P, mp)) != MP_OKAY)               goto LBL_RES;
         continue;
      }

      /* else we add it to the window */
      bitbuf |= (y << (winsize - ++bitcpy));
      mode    = 2;

      if (bitcpy == winsize) {
         /* ok window is filled so square as required and multiply  */
         /* square first */
         for (x = 0; x < winsize; x++) {
            if ((err = mp_sqr(&res, &res)) != MP_OKAY)            goto LBL_RES;
            if ((err = redux(&res, P, mp)) != MP_OKAY)            goto LBL_RES;
         }

         /* then multiply */
         if ((err = mp_mul(&res, &M[bitbuf], &res)) != MP_OKAY)   goto LBL_RES;
         if ((err = redux(&res, P, mp)) != MP_OKAY)               goto LBL_RES;

         /* empty window and reset */
         bitcpy = 0;
         bitbuf = 0;
         mode   = 1;
      }
   }

   /* if bits remain then square/multiply */
   if ((mode == 2) && (bitcpy > 0)) {
      /* square then multiply if the bit is set */
      for (x = 0; x < bitcpy; x++) {
         if ((err = mp_sqr(&res, &res)) != MP_OKAY)               goto LBL_RES;
         if ((err = redux(&res, P, mp)) != MP_OKAY)               goto LBL_RES;

         /* get next bit of the window */
         bitbuf <<= 1;
         if ((bitbuf & (1 << winsize)) != 0) {
            /* then multiply */
            if ((err = mp_mul(&res, &M[1], &res)) != MP_OKAY)     goto LBL_RES;
            if ((err = redux(&res, P, mp)) != MP_OKAY)            goto LBL_RES;
         }
      }
   }

   if (redmode == 0) {
      /* fixup result if Montgomery reduction is used
       * recall that any value in a Montgomery system is
       * actually multiplied by R mod n.  So we have
       * to reduce one more time to cancel out the factor
       * of R.
       */
      if ((err = redux(&res, P, mp)) != MP_OKAY)                  goto LBL_RES;
   }

   /* swap res with Y */
   mp_exch(&res, Y);
   err = MP_OKAY;
LBL_RES:
   mp_clear(&res);
LBL_M:
   mp_clear(&M[1]);
   for (x = 1<<(winsize-1); x < (1 << winsize); x++) {
      mp_clear(&M[x]);
   }
   return err;
}



/* determines if a number is a valid DR modulus */
bool mp_dr_is_modulus(const mp_int *a)
{
   int ix;

   /* must be at least two digits */
   if (a->used < 2) {
      return false;
   }

   /* must be of the form b**k - a [a <= b] so all
    * but the first digit must be equal to -1 (mod b).
    */
   for (ix = 1; ix < a->used; ix++) {
      if (a->dp[ix] != MP_MASK) {
         return false;
      }
   }
   return true;
}



mp_err s_mp_exptmod(const mp_int *G, const mp_int *X, const mp_int *P, mp_int *Y, int redmode)
{
   mp_int  M[TAB_SIZE], res, mu;
   mp_digit buf;
   mp_err   err;
   int      bitbuf, bitcpy, bitcnt, mode, digidx, x, y, winsize;
   mp_err(*redux)(mp_int *x, const mp_int *m, const mp_int *mu);

   /* find window size */
   x = mp_count_bits(X);
   if (x <= 7) {
      winsize = 2;
   } else if (x <= 36) {
      winsize = 3;
   } else if (x <= 140) {
      winsize = 4;
   } else if (x <= 450) {
      winsize = 5;
   } else if (x <= 1303) {
      winsize = 6;
   } else if (x <= 3529) {
      winsize = 7;
   } else {
      winsize = 8;
   }

   winsize = MAX_WINSIZE ? MP_MIN(MAX_WINSIZE, winsize) : winsize;

   /* init M array */
   /* init first cell */
   if ((err = mp_init(&M[1])) != MP_OKAY) {
      return err;
   }

   /* now init the second half of the array */
   for (x = 1<<(winsize-1); x < (1 << winsize); x++) {
      if ((err = mp_init(&M[x])) != MP_OKAY) {
         for (y = 1<<(winsize-1); y < x; y++) {
            mp_clear(&M[y]);
         }
         mp_clear(&M[1]);
         return err;
      }
   }

   /* create mu, used for Barrett reduction */
   if ((err = mp_init(&mu)) != MP_OKAY)                           goto LBL_M;

   if (redmode == 0) {
      if ((err = mp_reduce_setup(&mu, P)) != MP_OKAY)             goto LBL_MU;
      redux = mp_reduce;
   } else {
      if ((err = mp_reduce_2k_setup_l(P, &mu)) != MP_OKAY)        goto LBL_MU;
      redux = mp_reduce_2k_l;
   }

   /* create M table
    *
    * The M table contains powers of the base,
    * e.g. M[x] = G**x mod P
    *
    * The first half of the table is not
    * computed though accept for M[0] and M[1]
    */
   if ((err = mp_mod(G, P, &M[1])) != MP_OKAY)                    goto LBL_MU;

   /* compute the value at M[1<<(winsize-1)] by squaring
    * M[1] (winsize-1) times
    */
   if ((err = mp_copy(&M[1], &M[(size_t)1 << (winsize - 1)])) != MP_OKAY) goto LBL_MU;

   for (x = 0; x < (winsize - 1); x++) {
      /* square it */
      if ((err = mp_sqr(&M[(size_t)1 << (winsize - 1)],
                        &M[(size_t)1 << (winsize - 1)])) != MP_OKAY) goto LBL_MU;

      /* reduce modulo P */
      if ((err = redux(&M[(size_t)1 << (winsize - 1)], P, &mu)) != MP_OKAY) goto LBL_MU;
   }

   /* create upper table, that is M[x] = M[x-1] * M[1] (mod P)
    * for x = (2**(winsize - 1) + 1) to (2**winsize - 1)
    */
   for (x = (1 << (winsize - 1)) + 1; x < (1 << winsize); x++) {
      if ((err = mp_mul(&M[x - 1], &M[1], &M[x])) != MP_OKAY)     goto LBL_MU;
      if ((err = redux(&M[x], P, &mu)) != MP_OKAY)                goto LBL_MU;
   }

   /* setup result */
   if ((err = mp_init(&res)) != MP_OKAY)                          goto LBL_MU;
   mp_set(&res, 1uL);

   /* set initial mode and bit cnt */
   mode   = 0;
   bitcnt = 1;
   buf    = 0;
   digidx = X->used - 1;
   bitcpy = 0;
   bitbuf = 0;

   for (;;) {
      /* grab next digit as required */
      if (--bitcnt == 0) {
         /* if digidx == -1 we are out of digits */
         if (digidx == -1) {
            break;
         }
         /* read next digit and reset the bitcnt */
         buf    = X->dp[digidx--];
         bitcnt = (int)MP_DIGIT_BIT;
      }

      /* grab the next msb from the exponent */
      y     = (buf >> (mp_digit)(MP_DIGIT_BIT - 1)) & 1uL;
      buf <<= (mp_digit)1;

      /* if the bit is zero and mode == 0 then we ignore it
       * These represent the leading zero bits before the first 1 bit
       * in the exponent.  Technically this opt is not required but it
       * does lower the # of trivial squaring/reductions used
       */
      if ((mode == 0) && (y == 0)) {
         continue;
      }

      /* if the bit is zero and mode == 1 then we square */
      if ((mode == 1) && (y == 0)) {
         if ((err = mp_sqr(&res, &res)) != MP_OKAY)               goto LBL_RES;
         if ((err = redux(&res, P, &mu)) != MP_OKAY)              goto LBL_RES;
         continue;
      }

      /* else we add it to the window */
      bitbuf |= (y << (winsize - ++bitcpy));
      mode    = 2;

      if (bitcpy == winsize) {
         /* ok window is filled so square as required and multiply  */
         /* square first */
         for (x = 0; x < winsize; x++) {
            if ((err = mp_sqr(&res, &res)) != MP_OKAY)            goto LBL_RES;
            if ((err = redux(&res, P, &mu)) != MP_OKAY)           goto LBL_RES;
         }

         /* then multiply */
         if ((err = mp_mul(&res, &M[bitbuf], &res)) != MP_OKAY)  goto LBL_RES;
         if ((err = redux(&res, P, &mu)) != MP_OKAY)             goto LBL_RES;

         /* empty window and reset */
         bitcpy = 0;
         bitbuf = 0;
         mode   = 1;
      }
   }

   /* if bits remain then square/multiply */
   if ((mode == 2) && (bitcpy > 0)) {
      /* square then multiply if the bit is set */
      for (x = 0; x < bitcpy; x++) {
         if ((err = mp_sqr(&res, &res)) != MP_OKAY)               goto LBL_RES;
         if ((err = redux(&res, P, &mu)) != MP_OKAY)              goto LBL_RES;

         bitbuf <<= 1;
         if ((bitbuf & (1 << winsize)) != 0) {
            /* then multiply */
            if ((err = mp_mul(&res, &M[1], &res)) != MP_OKAY)     goto LBL_RES;
            if ((err = redux(&res, P, &mu)) != MP_OKAY)           goto LBL_RES;
         }
      }
   }

   mp_exch(&res, Y);
   err = MP_OKAY;
LBL_RES:
   mp_clear(&res);
LBL_MU:
   mp_clear(&mu);
LBL_M:
   mp_clear(&M[1]);
   for (x = 1<<(winsize-1); x < (1 << winsize); x++) {
      mp_clear(&M[x]);
   }
   return err;
}



/* determines if reduce_2k_l can be used */
bool mp_reduce_is_2k_l(const mp_int *a)
{
   if (mp_iszero(a)) {
      return false;
   } else if (a->used == 1) {
      return true;
   } else if (a->used > 1) {
      /* if more than half of the digits are -1 we're sold */
      int ix, iy;
      for (iy = ix = 0; ix < a->used; ix++) {
         if (a->dp[ix] == MP_DIGIT_MAX) {
            ++iy;
         }
      }
      return (iy >= (a->used/2));
   } else {
      return false;
   }
}


/* b = |a|
 *
 * Simple function copies the input and fixes the sign to positive
 */
mp_err mp_abs(const mp_int *a, mp_int *b)
{
   mp_err err;

   /* copy a to b */
   if ((err = mp_copy(a, b)) != MP_OKAY) {
      return err;
   }

   /* force the sign of b to positive */
   b->sign = MP_ZPOS;

   return MP_OKAY;
}

/* hac 14.61, pp608 */
mp_err mp_invmod(const mp_int *a, const mp_int *b, mp_int *c)
{
   /* for all n in N and n > 0, n = 0 mod 1 */
   if (!mp_isneg(a) && mp_cmp_d(b, 1uL) == MP_EQ) {
      mp_zero(c);
      return MP_OKAY;
   }

   /* b cannot be negative and has to be >1 */
   if (mp_isneg(b) || (mp_cmp_d(b, 1uL) != MP_GT)) {
      return MP_VAL;
   }

   /* if the modulus is odd we can use a faster routine instead */
   if (MP_HAS(S_MP_INVMOD_ODD) && mp_isodd(b)) {
      return s_mp_invmod_odd(a, b, c);
   }

   return MP_HAS(S_MP_INVMOD)
          ? s_mp_invmod(a, b, c)
          : MP_VAL;
}



/* c = a mod b, 0 <= c < b if b > 0, b < c <= 0 if b < 0 */
mp_err mp_mod(const mp_int *a, const mp_int *b, mp_int *c)
{
   mp_err err;
   if ((err = mp_div(a, b, NULL, c)) != MP_OKAY) {
      return err;
   }
   return mp_iszero(c) || (c->sign == b->sign) ? MP_OKAY : mp_add(b, c, c);
}



/* c = a * a (mod b) */
mp_err mp_sqrmod(const mp_int *a, const mp_int *b, mp_int *c)
{
   mp_err err;
   if ((err = mp_sqr(a, c)) != MP_OKAY) {
      return err;
   }
   return mp_mod(c, b, c);
}

/* this is a shell function that calls either the normal or Montgomery
 * exptmod functions.  Originally the call to the montgomery code was
 * embedded in the normal function but that wasted alot of stack space
 * for nothing (since 99% of the time the Montgomery code would be called)
 */
mp_err mp_exptmod(const mp_int *G, const mp_int *X, const mp_int *P, mp_int *Y)
{
   int dr;
   //char buff[4097];
   //////////////
   /*
   printf("printing mp_int\n");
   mp_to_decimal(G,buff,sizeof(buff));
   printf("\nG===\n%s\n",buff);

   mp_to_decimal(X,buff,sizeof(buff));
   printf("\nX===\n%s\n",buff);

   mp_to_decimal(P,buff,sizeof(buff));
   printf("\nP===\n%s\n",buff);
   */
   ////////////////

   /* modulus P must be positive */
   if (mp_isneg(P)) {
      //printf("\nmodulus p is negative\n");
      return MP_VAL;
   }

   /* if exponent X is negative we have to recurse */
   if (mp_isneg(X)) {
      mp_int tmpG, tmpX;
      mp_err err;

      if (!MP_HAS(MP_INVMOD)) {
         return MP_VAL;
      }

      if ((err = mp_init_multi(&tmpG, &tmpX, NULL)) != MP_OKAY) {
         return err;
      }

      /* first compute 1/G mod P */
      if ((err = mp_invmod(G, P, &tmpG)) != MP_OKAY) {
         goto LBL_ERR;
      }

      /* now get |X| */
      if ((err = mp_abs(X, &tmpX)) != MP_OKAY) {
         goto LBL_ERR;
      }

      /* and now compute (1/G)**|X| instead of G**X [X < 0] */
      err = mp_exptmod(&tmpG, &tmpX, P, Y);
LBL_ERR:
      mp_clear_multi(&tmpG, &tmpX, NULL);
      return err;
   }

   /* modified diminished radix reduction */
   if (MP_HAS(MP_REDUCE_IS_2K_L) && MP_HAS(MP_REDUCE_2K_L) && MP_HAS(S_MP_EXPTMOD) &&
       mp_reduce_is_2k_l(P)) {
      return s_mp_exptmod(G, X, P, Y, 1);
   }
  // printf("\nthis part didnt run\n");
   /* is it a DR modulus? default to no */
   dr = (MP_HAS(MP_DR_IS_MODULUS) && mp_dr_is_modulus(P)) ? 1 : 0;

   /* if not, is it a unrestricted DR modulus? */
   if (MP_HAS(MP_REDUCE_IS_2K) && (dr == 0)) {
      dr = (mp_reduce_is_2k(P)) ? 2 : 0;
   }

   /* if the modulus is odd or dr != 0 use the montgomery method */
   if (MP_HAS(S_MP_EXPTMOD_FAST) && (mp_isodd(P) || (dr != 0))) {
      return s_mp_exptmod_fast(G, X, P, Y, dr);
   }
 //printf("\nthis part didnt run\n");
   /* otherwise use the generic Barrett reduction technique */
   if (MP_HAS(S_MP_EXPTMOD)) {
      return s_mp_exptmod(G, X, P, Y, 0);
   }
// printf("\nthis part didnt run\n");
   /* no exptmod for evens */
   return MP_VAL;
}

static const char lnz[16] = {
   4, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0
};

/* Counts the number of lsbs which are zero before the first zero bit */
int mp_cnt_lsb(const mp_int *a)
{
   int x;
   mp_digit q;

   /* easy out */
   if (mp_iszero(a)) {
      return 0;
   }

   /* scan lower digits until non-zero */
   for (x = 0; (x < a->used) && (a->dp[x] == 0u); x++) {}
   q = a->dp[x];
   x *= MP_DIGIT_BIT;

   /* now scan this digit until a 1 is found */
   if ((q & 1u) == 0u) {
      mp_digit p;
      do {
         p = q & 15u;
         x += lnz[p];
         q >>= 4;
      } while (p == 0u);
   }
   return x;
}


/* creates "a" then copies b into it */
mp_err mp_init_copy(mp_int *a, const mp_int *b)
{
   mp_err     err;

   if ((err = mp_init_size(a, b->used)) != MP_OKAY) {
      return err;
   }

   if ((err = mp_copy(b, a)) != MP_OKAY) {
      mp_clear(a);
   }

   return err;
}


/* init a new mp_int */
mp_err mp_init(mp_int *a)
{
   /* allocate memory required and clear it */
   a->dp = (mp_digit *) MP_CALLOC((size_t)MP_DEFAULT_DIGIT_COUNT, sizeof(mp_digit));
   if (a->dp == NULL) {
      return MP_MEM;
   }

   /* set the used to zero, allocated digits to the default precision
    * and sign to positive */
   a->used  = 0;
   a->alloc = MP_DEFAULT_DIGIT_COUNT;//32
   a->sign  = MP_ZPOS;

   return MP_OKAY;
}


/* single digit division (based on routine from MPI) */
mp_err mp_div_d(const mp_int *a, mp_digit b, mp_int *c, mp_digit *d)
{
   mp_int  q;
   mp_word w;
   mp_err err;
   int ix;

   /* cannot divide by zero */
   if (b == 0u) {
      return MP_VAL;
   }

   /* quick outs */
   if ((b == 1u) || mp_iszero(a)) {
      if (d != NULL) {
         *d = 0;
      }
      if (c != NULL) {
         return mp_copy(a, c);
      }
      return MP_OKAY;
   }

   /* power of two ? */
   if (MP_HAS(MP_DIV_2) && (b == 2u)) {
      if (d != NULL) {
         *d = mp_isodd(a) ? 1u : 0u;
      }
      return (c == NULL) ? MP_OKAY : mp_div_2(a, c);
   }
   if (MP_HAS(MP_DIV_2D) && MP_IS_2EXPT(b)) {
      ix = 1;
      while ((ix < MP_DIGIT_BIT) && (b != (((mp_digit)1)<<ix))) {
         ix++;
      }
      if (d != NULL) {
         *d = a->dp[0] & (((mp_digit)1<<(mp_digit)ix) - 1uL);
      }
      return (c == NULL) ? MP_OKAY : mp_div_2d(a, ix, c, NULL);
   }

   /* three? */
   if (MP_HAS(S_MP_DIV_3) && (b == 3u)) {
      return s_mp_div_3(a, c, d);
   }

   /* no easy answer [c'est la vie].  Just division */
   if ((err = mp_init_size(&q, a->used)) != MP_OKAY) {
      return err;
   }

   q.used = a->used;
   q.sign = a->sign;
   w = 0;
   for (ix = a->used; ix --> 0;) {
      mp_digit t = 0;
      w = (w << (mp_word)MP_DIGIT_BIT) | (mp_word)a->dp[ix];
      if (w >= b) {
         t = (mp_digit)(w / b);
         w -= (mp_word)t * (mp_word)b;
      }
      q.dp[ix] = t;
   }

   if (d != NULL) {
      *d = (mp_digit)w;
   }

   if (c != NULL) {
      mp_clamp(&q);
      mp_exch(&q, c);
   }
   mp_clear(&q);

   return MP_OKAY;
}


/* c = a mod b, 0 <= c < b  */
#define mp_mod_d(a, b, c) mp_div_d((a), (b), NULL, (c))
/* clear one (frees)  */
void mp_clear(mp_int *a)
{
   /* only do anything if a hasn't been freed previously */
   if (a->dp != NULL) {
      /* free ram */
      MP_FREE_DIGS(a->dp, a->alloc);

      /* reset members to make debugging easier */
      a->dp    = NULL;
      a->alloc = a->used = 0;
      a->sign  = MP_ZPOS;
   }
}


/* compare a digit */
mp_ord mp_cmp_d(const mp_int *a, mp_digit b)
{
   /* compare based on sign */
   if (mp_isneg(a)) {
      return MP_LT;
   }

   /* compare based on magnitude */
   if (a->used > 1) {
      return MP_GT;
   }

   /* compare the only digit of a to b */
   if (a->dp[0] != b) {
      return a->dp[0] > b ? MP_GT : MP_LT;
   }

   return MP_EQ;
}

/* shift right by a certain bit count (store quotient in c, optional remainder in d) */
mp_err mp_div_2d(const mp_int *a, int b, mp_int *c, mp_int *d)
{
   mp_err err;

   if (b < 0) {
      return MP_VAL;
   }

   if ((err = mp_copy(a, c)) != MP_OKAY) {
      return err;
   }

   /* 'a' should not be used after here - it might be the same as d */

   /* get the remainder */
   if (d != NULL) {
      if ((err = mp_mod_2d(a, b, d)) != MP_OKAY) {
         return err;
      }
   }

   /* shift by as many digits in the bit count */
   if (b >= MP_DIGIT_BIT) {
      mp_rshd(c, b / MP_DIGIT_BIT);
   }

   /* shift any bit count < MP_DIGIT_BIT */
   b %= MP_DIGIT_BIT;
   if (b != 0u) {
      int x;
      mp_digit r, mask, shift;

      /* mask */
      mask = ((mp_digit)1 << b) - 1uL;

      /* shift for lsb */
      shift = (mp_digit)(MP_DIGIT_BIT - b);

      /* carry */
      r = 0;
      for (x = c->used; x --> 0;) {
         /* get the lower  bits of this word in a temp */
         mp_digit rr = c->dp[x] & mask;

         /* shift the current word and mix in the carry bits from the previous word */
         c->dp[x] = (c->dp[x] >> b) | (r << shift);

         /* set the carry to the carry bits of the current word found above */
         r = rr;
      }
   }
   mp_clamp(c);
   return MP_OKAY;
}



/* returns the number of bits in an int */
int mp_count_bits(const mp_int *a)
{
   int     r;
   mp_digit q;

   /* shortcut */
   if (mp_iszero(a)) {
      return 0;
   }

   /* get number of digits and add that */
   r = (a->used - 1) * MP_DIGIT_BIT;

   /* take the last digit and count the bits in it */
   q = a->dp[a->used - 1];
   while (q > 0u) {
      ++r;
      q >>= 1u;
   }
   return r;
}


/* set to a digit */
void mp_set(mp_int *a, mp_digit b)
{
   int oldused = a->used;
   a->dp[0] = b & MP_MASK;
   a->sign  = MP_ZPOS;
   a->used  = (a->dp[0] != 0u) ? 1 : 0;
   s_mp_zero_digs(a->dp + a->used, oldused - a->used);
}


/* compare two ints (signed)*/
mp_ord mp_cmp(const mp_int *a, const mp_int *b)
{
   /* compare based on sign */
   if (a->sign != b->sign) {
      return mp_isneg(a) ? MP_LT : MP_GT;
   }

   /* if negative compare opposite direction */
   if (mp_isneg(a)) {
      MP_EXCH(const mp_int *, a, b);
   }

   return mp_cmp_mag(a, b);
}


mp_err mp_read_radix(mp_int *a, const char *str, int radix)
{
   mp_err   err;
   mp_sign  sign = MP_ZPOS;

   /* make sure the radix is ok */
   if ((radix < 2) || (radix > 64)) {
      return MP_VAL;
   }

   /* if the leading digit is a
    * minus set the sign to negative.
    */
   if (*str == '-') {
      ++str;
      sign = MP_NEG;
   }

   /* set the integer to the default of zero */
   mp_zero(a);

   /* process each digit of the string */
   while (*str != '\0') {
      /* if the radix <= 36 the conversion is case insensitive
       * this allows numbers like 1AB and 1ab to represent the same  value
       * [e.g. in hex]
       */
      uint8_t y;
      char ch = (radix <= 36) ? (char)MP_TOUPPER((int)*str) : *str;
      unsigned pos = (unsigned)(ch - '+');
      if (MP_RADIX_MAP_REVERSE_SIZE <= pos) {
         break;
      }
      y = s_mp_radix_map_reverse[pos];

      /* if the char was found in the map
       * and is less than the given radix add it
       * to the number, otherwise exit the loop.
       */
      if (y >= radix) {
         break;
      }
      if ((err = mp_mul_d(a, (mp_digit)radix, a)) != MP_OKAY) {
         return err;
      }
      if ((err = mp_add_d(a, y, a)) != MP_OKAY) {
         return err;
      }
      ++str;
   }

   /* if an illegal character was found, fail. */
   if ((*str != '\0') && (*str != '\r') && (*str != '\n')) {
      return MP_VAL;
   }

   /* set the sign only if a != 0 */
   if (!mp_iszero(a)) {
      a->sign = sign;
   }
   return MP_OKAY;
}


/*
 * multiply bigint a with int d and put the result in c
 * Like mp_mul_d() but with a signed long as the small input
 */
static mp_err s_mul_si(const mp_int *a, int32_t d, mp_int *c)
{
   mp_int t;
   mp_err err;

   if ((err = mp_init(&t)) != MP_OKAY) {
      return err;
   }

   /*
    * mp_digit might be smaller than a long, which excludes
    * the use of mp_mul_d() here.
    */
   mp_set_i32(&t, d);
   err = mp_mul(a, &t, c);
   mp_clear(&t);
   return err;
}
/*
    Strong Lucas-Selfridge test.
    returns true if it is a strong L-S prime, false if it is composite

    Code ported from  Thomas Ray Nicely's implementation of the BPSW test
    at http://www.trnicely.net/misc/bpsw.html

    Freeware copyright (C) 2016 Thomas R. Nicely <http://www.trnicely.net>.
    Released into the public domain by the author, who disclaims any legal
    liability arising from its use

    The multi-line comments are made by Thomas R. Nicely and are copied verbatim.
    Additional comments marked "CZ" (without the quotes) are by the code-portist.

    (If that name sounds familiar, he is the guy who found the fdiv bug in the
     Pentium (P5x, I think) Intel processor)
*/
mp_err mp_prime_strong_lucas_selfridge(const mp_int *a, bool *result)
{
   /* CZ TODO: choose better variable names! */
   mp_int Dz, gcd, Np1, Uz, Vz, U2mz, V2mz, Qmz, Q2mz, Qkdz, T1z, T2z, T3z, T4z, Q2kdz;
   int32_t D, Ds, sign, P, Q, r, s, u, Nbits;
   int J;
   mp_err err;
   bool oddness;

   *result = false;
   /*
   Find the first element D in the sequence {5, -7, 9, -11, 13, ...}
   such that Jacobi(D,N) = -1 (Selfridge's algorithm). Theory
   indicates that, if N is not a perfect square, D will "nearly
   always" be "small." Just in case, an overflow trap for D is
   included.
   */

   if ((err = mp_init_multi(&Dz, &gcd, &Np1, &Uz, &Vz, &U2mz, &V2mz, &Qmz, &Q2mz, &Qkdz, &T1z, &T2z, &T3z, &T4z, &Q2kdz,
                            NULL)) != MP_OKAY) {
      return err;
   }

   D = 5;
   sign = 1;

   for (;;) {
      Ds   = sign * D;
      sign = -sign;
      mp_set_u32(&Dz, (uint32_t)D);
      if ((err = mp_gcd(a, &Dz, &gcd)) != MP_OKAY)                goto LBL_LS_ERR;

      /* if 1 < GCD < N then N is composite with factor "D", and
         Jacobi(D,N) is technically undefined (but often returned
         as zero). */
      if ((mp_cmp_d(&gcd, 1uL) == MP_GT) && (mp_cmp(&gcd, a) == MP_LT)) {
         goto LBL_LS_ERR;
      }
      if (Ds < 0) {
         Dz.sign = MP_NEG;
      }
      if ((err = mp_kronecker(&Dz, a, &J)) != MP_OKAY)            goto LBL_LS_ERR;

      if (J == -1) {
         break;
      }
      D += 2;

      if (D > (INT_MAX - 2)) {
         err = MP_VAL;
         goto LBL_LS_ERR;
      }
   }



   P = 1;              /* Selfridge's choice */
   Q = (1 - Ds) / 4;   /* Required so D = P*P - 4*Q */

   /* NOTE: The conditions (a) N does not divide Q, and
      (b) D is square-free or not a perfect square, are included by
      some authors; e.g., "Prime numbers and computer methods for
      factorization," Hans Riesel (2nd ed., 1994, Birkhauser, Boston),
      p. 130. For this particular application of Lucas sequences,
      these conditions were found to be immaterial. */

   /* Now calculate N - Jacobi(D,N) = N + 1 (even), and calculate the
      odd positive integer d and positive integer s for which
      N + 1 = 2^s*d (similar to the step for N - 1 in Miller's test).
      The strong Lucas-Selfridge test then returns N as a strong
      Lucas probable prime (slprp) if any of the following
      conditions is met: U_d=0, V_d=0, V_2d=0, V_4d=0, V_8d=0,
      V_16d=0, ..., etc., ending with V_{2^(s-1)*d}=V_{(N+1)/2}=0
      (all equalities mod N). Thus d is the highest index of U that
      must be computed (since V_2m is independent of U), compared
      to U_{N+1} for the standard Lucas-Selfridge test; and no
      index of V beyond (N+1)/2 is required, just as in the
      standard Lucas-Selfridge test. However, the quantity Q^d must
      be computed for use (if necessary) in the latter stages of
      the test. The result is that the strong Lucas-Selfridge test
      has a running time only slightly greater (order of 10 %) than
      that of the standard Lucas-Selfridge test, while producing
      only (roughly) 30 % as many pseudoprimes (and every strong
      Lucas pseudoprime is also a standard Lucas pseudoprime). Thus
      the evidence indicates that the strong Lucas-Selfridge test is
      more effective than the standard Lucas-Selfridge test, and a
      Baillie-PSW test based on the strong Lucas-Selfridge test
      should be more reliable. */

   if ((err = mp_add_d(a, 1uL, &Np1)) != MP_OKAY)                 goto LBL_LS_ERR;
   s = mp_cnt_lsb(&Np1);

   /* CZ
    * This should round towards zero because
    * Thomas R. Nicely used GMP's mpz_tdiv_q_2exp()
    * and mp_div_2d() is equivalent. Additionally:
    * dividing an even number by two does not produce
    * any leftovers.
    */
   if ((err = mp_div_2d(&Np1, s, &Dz, NULL)) != MP_OKAY)          goto LBL_LS_ERR;
   /* We must now compute U_d and V_d. Since d is odd, the accumulated
      values U and V are initialized to U_1 and V_1 (if the target
      index were even, U and V would be initialized instead to U_0=0
      and V_0=2). The values of U_2m and V_2m are also initialized to
      U_1 and V_1; the FOR loop calculates in succession U_2 and V_2,
      U_4 and V_4, U_8 and V_8, etc. If the corresponding bits
      (1, 2, 3, ...) of t are on (the zero bit having been accounted
      for in the initialization of U and V), these values are then
      combined with the previous totals for U and V, using the
      composition formulas for addition of indices. */

   mp_set(&Uz, 1uL);    /* U=U_1 */
   mp_set(&Vz, (mp_digit)P);    /* V=V_1 */
   mp_set(&U2mz, 1uL);  /* U_1 */
   mp_set(&V2mz, (mp_digit)P);  /* V_1 */

   mp_set_i32(&Qmz, Q);
   if ((err = mp_mul_2(&Qmz, &Q2mz)) != MP_OKAY)                  goto LBL_LS_ERR;
   /* Initializes calculation of Q^d */
   mp_set_i32(&Qkdz, Q);

   Nbits = mp_count_bits(&Dz);

   for (u = 1; u < Nbits; u++) { /* zero bit off, already accounted for */
      /* Formulas for doubling of indices (carried out mod N). Note that
       * the indices denoted as "2m" are actually powers of 2, specifically
       * 2^(ul-1) beginning each loop and 2^ul ending each loop.
       *
       * U_2m = U_m*V_m
       * V_2m = V_m*V_m - 2*Q^m
       */

      if ((err = mp_mul(&U2mz, &V2mz, &U2mz)) != MP_OKAY)         goto LBL_LS_ERR;
      if ((err = mp_mod(&U2mz, a, &U2mz)) != MP_OKAY)             goto LBL_LS_ERR;
      if ((err = mp_sqr(&V2mz, &V2mz)) != MP_OKAY)                goto LBL_LS_ERR;
      if ((err = mp_sub(&V2mz, &Q2mz, &V2mz)) != MP_OKAY)         goto LBL_LS_ERR;
      if ((err = mp_mod(&V2mz, a, &V2mz)) != MP_OKAY)             goto LBL_LS_ERR;

      /* Must calculate powers of Q for use in V_2m, also for Q^d later */
      if ((err = mp_sqr(&Qmz, &Qmz)) != MP_OKAY)                  goto LBL_LS_ERR;

      /* prevents overflow */ /* CZ  still necessary without a fixed prealloc'd mem.? */
      if ((err = mp_mod(&Qmz, a, &Qmz)) != MP_OKAY)               goto LBL_LS_ERR;
      if ((err = mp_mul_2(&Qmz, &Q2mz)) != MP_OKAY)               goto LBL_LS_ERR;

      if (s_mp_get_bit(&Dz, u)) {
         /* Formulas for addition of indices (carried out mod N);
          *
          * U_(m+n) = (U_m*V_n + U_n*V_m)/2
          * V_(m+n) = (V_m*V_n + D*U_m*U_n)/2
          *
          * Be careful with division by 2 (mod N)!
          */
         if ((err = mp_mul(&U2mz, &Vz, &T1z)) != MP_OKAY)         goto LBL_LS_ERR;
         if ((err = mp_mul(&Uz, &V2mz, &T2z)) != MP_OKAY)         goto LBL_LS_ERR;
         if ((err = mp_mul(&V2mz, &Vz, &T3z)) != MP_OKAY)         goto LBL_LS_ERR;
         if ((err = mp_mul(&U2mz, &Uz, &T4z)) != MP_OKAY)         goto LBL_LS_ERR;
         if ((err = s_mul_si(&T4z, Ds, &T4z)) != MP_OKAY)      goto LBL_LS_ERR;
         if ((err = mp_add(&T1z, &T2z, &Uz)) != MP_OKAY)          goto LBL_LS_ERR;
         if (mp_isodd(&Uz)) {
            if ((err = mp_add(&Uz, a, &Uz)) != MP_OKAY)           goto LBL_LS_ERR;
         }
         /* CZ
          * This should round towards negative infinity because
          * Thomas R. Nicely used GMP's mpz_fdiv_q_2exp().
          * But mp_div_2() does not do so, it is truncating instead.
          */
         oddness = mp_isodd(&Uz);
         if ((err = mp_div_2(&Uz, &Uz)) != MP_OKAY)               goto LBL_LS_ERR;
         if (mp_isneg(&Uz) && oddness) {
            if ((err = mp_sub_d(&Uz, 1uL, &Uz)) != MP_OKAY)       goto LBL_LS_ERR;
         }
         if ((err = mp_add(&T3z, &T4z, &Vz)) != MP_OKAY)          goto LBL_LS_ERR;
         if (mp_isodd(&Vz)) {
            if ((err = mp_add(&Vz, a, &Vz)) != MP_OKAY)           goto LBL_LS_ERR;
         }
         oddness = mp_isodd(&Vz);
         if ((err = mp_div_2(&Vz, &Vz)) != MP_OKAY)               goto LBL_LS_ERR;
         if (mp_isneg(&Vz) && oddness) {
            if ((err = mp_sub_d(&Vz, 1uL, &Vz)) != MP_OKAY)       goto LBL_LS_ERR;
         }
         if ((err = mp_mod(&Uz, a, &Uz)) != MP_OKAY)              goto LBL_LS_ERR;
         if ((err = mp_mod(&Vz, a, &Vz)) != MP_OKAY)              goto LBL_LS_ERR;

         /* Calculating Q^d for later use */
         if ((err = mp_mul(&Qkdz, &Qmz, &Qkdz)) != MP_OKAY)       goto LBL_LS_ERR;
         if ((err = mp_mod(&Qkdz, a, &Qkdz)) != MP_OKAY)          goto LBL_LS_ERR;
      }
   }

   /* If U_d or V_d is congruent to 0 mod N, then N is a prime or a
      strong Lucas pseudoprime. */
   if (mp_iszero(&Uz) || mp_iszero(&Vz)) {
      *result = true;
      goto LBL_LS_ERR;
   }

   /* NOTE: Ribenboim ("The new book of prime number records," 3rd ed.,
      1995/6) omits the condition V0 on p.142, but includes it on
      p. 130. The condition is NECESSARY; otherwise the test will
      return false negatives---e.g., the primes 29 and 2000029 will be
      returned as composite. */

   /* Otherwise, we must compute V_2d, V_4d, V_8d, ..., V_{2^(s-1)*d}
      by repeated use of the formula V_2m = V_m*V_m - 2*Q^m. If any of
      these are congruent to 0 mod N, then N is a prime or a strong
      Lucas pseudoprime. */

   /* Initialize 2*Q^(d*2^r) for V_2m */
   if ((err = mp_mul_2(&Qkdz, &Q2kdz)) != MP_OKAY)                goto LBL_LS_ERR;

   for (r = 1; r < s; r++) {
      if ((err = mp_sqr(&Vz, &Vz)) != MP_OKAY)                    goto LBL_LS_ERR;
      if ((err = mp_sub(&Vz, &Q2kdz, &Vz)) != MP_OKAY)            goto LBL_LS_ERR;
      if ((err = mp_mod(&Vz, a, &Vz)) != MP_OKAY)                 goto LBL_LS_ERR;
      if (mp_iszero(&Vz)) {
         *result = true;
         goto LBL_LS_ERR;
      }
      /* Calculate Q^{d*2^r} for next r (final iteration irrelevant). */
      if (r < (s - 1)) {
         if ((err = mp_sqr(&Qkdz, &Qkdz)) != MP_OKAY)             goto LBL_LS_ERR;
         if ((err = mp_mod(&Qkdz, a, &Qkdz)) != MP_OKAY)          goto LBL_LS_ERR;
         if ((err = mp_mul_2(&Qkdz, &Q2kdz)) != MP_OKAY)          goto LBL_LS_ERR;
      }
   }
LBL_LS_ERR:
   mp_clear_multi(&Q2kdz, &T4z, &T3z, &T2z, &T1z, &Qkdz, &Q2mz, &Qmz, &V2mz, &U2mz, &Vz, &Uz, &Np1, &gcd, &Dz, NULL);
   return err;
}



/* Miller-Rabin test of "a" to the base of "b" as described in
 * HAC pp. 139 Algorithm 4.24
 *
 * Sets result to 0 if definitely composite or 1 if probably prime.
 * Randomly the chance of error is no more than 1/4 and often
 * very much lower.
 */
mp_err mp_prime_miller_rabin(const mp_int *a, const mp_int *b, bool *result)
{
   mp_int  n1, y, r;
   mp_err  err;
   int     s, j;

   /* ensure b > 1 */
   if (mp_cmp_d(b, 1uL) != MP_GT) {
      return MP_VAL;
   }

   /* get n1 = a - 1 */
   if ((err = mp_init_copy(&n1, a)) != MP_OKAY) {
     // printf("\nerror in mp_init_copy line 5438 %d\n",err);
      return err;
   }
   if ((err = mp_sub_d(&n1, 1uL, &n1)) != MP_OKAY) {
    //  printf("\nerror in mp_sub_d line 5442 %d\n",err);
      goto LBL_ERR1;
   }

   /* set 2**s * r = n1 */
   if ((err = mp_init_copy(&r, &n1)) != MP_OKAY) {
     // printf("\nerror in mp_init_copy line 5448 %d\n",err);
      goto LBL_ERR1;
   }

   /* count the number of least significant bits
    * which are zero
    */
   s = mp_cnt_lsb(&r);

   /* now divide n - 1 by 2**s */
   if ((err = mp_div_2d(&r, s, &r, NULL)) != MP_OKAY) {
      //printf("\nerror in mp_div_2d line 5459 %d\n",err);
      goto LBL_ERR2;
   }

   /* compute y = b**r mod a */
   if ((err = mp_init(&y)) != MP_OKAY) {
     // printf("\nerror in mp_init line 5465 %d\n",err);
      goto LBL_ERR2;
   }
   if ((err = mp_exptmod(b, &r, a, &y)) != MP_OKAY) {
     // printf("\nerror in mp_exptmod line 5469 %d\n",err);
      goto LBL_END;
   }

   /* if y != 1 and y != n1 do */
   if ((mp_cmp_d(&y, 1uL) != MP_EQ) && (mp_cmp(&y, &n1) != MP_EQ)) {
      j = 1;
      /* while j <= s-1 and y != n1 */
      while ((j <= (s - 1)) && (mp_cmp(&y, &n1) != MP_EQ)) {
         if ((err = mp_sqrmod(&y, a, &y)) != MP_OKAY) {
            goto LBL_END;
         }

         /* if y == 1 then composite */
         if (mp_cmp_d(&y, 1uL) == MP_EQ) {
            *result = false;
            goto LBL_END;
         }

         ++j;
      }

      /* if y != n1 then composite */
      if (mp_cmp(&y, &n1) != MP_EQ) {
         *result = false;
         goto LBL_END;
      }
   }

   /* probably prime now */
   *result = true;

LBL_END:
   mp_clear(&y);
LBL_ERR2:
   mp_clear(&r);
LBL_ERR1:
   mp_clear(&n1);
   return err;
}



/* initialize and set a digit */
mp_err mp_init_set(mp_int *a, mp_digit b)
{
   mp_err err;
   if ((err = mp_init(a)) != MP_OKAY) {
      return err;
   }
   mp_set(a, b);
   return err;
}



/* determines if an integers is divisible by one
 * of the first PRIME_SIZE primes or not
 *
 * sets result to 0 if not, 1 if yes
 */
int MP_PRIME_TAB_SIZE=256;

const mp_digit s_mp_prime_tab[] = {
   0x0002, 0x0003, 0x0005, 0x0007, 0x000B, 0x000D, 0x0011, 0x0013,
   0x0017, 0x001D, 0x001F, 0x0025, 0x0029, 0x002B, 0x002F, 0x0035,
   0x003B, 0x003D, 0x0043, 0x0047, 0x0049, 0x004F, 0x0053, 0x0059,
   0x0061, 0x0065, 0x0067, 0x006B, 0x006D, 0x0071, 0x007F, 0x0083,
   0x0089, 0x008B, 0x0095, 0x0097, 0x009D, 0x00A3, 0x00A7, 0x00AD,
   0x00B3, 0x00B5, 0x00BF, 0x00C1, 0x00C5, 0x00C7, 0x00D3, 0x00DF,
   0x00E3, 0x00E5, 0x00E9, 0x00EF, 0x00F1, 0x00FB, 0x0101, 0x0107,
   0x010D, 0x010F, 0x0115, 0x0119, 0x011B, 0x0125, 0x0133, 0x0137,

   0x0139, 0x013D, 0x014B, 0x0151, 0x015B, 0x015D, 0x0161, 0x0167,
   0x016F, 0x0175, 0x017B, 0x017F, 0x0185, 0x018D, 0x0191, 0x0199,
   0x01A3, 0x01A5, 0x01AF, 0x01B1, 0x01B7, 0x01BB, 0x01C1, 0x01C9,
   0x01CD, 0x01CF, 0x01D3, 0x01DF, 0x01E7, 0x01EB, 0x01F3, 0x01F7,
   0x01FD, 0x0209, 0x020B, 0x021D, 0x0223, 0x022D, 0x0233, 0x0239,
   0x023B, 0x0241, 0x024B, 0x0251, 0x0257, 0x0259, 0x025F, 0x0265,
   0x0269, 0x026B, 0x0277, 0x0281, 0x0283, 0x0287, 0x028D, 0x0293,
   0x0295, 0x02A1, 0x02A5, 0x02AB, 0x02B3, 0x02BD, 0x02C5, 0x02CF,

   0x02D7, 0x02DD, 0x02E3, 0x02E7, 0x02EF, 0x02F5, 0x02F9, 0x0301,
   0x0305, 0x0313, 0x031D, 0x0329, 0x032B, 0x0335, 0x0337, 0x033B,
   0x033D, 0x0347, 0x0355, 0x0359, 0x035B, 0x035F, 0x036D, 0x0371,
   0x0373, 0x0377, 0x038B, 0x038F, 0x0397, 0x03A1, 0x03A9, 0x03AD,
   0x03B3, 0x03B9, 0x03C7, 0x03CB, 0x03D1, 0x03D7, 0x03DF, 0x03E5,
   0x03F1, 0x03F5, 0x03FB, 0x03FD, 0x0407, 0x0409, 0x040F, 0x0419,
   0x041B, 0x0425, 0x0427, 0x042D, 0x043F, 0x0443, 0x0445, 0x0449,
   0x044F, 0x0455, 0x045D, 0x0463, 0x0469, 0x047F, 0x0481, 0x048B,

   0x0493, 0x049D, 0x04A3, 0x04A9, 0x04B1, 0x04BD, 0x04C1, 0x04C7,
   0x04CD, 0x04CF, 0x04D5, 0x04E1, 0x04EB, 0x04FD, 0x04FF, 0x0503,
   0x0509, 0x050B, 0x0511, 0x0515, 0x0517, 0x051B, 0x0527, 0x0529,
   0x052F, 0x0551, 0x0557, 0x055D, 0x0565, 0x0577, 0x0581, 0x058F,
   0x0593, 0x0595, 0x0599, 0x059F, 0x05A7, 0x05AB, 0x05AD, 0x05B3,
   0x05BF, 0x05C9, 0x05CB, 0x05CF, 0x05D1, 0x05D5, 0x05DB, 0x05E7,
   0x05F3, 0x05FB, 0x0607, 0x060D, 0x0611, 0x0617, 0x061F, 0x0623,
   0x062B, 0x062F, 0x063D, 0x0641, 0x0647, 0x0649, 0x064D, 0x0653
};



mp_err s_mp_prime_is_divisible(const mp_int *a, bool *result)
{
   int i;
   for (i = 0; i < MP_PRIME_TAB_SIZE; i++) {
      /* what is a mod LBL_prime_tab[i] */
      mp_err err;
      mp_digit res;
      if ((err = mp_mod_d(a, s_mp_prime_tab[i], &res)) != MP_OKAY) {
         return err;
      }

      /* is the residue zero? */
      if (res == 0u) {
         *result = true;
         return MP_OKAY;
      }
   }

   /* default to not */
   *result = false;
   return MP_OKAY;
}





/* Check if remainders are possible squares - fast exclude non-squares */
static const char rem_128[128] = {
   0, 0, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1,
   0, 0, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1,
   1, 0, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1,
   1, 0, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1,
   0, 0, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1,
   1, 0, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1,
   1, 0, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1,
   1, 0, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1
};

static const char rem_105[105] = {
   0, 0, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1,
   0, 0, 1, 1, 1, 1, 0, 1, 1, 1, 0, 1, 1, 1, 1,
   0, 1, 1, 1, 1, 1, 0, 1, 1, 0, 1, 1, 1, 1, 1,
   1, 0, 1, 1, 0, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1,
   0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1,
   1, 1, 1, 1, 0, 1, 0, 1, 1, 0, 0, 1, 1, 1, 1,
   1, 0, 1, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1
};

/* Store non-zero to ret if arg is square, and zero if not */
mp_err mp_is_square(const mp_int *arg, bool *ret)
{
   mp_err   err;
   mp_digit c;
   mp_int   t;
   uint32_t r;

   /* Default to Non-square :) */
   *ret = false;

   if (mp_isneg(arg)) {
      return MP_VAL;
   }

   if (mp_iszero(arg)) {
      return MP_OKAY;
   }

   /* First check mod 128 (suppose that MP_DIGIT_BIT is at least 7) */
   if (rem_128[127u & arg->dp[0]] == (char)1) {
      return MP_OKAY;
   }

   /* Next check mod 105 (3*5*7) */
   if ((err = mp_mod_d(arg, 105uL, &c)) != MP_OKAY) {
      return err;
   }
   if (rem_105[c] == (char)1) {
      return MP_OKAY;
   }


   if ((err = mp_init_u32(&t, 11u*13u*17u*19u*23u*29u*31u)) != MP_OKAY) {
      return err;
   }
   if ((err = mp_mod(arg, &t, &t)) != MP_OKAY) {
      goto LBL_ERR;
   }
   r = mp_get_u32(&t);
   /* Check for other prime modules, note it's not an ERROR but we must
    * free "t" so the easiest way is to goto LBL_ERR.  We know that err
    * is already equal to MP_OKAY from the mp_mod call
    */
   if (((1uL<<(r%11uL)) & 0x5C4uL) != 0uL)         goto LBL_ERR;
   if (((1uL<<(r%13uL)) & 0x9E4uL) != 0uL)         goto LBL_ERR;
   if (((1uL<<(r%17uL)) & 0x5CE8uL) != 0uL)        goto LBL_ERR;
   if (((1uL<<(r%19uL)) & 0x4F50CuL) != 0uL)       goto LBL_ERR;
   if (((1uL<<(r%23uL)) & 0x7ACCA0uL) != 0uL)      goto LBL_ERR;
   if (((1uL<<(r%29uL)) & 0xC2EDD0CuL) != 0uL)     goto LBL_ERR;
   if (((1uL<<(r%31uL)) & 0x6DE2B848uL) != 0uL)    goto LBL_ERR;

   /* Final check - is sqr(sqrt(arg)) == arg ? */
   if ((err = mp_sqrt(arg, &t)) != MP_OKAY) {
      goto LBL_ERR;
   }
   if ((err = mp_sqr(&t, &t)) != MP_OKAY) {
      goto LBL_ERR;
   }

   *ret = (mp_cmp_mag(&t, arg) == MP_EQ);
LBL_ERR:
   mp_clear(&t);
   return err;
}


/* shift left a certain amount of digits */
mp_err mp_lshd(mp_int *a, int b)
{
   mp_err err;
   int x;

   /* if its less than zero return */
   if (b <= 0) {
      return MP_OKAY;
   }
   /* no need to shift 0 around */
   if (mp_iszero(a)) {
      return MP_OKAY;
   }

   /* grow to fit the new digits */
   if ((err = mp_grow(a, a->used + b)) != MP_OKAY) {
      return err;
   }

   /* increment the used by the shift amount then copy upwards */
   a->used += b;

   /* much like mp_rshd this is implemented using a sliding window
    * except the window goes the otherway around.  Copying from
    * the bottom to the top.  see mp_rshd.c for more info.
    */
   for (x = a->used; x --> b;) {
      a->dp[x] = a->dp[x - b];
   }

   /* zero the lower digits */
   s_mp_zero_digs(a->dp, b);

   return MP_OKAY;
}




/* trim unused digits
 *
 * This is used to ensure that leading zero digits are
 * trimed and the leading "used" digit will be non-zero
 * Typically very fast.  Also fixes the sign if there
 * are no more leading digits
 */
void mp_clamp(mp_int *a)
{
   /* decrease used while the most significant digit is
    * zero.
    */
   while ((a->used > 0) && (a->dp[a->used - 1] == 0u)) {
      --(a->used);
   }

   /* reset the sign flag if zero */
   if (mp_iszero(a)) {
      a->sign = MP_ZPOS;
   }
}


/* set to zero */
void mp_zero(mp_int *a)
{
   a->sign = MP_ZPOS;
   s_mp_zero_digs(a->dp, a->used);
   a->used = 0;
}


void s_mp_copy_digs(mp_digit *d, const mp_digit *s, int digits)
{

   while (digits-- > 0) {
      *d++ = *s++;
   }

}


/* copy, b = a */
mp_err mp_copy(const mp_int *a, mp_int *b)
{
   mp_err err;

   /* if dst == src do nothing */
   if (a == b) {
      return MP_OKAY;
   }

   /* grow dest */
   if ((err = mp_grow(b, a->used)) != MP_OKAY) {
      return err;
   }

   /* copy everything over and zero high digits */
   s_mp_copy_digs(b->dp, a->dp, a->used);
   s_mp_zero_digs(b->dp + a->used, b->used - a->used);
   b->used = a->used;
   b->sign = a->sign;

   return MP_OKAY;
}

/* shift left by a certain bit count */
mp_err mp_mul_2d(const mp_int *a, int b, mp_int *c)
{
   mp_err err;

   if (b < 0) {
      return MP_VAL;
   }

   if ((err = mp_copy(a, c)) != MP_OKAY) {
      return err;
   }

   if ((err = mp_grow(c, c->used + (b / MP_DIGIT_BIT) + 1)) != MP_OKAY) {
      return err;
   }

   /* shift by as many digits in the bit count */
   if (b >= MP_DIGIT_BIT) {
      if ((err = mp_lshd(c, b / MP_DIGIT_BIT)) != MP_OKAY) {
         return err;
      }
   }

   /* shift any bit count < MP_DIGIT_BIT */
   b %= MP_DIGIT_BIT;
   if (b != 0u) {
      mp_digit shift, mask, r;
      int x;

      /* bitmask for carries */
      mask = ((mp_digit)1 << b) - (mp_digit)1;

      /* shift for msbs */
      shift = (mp_digit)(MP_DIGIT_BIT - b);

      /* carry */
      r    = 0;
      for (x = 0; x < c->used; x++) {
         /* get the higher bits of the current word */
         mp_digit rr = (c->dp[x] >> shift) & mask;

         /* shift the current word and OR in the carry */
         c->dp[x] = ((c->dp[x] << b) | r) & MP_MASK;

         /* set the carry to the carry bits of the current word */
         r = rr;
      }

      /* set final carry */
      if (r != 0u) {
         c->dp[(c->used)++] = r;
      }
   }
   mp_clamp(c);
   return MP_OKAY;
}



/* grow as required */
mp_err mp_grow(mp_int *a, int size)
{
   /* if the alloc size is smaller alloc more ram */
   if (a->alloc < size) {
      mp_digit *dp;

      if (size > MP_MAX_DIGIT_COUNT) {
         return MP_OVF;
      }

      /* reallocate the array a->dp
       *
       * We store the return in a temporary variable
       * in case the operation failed we don't want
       * to overwrite the dp member of a.
       */
      dp = (mp_digit *) MP_REALLOC(a->dp,
                                   (size_t)a->alloc * sizeof(mp_digit),
                                   (size_t)size * sizeof(mp_digit));
      if (dp == NULL) {
         /* reallocation failed but "a" is still valid [can be freed] */
         return MP_MEM;
      }

      /* reallocation succeeded so set a->dp */
      a->dp = dp;

      /* zero excess digits */
      s_mp_zero_digs(a->dp + a->alloc, size - a->alloc);
      a->alloc = size;
   }
   return MP_OKAY;
}



static mp_err s_read_urandom(void *p, size_t n)
{
   int fd;
   char *q = (char *)p;

   do {
      fd = open(MP_DEV_URANDOM, O_RDONLY);
   } while ((fd == -1) && (errno == EINTR));
   if (fd == -1) return MP_ERR;

   while (n > 0u) {
      ssize_t ret = read(fd, p, n);
      if (ret < 0) {
         if (errno == EINTR) {
            continue;
         }
         close(fd);
         return MP_ERR;
      }
      q += ret;
      n -= (size_t)ret;
   }

   close(fd);
   return MP_OKAY;
}


//
mp_err mp_rand(mp_int *a, int digits)
{
   int i;
   mp_err err;

   mp_zero(a);

   if (digits <= 0) {
      return MP_OKAY;
   }

   if ((err = mp_grow(a, digits)) != MP_OKAY) {
      return err;
   }

   if ((err = s_read_urandom(a->dp, (size_t)digits * sizeof(mp_digit))) != MP_OKAY) {
      return err;
   }

   /* TODO: We ensure that the highest digit is nonzero. Should this be removed? */
   while ((a->dp[digits - 1] & MP_MASK) == 0u) {
      if ((err = s_read_urandom(a->dp + digits - 1, sizeof(mp_digit))) != MP_OKAY) {
         return err;
      }
   }

   a->used = digits;
   for (i = 0; i < digits; ++i) {
      a->dp[i] &= MP_MASK;
   }

   return MP_OKAY;
}

/* b = a*2 */
mp_err mp_mul_2(const mp_int *a, mp_int *b)
{
   mp_err err;
   int x, oldused;
   mp_digit r;

   /* grow to accomodate result */
   if ((err = mp_grow(b, a->used + 1)) != MP_OKAY) {
      return err;
   }

   oldused = b->used;
   b->used = a->used;

   /* carry */
   r = 0;
   for (x = 0; x < a->used; x++) {

      /* get what will be the *next* carry bit from the
       * MSB of the current digit
       */
      mp_digit rr = a->dp[x] >> (mp_digit)(MP_DIGIT_BIT - 1);

      /* now shift up this digit, add in the carry [from the previous] */
      b->dp[x] = ((a->dp[x] << 1uL) | r) & MP_MASK;

      /* copy the carry that would be from the source
       * digit into the next iteration
       */
      r = rr;
   }

   /* new leading digit? */
   if (r != 0u) {
      /* add a MSB which is always 1 at this point */
      b->dp[b->used++] = 1;
   }

   /* now zero any excess digits on the destination
    * that we didn't write to
    */
   s_mp_zero_digs(b->dp + b->used, oldused - b->used);

   b->sign = a->sign;
   return MP_OKAY;
}



/* b = a/2 */
mp_err mp_div_2(const mp_int *a, mp_int *b)
{
   mp_err err;
   int x, oldused;
   mp_digit r;

   if ((err = mp_grow(b, a->used)) != MP_OKAY) {
      return err;
   }

   oldused = b->used;
   b->used = a->used;

   /* carry */
   r = 0;
   for (x = b->used; x --> 0;) {
      /* get the carry for the next iteration */
      mp_digit rr = a->dp[x] & 1u;

      /* shift the current digit, add in carry and store */
      b->dp[x] = (a->dp[x] >> 1) | (r << (MP_DIGIT_BIT - 1));

      /* forward carry to next iteration */
      r = rr;
   }

   /* zero excess digits */
   s_mp_zero_digs(b->dp + b->used, oldused - b->used);

   b->sign = a->sign;
   mp_clamp(b);
   return MP_OKAY;
}


/* single digit subtraction */
mp_err mp_sub_d(const mp_int *a, mp_digit b, mp_int *c)
{
   mp_err err;
   int oldused;

   /* fast path for a == c */
   if (a == c) {
      if ((c->sign == MP_NEG) &&
          ((c->dp[0] + b) < MP_DIGIT_MAX)) {
         c->dp[0] += b;
         return MP_OKAY;
      }
      if ((c->sign == MP_ZPOS) &&
          (c->dp[0] > b)) {
         c->dp[0] -= b;
         return MP_OKAY;
      }
   }

   /* grow c as required */
   if ((err = mp_grow(c, a->used + 1)) != MP_OKAY) {
      return err;
   }

   /* if a is negative just do an unsigned
    * addition [with fudged signs]
    */
   if (a->sign == MP_NEG) {
      mp_int a_ = *a;
      a_.sign = MP_ZPOS;
      err     = mp_add_d(&a_, b, c);
      c->sign = MP_NEG;

      /* clamp */
      mp_clamp(c);

      return err;
   }

   oldused = c->used;

   /* if a <= b simply fix the single digit */
   if (((a->used == 1) && (a->dp[0] <= b)) || mp_iszero(a)) {
      c->dp[0] = (a->used == 1) ? b - a->dp[0] : b;

      /* negative/1digit */
      c->sign = MP_NEG;
      c->used = 1;
   } else {
      int i;
      mp_digit mu = b;

      /* positive/size */
      c->sign = MP_ZPOS;
      c->used = a->used;

      /* subtract digits, mu is carry */
      for (i = 0; i < a->used; i++) {
         c->dp[i] = a->dp[i] - mu;
         mu = c->dp[i] >> (MP_SIZEOF_BITS(mp_digit) - 1u);
         c->dp[i] &= MP_MASK;
      }
   }

   /* zero excess digits */
   s_mp_zero_digs(c->dp + c->used, oldused - c->used);

   mp_clamp(c);
   return MP_OKAY;
}




/* single digit addition */
mp_err mp_add_d(const mp_int *a, mp_digit b, mp_int *c)
{
   mp_err err;
   int oldused;

   /* fast path for a == c */
   if (a == c) {
      if (!mp_isneg(c) &&
          !mp_iszero(c) &&
          ((c->dp[0] + b) < MP_DIGIT_MAX)) {
         c->dp[0] += b;
         return MP_OKAY;
      }
      if (mp_isneg(c) &&
          (c->dp[0] > b)) {
         c->dp[0] -= b;
         return MP_OKAY;
      }
   }

   /* grow c as required */
   if ((err = mp_grow(c, a->used + 1)) != MP_OKAY) {
      return err;
   }

   /* if a is negative and |a| >= b, call c = |a| - b */
   if (mp_isneg(a) && ((a->used > 1) || (a->dp[0] >= b))) {
      mp_int a_ = *a;
      /* temporarily fix sign of a */
      a_.sign = MP_ZPOS;

      /* c = |a| - b */
      err = mp_sub_d(&a_, b, c);

      /* fix sign  */
      c->sign = MP_NEG;

      /* clamp */
      mp_clamp(c);

      return err;
   }

   /* old number of used digits in c */
   oldused = c->used;

   /* if a is positive */
   if (!mp_isneg(a)) {
      /* add digits, mu is carry */
      int i;
      mp_digit mu = b;
      for (i = 0; i < a->used; i++) {
         c->dp[i] = a->dp[i] + mu;
         mu = c->dp[i] >> MP_DIGIT_BIT;
         c->dp[i] &= MP_MASK;
      }
      /* set final carry */
      c->dp[i] = mu;

      /* setup size */
      c->used = a->used + 1;
   } else {
      /* a was negative and |a| < b */
      c->used = 1;

      /* the result is a single digit */
      c->dp[0] = (a->used == 1) ? b - a->dp[0] : b;
   }

   /* sign always positive */
   c->sign = MP_ZPOS;

   /* now zero to oldused */
   s_mp_zero_digs(c->dp + c->used, oldused - c->used);
   mp_clamp(c);

   return MP_OKAY;
}


/* portable integer log of two with small footprint */
static unsigned int s_floor_ilog2(int value)
{
   unsigned int r = 0;
   while ((value >>= 1) != 0) {
      r++;
   }
   return r;
}



mp_err mp_prime_is_prime(const mp_int *a, int t, bool *result)
{
   mp_int  b;
   int     ix;
   bool    res;
   mp_err  err;

   /* default to no */
   *result = false;

   /* Some shortcuts */
   /* N > 3 */
   if (a->used == 1) {
      if ((a->dp[0] == 0u) || (a->dp[0] == 1u)) {
         *result = false;
         return MP_OKAY;
      }
      if (a->dp[0] == 2u) {
         *result = true;
         return MP_OKAY;
      }
   }

   /* N must be odd */
   if (mp_iseven(a)) {
      return MP_OKAY;
   }
   /* N is not a perfect square: floor(sqrt(N))^2 != N */
   if ((err = mp_is_square(a, &res)) != MP_OKAY) {
     // printf("\nerror is at line number 6253 function mp_is_square %d\n",err);
      return err;
   }
   if (res) {
      return MP_OKAY;
   }

   /* is the input equal to one of the primes in the table? */
   for (ix = 0; ix < MP_PRIME_TAB_SIZE; ix++) {
      if (mp_cmp_d(a, s_mp_prime_tab[ix]) == MP_EQ) {
         *result = true;
         return MP_OKAY;
      }
   }
   /* first perform trial division */
   if ((err = s_mp_prime_is_divisible(a, &res)) != MP_OKAY) {
     // printf("\nerror is at line number 6269 function s_mp_prime_is_divisible %d\n",err);
      return err;
   }

   /* return if it was trivially divisible */
   if (res) {
      return MP_OKAY;
   }

   /*
       Run the Miller-Rabin test with base 2 for the BPSW test.
    */
   if ((err = mp_init_set(&b, 2uL)) != MP_OKAY) {
     // printf("\nerror is at line number 6282 function mp_init_set %d \n",err);
      return err;
   }

   if ((err = mp_prime_miller_rabin(a, &b, &res)) != MP_OKAY) {
     //  printf("\nerror is at line number 6287 function mp_prime_miller_rabin %d \n",err);
      goto LBL_B;
   }
   if (!res) {
      goto LBL_B;
   }
   /*
      Rumours have it that Mathematica does a second M-R test with base 3.
      Other rumours have it that their strong L-S test is slightly different.
      It does not hurt, though, beside a bit of extra runtime.
   */
   b.dp[0]++;
   if ((err = mp_prime_miller_rabin(a, &b, &res)) != MP_OKAY) {
      goto LBL_B;
   }
   if (!res) {
      goto LBL_B;
   }

   /*
    * Both, the Frobenius-Underwood test and the the Lucas-Selfridge test are quite
    * slow so if speed is an issue, define LTM_USE_ONLY_MR to use M-R tests with
    * bases 2, 3 and t random bases.
    */
#ifndef LTM_USE_ONLY_MR
   if (t >= 0) {
#ifdef LTM_USE_FROBENIUS_TEST
      err = mp_prime_frobenius_underwood(a, &res);
      if ((err != MP_OKAY) && (err != MP_ITER)) {
         goto LBL_B;
      }
      if (!res) {
         goto LBL_B;
      }
#else
      if ((err = mp_prime_strong_lucas_selfridge(a, &res)) != MP_OKAY) {
         goto LBL_B;
      }
      if (!res) {
         goto LBL_B;
      }
#endif
   }
#endif

   /* run at least one Miller-Rabin test with a random base */
   if (t == 0) {
      t = 1;
   }

   /*
      Only recommended if the input range is known to be < 3317044064679887385961981

      It uses the bases necessary for a deterministic M-R test if the input is
      smaller than  3317044064679887385961981
      The caller has to check the size.
      TODO: can be made a bit finer grained but comparing is not free.
   */
   if (t < 0) {
      int p_max = 0;

      /*
          Sorenson, Jonathan; Webster, Jonathan (2015).
           "Strong Pseudoprimes to Twelve Prime Bases".
       */
      /* 0x437ae92817f9fc85b7e5 = 318665857834031151167461 */
      if ((err =   mp_read_radix(&b, "437ae92817f9fc85b7e5", 16)) != MP_OKAY) {
         goto LBL_B;
      }

      if (mp_cmp(a, &b) == MP_LT) {
         p_max = 12;
      } else {
         /* 0x2be6951adc5b22410a5fd = 3317044064679887385961981 */
         if ((err = mp_read_radix(&b, "2be6951adc5b22410a5fd", 16)) != MP_OKAY) {
            goto LBL_B;
         }

         if (mp_cmp(a, &b) == MP_LT) {
            p_max = 13;
         } else {
            err = MP_VAL;
            goto LBL_B;
         }
      }

      /* we did bases 2 and 3  already, skip them */
      for (ix = 2; ix < p_max; ix++) {
         mp_set(&b, s_mp_prime_tab[ix]);
         if ((err = mp_prime_miller_rabin(a, &b, &res)) != MP_OKAY) {
            goto LBL_B;
         }
         if (!res) {
            goto LBL_B;
         }
      }
   }
   /*
       Do "t" M-R tests with random bases between 3 and "a".
       See Fips 186.4 p. 126ff
   */
   else if (t > 0) {
      unsigned int mask;
      int size_a;

      /*
       * The mp_digit's have a defined bit-size but the size of the
       * array a.dp is a simple 'int' and this library can not assume full
       * compliance to the current C-standard (ISO/IEC 9899:2011) because
       * it gets used for small embeded processors, too. Some of those MCUs
       * have compilers that one cannot call standard compliant by any means.
       * Hence the ugly type-fiddling in the following code.
       */
      size_a = mp_count_bits(a);
      mask = (1u << s_floor_ilog2(size_a)) - 1u;
      /*
         Assuming the General Rieman hypothesis (never thought to write that in a
         comment) the upper bound can be lowered to  2*(log a)^2.
         E. Bach, "Explicit bounds for primality testing and related problems,"
         Math. Comp. 55 (1990), 355-380.

            size_a = (size_a/10) * 7;
            len = 2 * (size_a * size_a);

         E.g.: a number of size 2^2048 would be reduced to the upper limit

            floor(2048/10)*7 = 1428
            2 * 1428^2       = 4078368

         (would have been ~4030331.9962 with floats and natural log instead)
         That number is smaller than 2^28, the default bit-size of mp_digit.
      */

      /*
        How many tests, you might ask? Dana Jacobsen of Math::Prime::Util fame
        does exactly 1. In words: one. Look at the end of _GMP_is_prime() in
        Math-Prime-Util-GMP-0.50/primality.c if you do not believe it.

        The function mp_rand() goes to some length to use a cryptographically
        good PRNG. That also means that the chance to always get the same base
        in the loop is non-zero, although very low.
        If the BPSW test and/or the addtional Frobenious test have been
        performed instead of just the Miller-Rabin test with the bases 2 and 3,
        a single extra test should suffice, so such a very unlikely event
        will not do much harm.

        To preemptivly answer the dangling question: no, a witness does not
        need to be prime.
      */
      for (ix = 0; ix < t; ix++) {
         unsigned int fips_rand;
         int len;

         /* mp_rand() guarantees the first digit to be non-zero */
         if ((err = mp_rand(&b, 1)) != MP_OKAY) {
            goto LBL_B;
         }
         /*
          * Reduce digit before casting because mp_digit might be bigger than
          * an unsigned int and "mask" on the other side is most probably not.
          */
         fips_rand = (unsigned int)(b.dp[0] & (mp_digit) mask);
         if (fips_rand > (unsigned int)(INT_MAX - MP_DIGIT_BIT)) {
            len = INT_MAX / MP_DIGIT_BIT;
         } else {
            len = (((int)fips_rand + MP_DIGIT_BIT) / MP_DIGIT_BIT);
         }
         /*  Unlikely. */
         if (len < 0) {
            ix--;
            continue;
         }
         if ((err = mp_rand(&b, len)) != MP_OKAY) {
            goto LBL_B;
         }
         /*
          * That number might got too big and the witness has to be
          * smaller than "a"
          */
         len = mp_count_bits(&b);
         if (len >= size_a) {
            len = (len - size_a) + 1;
            if ((err = mp_div_2d(&b, len, &b, NULL)) != MP_OKAY) {
               goto LBL_B;
            }
         }
         /* Although the chance for b <= 3 is miniscule, try again. */
         if (mp_cmp_d(&b, 3uL) != MP_GT) {
            ix--;
            continue;
         }
         if ((err = mp_prime_miller_rabin(a, &b, &res)) != MP_OKAY) {
            goto LBL_B;
         }
         if (!res) {
            goto LBL_B;
         }
      }
   }

   /* passed the test */
   *result = true;
LBL_B:
   mp_clear(&b);
   return err;
}



/* reads a uint8_t array, assumes the msb is stored first [big endian] */
mp_err mp_from_ubin(mp_int *a, const uint8_t *buf, size_t size)
{
   mp_err err;

   /* make sure there are at least two digits */
   if ((err = mp_grow(a, 2)) != MP_OKAY) {
      return err;
   }

   /* zero the int */
   mp_zero(a);

   /* read the bytes in */
   while (size-- > 0u) {
      if ((err = mp_mul_2d(a, 8, a)) != MP_OKAY) {
         return err;
      }
      a->dp[0] |= *buf++;
      a->used += 1;
   }
   mp_clamp(a);
   return MP_OKAY;
}





///////////////////////////////////////////




//function 1 (to generate prime)
mp_err mp_prime_rand(mp_int *a, int t, int size, int flags)
{
   uint8_t *tmp, maskAND, maskOR_msb, maskOR_lsb;
   int bsize, maskOR_msb_offset;
   bool res;
   mp_err err;

   /* sanity check the input */
   if ((size <= 1) || (t <= 0)) {
      return MP_VAL;
   }

   /* MP_PRIME_SAFE implies MP_PRIME_BBS */
   if ((flags & MP_PRIME_SAFE) != 0) {
      flags |= MP_PRIME_BBS;
   }

   /* calc the byte size */
   bsize = (size>>3) + ((size&7)?1:0);
  // printf("printing size in bytes %d",bsize);

   /* we need a buffer of bsize bytes */
   tmp = (uint8_t *) MP_MALLOC((size_t)bsize);
   if (tmp == NULL) {
      return MP_MEM;
   }

   /* calc the maskAND value for the MSbyte*/
   maskAND = ((size&7) == 0) ? 0xFFu : (uint8_t)(0xFFu >> (8 - (size & 7)));

   /* calc the maskOR_msb */
   maskOR_msb        = 0;
   maskOR_msb_offset = ((size & 7) == 1) ? 1 : 0;

   if ((flags & MP_PRIME_2MSB_ON) != 0) {
      maskOR_msb       |= (uint8_t)(0x80 >> ((9 - size) & 7));
   }

   /* get the maskOR_lsb */
   maskOR_lsb         = 1u;
   if ((flags & MP_PRIME_BBS) != 0) {
      maskOR_lsb     |= 3u;
   }

   do {
      /* read the bytes */
      if ((err = s_read_urandom(tmp, (size_t)bsize)) != MP_OKAY) {
        // printf("\nproblem lies in generating random number\n");
         goto LBL_ERR;
      }
     // printf("\n%s\n",tmp);
      /* work over the MSbyte */
      tmp[0]    &= maskAND;
      tmp[0]    |= (uint8_t)(1 << ((size - 1) & 7));

      /* mix in the maskORs */
      tmp[maskOR_msb_offset]   |= maskOR_msb;
      tmp[bsize-1]             |= maskOR_lsb;

      /* read it in */
      /* TODO: casting only for now until all lengths have been changed to the type "size_t"*/
      if ((err = mp_from_ubin(a, tmp, (size_t)bsize)) != MP_OKAY) {
        // printf("\nthe problem is in mp_from_ubin\n");
         goto LBL_ERR;
      }

      /* is it prime? */
      if ((err = mp_prime_is_prime(a, t, &res)) != MP_OKAY) {
        // printf("\nthe error is in mp_prime_is_prime function line 6603  %d\n",err);
         goto LBL_ERR;
      }
      if (!res) {
         continue;
      }

      if ((flags & MP_PRIME_SAFE) != 0) {
         /* see if (a-1)/2 is prime */
         if ((err = mp_sub_d(a, 1uL, a)) != MP_OKAY) {
            goto LBL_ERR;
         }
         if ((err = mp_div_2(a, a)) != MP_OKAY) {
            goto LBL_ERR;
         }

         /* is it prime? */
         if ((err = mp_prime_is_prime(a, t, &res)) != MP_OKAY) {
            goto LBL_ERR;
         }
      }
   } while (!res);

   if ((flags & MP_PRIME_SAFE) != 0) {
      /* restore a to the original value */
      if ((err = mp_mul_2(a, a)) != MP_OKAY) {
         goto LBL_ERR;
      }
      if ((err = mp_add_d(a, 1uL, a)) != MP_OKAY) {
         goto LBL_ERR;
      }
   }

   err = MP_OKAY;
LBL_ERR:
   MP_FREE_BUF(tmp, (size_t)bsize);
   return err;
}



/*
Class for SHA256 implementation
*/
///////////////////////////////////////////

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

SHA256::SHA256(): m_blocklen(0), m_bitlen(0) {
	m_state[0] = 0x6a09e667;
	m_state[1] = 0xbb67ae85;
	m_state[2] = 0x3c6ef372;
	m_state[3] = 0xa54ff53a;
	m_state[4] = 0x510e527f;
	m_state[5] = 0x9b05688c;
	m_state[6] = 0x1f83d9ab;
	m_state[7] = 0x5be0cd19;
}
void SHA256::update( uint8_t * data, size_t length) {
	for (size_t i = 0 ; i < length ; i++) {
		m_data[m_blocklen++] = data[i];
	   // printf("%c",data[i]);
		if (m_blocklen == 64) {
			transform();

			// End of the block
			m_bitlen += 512;
			m_blocklen = 0;
		}
	}
}
unsigned char digest_result[32];
uint8_t * SHA256::digest() {
	uint8_t * hash = new uint8_t[32];
   // printf("\n\ninside digest\n\n ");
	pad();
   // printf("\n\npad pass\n\n ");
	revert(hash);

   for (int i=0;i<32;i++){
      digest_result[i]=*(hash+i);
   }

   delete hash;

   return &digest_result[0];
	//return hash;
}

uint32_t SHA256::rotr(uint32_t x, uint32_t n) {
	return (x >> n) | (x << (32 - n));
}

uint32_t SHA256::choose(uint32_t e, uint32_t f, uint32_t g) {
	return (e & f) ^ (~e & g);
}

uint32_t SHA256::majority(uint32_t a, uint32_t b, uint32_t c) {
	return (a & (b | c)) | (b & c);
}

uint32_t SHA256::sig0(uint32_t x) {
	return SHA256::rotr(x, 7) ^ SHA256::rotr(x, 18) ^ (x >> 3);
}

uint32_t SHA256::sig1(uint32_t x) {
	return SHA256::rotr(x, 17) ^ SHA256::rotr(x, 19) ^ (x >> 10);
}

void SHA256::transform() {
	uint32_t maj, xorA, ch, xorE, sum, newA, newE, m[64];
	uint32_t state[8];

	for (uint8_t i = 0, j = 0; i < 16; i++, j += 4) { // Split data in 32 bit blocks for the 16 first words
		m[i] = (m_data[j] << 24) | (m_data[j + 1] << 16) | (m_data[j + 2] << 8) | (m_data[j + 3]);
	}

	for (uint8_t k = 16 ; k < 64; k++) { // Remaining 48 blocks
		m[k] = SHA256::sig1(m[k - 2]) + m[k - 7] + SHA256::sig0(m[k - 15]) + m[k - 16];
	}

	for(uint8_t i = 0 ; i < 8 ; i++) {
		state[i] = m_state[i];
	}

	for (uint8_t i = 0; i < 64; i++) {
		maj   = SHA256::majority(state[0], state[1], state[2]);
		xorA  = SHA256::rotr(state[0], 2) ^ SHA256::rotr(state[0], 13) ^ SHA256::rotr(state[0], 22);

		ch = choose(state[4], state[5], state[6]);

		xorE  = SHA256::rotr(state[4], 6) ^ SHA256::rotr(state[4], 11) ^ SHA256::rotr(state[4], 25);

		sum  = m[i] + K[i] + state[7] + ch + xorE;
		newA = xorA + maj + sum;
		newE = state[3] + sum;

		state[7] = state[6];
		state[6] = state[5];
		state[5] = state[4];
		state[4] = newE;
		state[3] = state[2];
		state[2] = state[1];
		state[1] = state[0];
		state[0] = newA;
	}

	for(uint8_t i = 0 ; i < 8 ; i++) {
		m_state[i] += state[i];
	}
}

void SHA256::pad() {

	uint64_t i = m_blocklen;
	uint8_t end = m_blocklen < 56 ? 56 : 64;

	m_data[i++] = 0x80; // Append a bit 1
	while (i < end) {
		m_data[i++] = 0x00; // Pad with zeros
	}

	if(m_blocklen >= 56) {
		transform();
		//memset(m_data, 0, 56);
		for(int g=0;g<56;g++){
			m_data[g]=0;
		}
	}

	// Append to the padding the total message's length in bits and transform.
//	printf("   %lu   ",m_bitlen);
	m_bitlen += m_blocklen * 8;
	//printf("   %lu   ",m_bitlen);
	m_data[63] = m_bitlen;
	m_data[62] = m_bitlen >> 8;
	m_data[61] = m_bitlen >> 16;
	m_data[60] = m_bitlen >> 24;
	m_data[59] = m_bitlen >> 32;
	m_data[58] = m_bitlen >> 40;
	m_data[57] = m_bitlen >> 48;
	m_data[56] = m_bitlen >> 56;




	transform();
}

void SHA256::revert(uint8_t * hash) {
	// SHA uses big endian byte ordering
	// Revert all bytes

   // printf("hello");

	for (uint8_t i = 0 ; i < 4 ; i++) {
		for(uint8_t j = 0 ; j < 8 ; j++) {
			hash[i + (j * 4)] = (m_state[j] >> (24 - i * 8)) & 0x000000ff;
            //printf("%d",(i+(j*4)));
		}
	}
}



///////////////////////////////////////////



