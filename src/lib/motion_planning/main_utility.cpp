#include "main_utility.hpp"
#include <stdio.h>
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
   mp_int A1, A2, B1;

   mp_int Q1, Q0, R1, R0, t;
   mp_int  B_0;
   int m = a->used - b->used, k = m/2;

   if (m < (MP_MUL_KARATSUBA_CUTOFF)) {
      return s_mp_div_school(a, b, q, r);
   }

   if ((err = mp_init_multi(&A1, &A2, &B1, &B_0, &Q1, &Q0, &R1, &R0, &t, NULL)) != MP_OKAY) {
      goto LBL_ERR;
   }

   /* B1 = b / beta^k, B_0 = b % beta^k*/
   if ((err = mp_div_2d(b, k * MP_DIGIT_BIT, &B1, &B_0)) != MP_OKAY)        goto LBL_ERR;

   /* (Q1, R1) =  RecursiveDivRem(A / beta^(2k), B1) */
   if ((err = mp_div_2d(a, 2*k * MP_DIGIT_BIT, &A1, &t)) != MP_OKAY)       goto LBL_ERR;
   if ((err = s_recursion(&A1, &B1, &Q1, &R1)) != MP_OKAY)                 goto LBL_ERR;

   /* A1 = (R1 * beta^(2k)) + (A % beta^(2k)) - (Q1 * B_0 * beta^k) */
   if ((err = mp_lshd(&R1, 2*k)) != MP_OKAY)                               goto LBL_ERR;
   if ((err = mp_add(&R1, &t, &A1)) != MP_OKAY)                            goto LBL_ERR;
   if ((err = mp_mul(&Q1, &B_0, &t)) != MP_OKAY)                            goto LBL_ERR;
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

   /* A2 = (R0*beta^k) +  (A1 % beta^k) - (Q0*B_0) */
   if ((err = mp_lshd(&R0, k)) != MP_OKAY)                                 goto LBL_ERR;
   if ((err = mp_add(&R0, &t, &A2)) != MP_OKAY)                            goto LBL_ERR;
   if ((err = mp_mul(&Q0, &B_0, &t)) != MP_OKAY)                            goto LBL_ERR;
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
   mp_clear_multi(&A1, &A2, &B1, &B_0, &Q1, &Q0, &R1, &R0, &t, NULL);
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

   static const signed char table[] = {0, 1, 0, -1, 0, -1, 0, 1};

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
  // return s_mp_exptmod_fast(G, X, P, Y, dr);
   return s_mp_exptmod(G, X, P, Y, 0);
   /* if the modulus is odd or dr != 0 use the montgomery method */
   if (MP_HAS(S_MP_EXPTMOD_FAST) && (mp_isodd(P) || (dr != 0))) {
      return s_mp_exptmod_fast(G, X, P, Y, dr);
   }
  // return MP_OKAY;//////////////////////////////debug
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
//overriding SHA256::update  for our custom made files.
void SHA256::update( FILE* data, int length) {
   unsigned char ch;
   FILE *fptrdef;
   fptrdef=fopen("/fs/microsd/debug_sha.txt","w");
	for (int i = 0 ; i < length ; i++)
   {
      ch=fgetc(data);
      if(ch!='\n')
      {   fprintf(fptrdef,"%c",ch);
         m_data[m_blocklen++] = ch;//data[i];
         // printf("%c",data[i]);
         if (m_blocklen == 64) {
            transform();

            // End of the block
            m_bitlen += 512;
            m_blocklen = 0;
         }

      }
	}
   fclose(fptrdef);
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

   delete [] hash;
   hash=NULL;

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




////////////////////// file_operations.h



int isSubstring(char *s1, char *s2,FILE *fptr)///s1 is the sub string ; s2 is the larger string
{
   int M = strlen(s1);
 //  printf("11 %s jkl\n",s1);

   int N=0;
   int j;
   if(fptr==NULL){
      N= strlen(s2);
      for (int i = 0; i <= N - M; i++)
      {


         for (j = 0; j < M; j++)
            if (s2[i + j] != s1[j])
               break;
         if (j == M)
            return i;
      }
   }else{
      int index_count=0;
      signed char ch;
      j=0;
      int match_count=0;
      while((ch=fgetc(fptr))!=EOF){

         if(s1[j]!=ch){
            match_count=0;
            j=0;
         }else{
            match_count=match_count+1;
            if(match_count==M){
                return (index_count-M+1);
            }
            j++;
         }
         index_count++;
      }
   }  // printf("22 %s jkl\n",s2);


   return -1;
}


int isSubstring(char *s1, char *s2)///s1 is the sub string ; s2 is the larger string
{
   int M = strlen(s1);
 //  printf("11 %s jkl\n",s1);

   int N=0;
   int j;

      N= strlen(s2);
      for (int i = 0; i <= N - M; i++)
      {


         for (j = 0; j < M; j++)
            if (s2[i + j] != s1[j])
               break;
         if (j == M)
            return i;
      }
   

   return -1;
}



int find_int_hex(char ptr){
    char char_set16[]="0123456789ABCDEF";
    int i=0;
    while(1){
        if(char_set16[i]==ptr){
            return i;
        }
        i++;
    }

}
void base64Encoder(char input_str[], int len_str, char *result)
{

    //character set of base 64 encoding scheme
    char char_set[]="ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";


  //  char *res_str=(char *) malloc(1000 * sizeof(char));//chr *res_str =new char(1000*1)

    if(len_str%2!=0){
        printf("the hex string is not valid for having conversion to base64 digits");
        char aux[1000]="000";
        strcat(aux,input_str);
        strcpy(input_str,aux);
       // exit(1);
    }

    int i=len_str;

    /* basic agenda
    an initial check that number of char in string are even or not (gouping of bytes(8 bits= 2 hex char) is done)
    first check if 6 continous hex digits could be taken or not
    if yes:
    1)  Take 6 hex digits :24 bits
    2)  24/6=4 then write 4 base64 characters
    if no:
    1)check how many hex char are available(<6)
    2)possible cases 2,4
    if 2 hex char available:
        -add 4 bits(valued 0) to make the bit length divisible by 6
    if 4 hex char available
        - add 2 bits(valued 0) to make the bit length divisible by 6
    */
    int j=0;
    int res_count=0;
    while(1){
        if(i>=6){//if yes code starts here

            int aux=0;
           // printf("%d\n",res_count);
            for(int index=j ; index<j+6;index++){
                aux=(aux<<4);
              //  printf("%d ",find_int_hex(input_str[index]));
                aux=aux|find_int_hex(input_str[index]);
              //  printf(" %08x ",aux);
            }
            j=j+6;
           // printf("\n");
            int aux1=0;
            for(int o=3;o>=0;o--){
                aux1=aux & 0x3f;
                result[res_count+o]=char_set[aux1];
                aux=aux>>6;
                //res_count--;
            }
            res_count=res_count+4;
            i=i-6;
        }else{

                if(i==0){
               //     printf("Done...");
                 //   printf("\n");
                    break;
                }
                if(i==2){
                    int aux=0;
                    for(int index=j ; index<j+2;index++){
                            aux=(aux<<4);
                        //  printf("%d ",find_int_hex(input_str[index]));
                            aux=aux|find_int_hex(input_str[index]);
                            //  printf(" %08x ",aux);
                        }
                        aux=aux<<4;
                    int aux1=0;
                    for(int o=1;o>=0;o--){
                        aux1=aux & 0x3f;
                        result[res_count+o]=char_set[aux1];
                        aux=aux>>6;
                    //res_count--;
                    }
                    res_count=res_count+2;
                    result[res_count]='=';
                    res_count++;
                    result[res_count]='=';
                    res_count++;
                    result[res_count]='\0';
                    break;

                }else{

                    int aux=0;
                    for(int index=j ; index<j+4;index++){
                            aux=(aux<<4);
                        //  printf("%d ",find_int_hex(input_str[index]));
                            aux=aux|find_int_hex(input_str[index]);
                            //  printf(" %08x ",aux);
                        }
                        aux=aux<<2;
                    int aux1=0;
                    for(int o=2;o>=0;o--){
                        aux1=aux & 0x3f;
                        result[res_count+o]=char_set[aux1];
                        aux=aux>>6;
                    //res_count--;
                    }
                    res_count=res_count+3;
                    result[res_count]='=';
                    res_count++;
                    result[res_count]='\0';
                    break;


                }

        }


    }


}





////////

// signs the content with signKey(a key) and puts the result inside result
void signing_support_0(key signKey,char *content,char *result,char *result1){
    printf("\n\n%s\n\n",content);
    char ch;
    unsigned char *st=new unsigned char[3000];
    unsigned int i_sha=0;
    //// assigning char to unsigned char, this has to be done to implement sha256
    int i_sha_a=0;
    while(1){

	   ch=content[i_sha_a];
      // printf("%c",ch);
	   if(ch=='\0'){
		   break;
	   }
       if(ch!='\n'){
        *(st+i_sha)=ch;
	//	printf("%c",*(st+i_sha));
		i_sha++;}
        i_sha_a++;

    }
    *(st+i_sha)='\0';
    printf("ooooLLL\n\n%s\n\n",st);

    SHA256 sha;

	sha.update(st,i_sha);

	uint8_t * digest = sha.digest();
    //int it=0;//just an iterator

    char HEX_format_Digest[65];

    // code for making hex string
    int hex_count=0;
    for(int yui=0;yui<32;yui++){
     //    printf("%d ",*(digest+yui)>>4||0x0f);
         HEX_format_Digest[hex_count]=HEX_ele[*(digest+yui)>>4 & 0x0f];
         hex_count++;
         HEX_format_Digest[hex_count]=HEX_ele[*(digest+yui) & 0x0f];
         hex_count++;
      }
   // delete[] digest;
    //digest=NULL;
    HEX_format_Digest[hex_count]='\0';
static int pass=0;
pass++;
    printf("\nhere is the string %d:%s\n\n",pass,HEX_format_Digest);
    strcpy(result1,HEX_format_Digest);

    //At this point we have got the string format for the digest


    delete []st;
    st=NULL;
    ////////////////////

    // modulus in public key

    mp_int modulus;
    mp_init(&modulus);
    mp_read_radix(&modulus,signKey.modulus,10);

    //private key exponent
    mp_int Private_key;
    mp_init(&Private_key);
    mp_read_radix(&Private_key,signKey.private_exponent,10);

   //making the signing process pkcs.1.15 compatible
   char padding_SHA256[446]="1ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff003031300d060960864801650304020105000420";
   strcat(padding_SHA256,HEX_format_Digest);
   //printf("\npadded sha digest :%s\n",padding_SHA256);
   mp_int hash_to_sign;
   mp_init(&hash_to_sign);
   mp_read_radix(&hash_to_sign,padding_SHA256,16);
   char aux_hash_ex[514];
   mp_to_hex(&hash_to_sign,aux_hash_ex,sizeof(aux_hash_ex));
   printf("\n   hash to sign ===\n%s\n\n",aux_hash_ex);

   ////Signing begins

   mp_int cipher;
   mp_init(&cipher);

   mp_exptmod(&hash_to_sign,&Private_key,&modulus,&cipher);

   char cipher_string_hex[513];
   mp_to_hex(&cipher,cipher_string_hex,sizeof(cipher_string_hex));
   printf("\n%s\n",cipher_string_hex);//this is the encrypted message


    base64Encoder(cipher_string_hex,strlen(cipher_string_hex),result);


mp_clear(&hash_to_sign);
mp_clear(&modulus);
mp_clear(&cipher);
mp_clear(&Private_key);

}

//////////









// signs the content with signKey(a key) and puts the result inside result
void signing_support(key signKey,char *content,char*result){
    printf("\n\n%s\n\n",content);
    char ch;
    unsigned char *st=new unsigned char[3000];
    unsigned int i_sha=0;
    //// assigning char to unsigned char, this has to be done to implement sha256
    int i_sha_a=0;
    while(1){

	   ch=content[i_sha_a];
      // printf("%c",ch);
	   if(ch=='\0'){
		   break;
	   }
       if(ch!='\n'){
        *(st+i_sha)=ch;
	//	printf("%c",*(st+i_sha));
		i_sha++;}
        i_sha_a++;

    }
    *(st+i_sha)='\0';
    printf("ooooLLL\n\n%s\n\n",st);

    SHA256 sha;

	sha.update(st,i_sha);

	uint8_t * digest = sha.digest();
    //int it=0;//just an iterator

    char HEX_format_Digest[65];

    // code for making hex string
    int hex_count=0;
    for(int yui=0;yui<32;yui++){
     //    printf("%d ",*(digest+yui)>>4||0x0f);
         HEX_format_Digest[hex_count]=HEX_ele[*(digest+yui)>>4 & 0x0f];
         hex_count++;
         HEX_format_Digest[hex_count]=HEX_ele[*(digest+yui) & 0x0f];
         hex_count++;
      }
   // delete[] digest;
    //digest=NULL;
    HEX_format_Digest[hex_count]='\0';
static int pass=0;
pass++;
    printf("\nhere is the string %d:%s\n\n",pass,HEX_format_Digest);

    //At this point we have got the string format for the digest


    delete []st;
    st=NULL;
    ////////////////////

    // modulus in public key

    mp_int modulus;
    mp_init(&modulus);
    mp_read_radix(&modulus,signKey.modulus,10);

    //private key exponent
    mp_int Private_key;
    mp_init(&Private_key);
    mp_read_radix(&Private_key,signKey.private_exponent,10);

   //making the signing process pkcs.1.15 compatible
   char padding_SHA256[446]="1ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff003031300d060960864801650304020105000420";
   strcat(padding_SHA256,HEX_format_Digest);
   //printf("\npadded sha digest :%s\n",padding_SHA256);
   mp_int hash_to_sign;
   mp_init(&hash_to_sign);
   mp_read_radix(&hash_to_sign,padding_SHA256,16);
   char aux_hash_ex[514];
   mp_to_hex(&hash_to_sign,aux_hash_ex,sizeof(aux_hash_ex));
   printf("\n   hash to sign ===\n%s\n\n",aux_hash_ex);

   ////Signing begins

   mp_int cipher;
   mp_init(&cipher);

   mp_exptmod(&hash_to_sign,&Private_key,&modulus,&cipher);

   char cipher_string_hex[513];
   mp_to_hex(&cipher,cipher_string_hex,sizeof(cipher_string_hex));
   printf("\n%s\n",cipher_string_hex);//this is the encrypted message


    base64Encoder(cipher_string_hex,strlen(cipher_string_hex),result);


mp_clear(&hash_to_sign);
mp_clear(&modulus);
mp_clear(&cipher);
mp_clear(&Private_key);

}

void pair_file_write(pair_set *ptr,int ptr_quant ,char *file_name , key Skey){


    char content[6000];//=(char*) malloc(sizeof(char)*6000);
    strcpy(content,"<content>\n");

    int i=0;
    while(i<ptr_quant){

    strcat(content,"<");
    strcat(content,ptr[i].tag);
    strcat(content,"=");
    strcat(content,ptr[i].value);
    strcat(content,">\n");
    i++;
    }
    strcat(content,"</content>\n");

    strcat(content,"\0");

    printf("\nINside pair_file_write %s \n",content);

    // now signing takes place
    char signature[2000];//=(char*) malloc(sizeof(char)*2000);
    signing_support( Skey,content,signature);
    content[int(strlen(content))-1]='\n';
    strcat(content,"<Sign>\n");
    strcat(content,"<Signature=");
    strcat(content,signature);
    strcat(content,">\n");
    strcat(content,"</Sign>");
    strcat(content,"\0");




    ////////////////////

    FILE *fptr1;
    fptr1=fopen(file_name,"w");
    fprintf(fptr1, "%s", content);
    fclose(fptr1);
    //free(content);
  //  free(signature);
}

//this function is for validating the txt file wrt  the given key

int inBase64(char *d){
    for(int i=0;i<64;i++){
        if(*d==BASE64_DIGITS[i]){
            return i;
        }
    }
    return 0;
}


void base64decoder(char *base64_ptr,char *hex){



    int ik=0;



    int aux1,aux2,aux3,num_bits;

    for(int i=0;base64_ptr[i]!='\0';i=i+4){
        aux1=0;num_bits=0;

        for(int j=0;j<4;j++){

        if(base64_ptr[i+j]!='='){
        aux1=aux1<<6;
        num_bits=num_bits+6;

        aux2=inBase64(&base64_ptr[i+j]);


        aux1=aux1|aux2;


        }else{
            aux1=aux1>>2;
            num_bits=num_bits-2;
        }

        }

      //  int count_bit=num_bits;
        while(num_bits!=0){
            num_bits=num_bits-4;
            aux3=(aux1>>num_bits) & 15;


           // printf("aux3:: %d ",aux3);
            hex[ik]=HEX_DIGITS[aux3];

          // printf("%c |",hex[ik]);
            ik++;
        }

       // printf("\n");


    }
    hex[ik]='\0';

   // printf("length %d\n",ik);
   // printf("\n");
    /*int i=0;

    while(hex[i]!='\0'){
        printf("%c",hex[i]);
        i++;
    }*/

   // printf("\n");

}

char small_letter(char a){
   char uset[]="0123456789";
   for(int i=0;i<10;i++){
      if(uset[i]==a){
         return a;
      }
   }
   return a+32;
}

int Validating_File(char *file , key key, FILE *fptr_de){

    char tag1[30]="</content>";
    char tag2[30]="Signature";
    printf("%s",tag2);
    FILE *fptr;
    fptr=fopen(file,"r");
/*    signed char ch;
    char file_content[3000];//=(char*) malloc(sizeof(char)*3000);
    int i=0;

    while((ch=fgetc(fptr))!=EOF)
    {
        file_content[i]=ch;
        i++;
    }

    file_content[i]='\0';
   fclose(fptr);*/

    int index=isSubstring(tag1,NULL,fptr);
    fclose(fptr);
    fptr=fopen(file,"r");

    index=index+int(strlen(tag1));


/// sha256 implemenrtion begins input should be a file pointer(this only for custom made files, for pa, scenario is different
//   second input should be  (index+1)
//  output should be 64 long hash.
    /*unsigned char *st=new unsigned char[2000];
    unsigned int i_sha=0;
    //// assigning char to unsigned char, this has to be done to implement sha256
    int i_sha_a=0;
    while(i_sha_a<index+1)
    {

	   ch=file_content[i_sha_a];
      // printf("%c",ch);

       if(ch!='\n'){
        *(st+i_sha)=ch;
	//	printf("%c",*(st+i_sha));
		i_sha++;}
        i_sha_a++;

    }



   SHA256 sha;

	sha.update(st,i_sha);
   delete [] st;
   st=NULL;

	uint8_t * digest = sha.digest();
   // int it=0;//just an iterator

    char HEX_format_Digest[65];

    // code for making hex string
    int hex_count=0;
    for(int yui=0;yui<32;yui++)
    {
        // printf("%d ",*(digest+yui)>>4||0x0f);
         HEX_format_Digest[hex_count]=HEX_ele[*(digest+yui)>>4 & 0x0f];
         hex_count++;
         HEX_format_Digest[hex_count]=HEX_ele[*(digest+yui) & 0x0f];
         hex_count++;
    }
    HEX_format_Digest[hex_count]='\0';
    */
    char HEX_format_Digest[65];

    Sha256_implement(NULL,HEX_format_Digest,fptr,index+1);// implemented using file pointer
    fclose(fptr);
    fprintf(fptr_de,"\n\nHEX format %s\n\n",HEX_format_Digest);
    
    static int pass=0;
    pass++;
    printf("\nhere is the     string %d: %s\n\n",pass,HEX_format_Digest);
   // delete []st;
    //st=NULL;

   char signature[520];//=(char*) malloc(sizeof(char)*500);
   fetch_tag(NULL,tag2,signature,file);
   
   fprintf(fptr_de,"\n\nfetched signature %s\n\n",signature);
   // free(file_content);
    printf("signature\n%s\n",signature);
    // base64 decoder


   char hex[1000];//=(char*) malloc(sizeof(char)*(1000));

   base64decoder(&signature[0],hex);

   char decrypted_hex[1000];//=(char*) malloc(sizeof(char)*1000);

 //  check_debug_mp(hex,decrypted_hex, &key);

   // printf("\n%s\n",hex);
   //mp_expmod function input 1)hex  output 1)useful_decrypted_hex
    mp_int HEX1;
    mp_init(&HEX1);
    //mp_read_radix(&HEX1,"6676CB59FC89868EB6F2EF269CEF076E265C963779DE44B9E2E234A3391043B10E7667892400753214A9B1FD51AB7F48A429BD6AE73B0EC894785CCE3E0EFD735C4BBD54D2B9F7709629BC6C5A635F1AF52BBEBA1352D876154EDFA8EF4F1C58D4EFF9ADAB3EB81AF329B35595BA94B98505B67EBB814963B71C35312CA2904BA56CC2A4DDBD53D161BB900A74B5CF647531476A343293895433F70A0A35E7110EC220299A9F685BF6A98685925C3DA603BFC11EE0BF6E2216F47873DEF58EDB0CFB4CAE158F70E60E6233B09542CAA1F21722CAC5F24A8C09E4A32AFD34B879C8CA68E0DBD4CBE65F30A793333D5983B006EB91CC5FD86549939D526EBB3CE6",16);
    mp_read_radix(&HEX1,hex,16);
    mp_int modulus;
    mp_init(&modulus);
    mp_read_radix(&modulus,key.modulus,10);
    mp_int public_key;
    mp_init(&public_key);
    mp_read_radix(&public_key,"65537",10);
    mp_int decrypted;
    mp_init(&decrypted);
    mp_exptmod(&HEX1,&public_key,&modulus,&decrypted);
    mp_to_hex(&decrypted,decrypted_hex,sizeof(decrypted_hex));

    printf("\n%s\n",decrypted_hex);//this is the encrypted message

    char useful_decrypted_hex[100];//=(char*) malloc(sizeof(char)*400);

    int i=445;
    int i_counter=0;
    while(i<int(strlen(decrypted_hex))){
        useful_decrypted_hex[i_counter]=small_letter(decrypted_hex[i]);
        i++;
        i_counter++;
    }
    useful_decrypted_hex[i_counter]='\0';

  //  printf("heloo\n\n%s\n\n %s\n\n",useful_decrypted_hex,HEX_format_Digest);

  //  free(decrypted_hex);
   // free(useful_decrypted_hex);
  // mp_clear(&HEX1);
   //mp_clear(&modulus);
   //mp_clear(&public_key);
   //mp_clear(&decrypted);
fprintf(fptr_de,"\n\nHEX format %s\n\n",HEX_format_Digest);
fprintf(fptr_de,"\n\n from mp_exptmod  %s\n\n",useful_decrypted_hex);
   if (strcmp(useful_decrypted_hex,HEX_format_Digest)==0){
       return 1;// valid
   }else{

    return 0;// not valid
   }


}




int Validating_File(char *file , key key){

    char tag1[30]="</content>";
    char tag2[30]="Signature";
    printf("%s",tag2);
    FILE *fptr;
    fptr=fopen(file,"r");
/*    signed char ch;
    char file_content[3000];//=(char*) malloc(sizeof(char)*3000);
    int i=0;

    while((ch=fgetc(fptr))!=EOF)
    {
        file_content[i]=ch;
        i++;
    }

    file_content[i]='\0';
   fclose(fptr);*/

    int index=isSubstring(tag1,NULL,fptr);
    fclose(fptr);
    fptr=fopen(file,"r");

    index=index+int(strlen(tag1));


/// sha256 implemenrtion begins input should be a file pointer(this only for custom made files, for pa, scenario is different
//   second input should be  (index+1)
//  output should be 64 long hash.
    /*unsigned char *st=new unsigned char[2000];
    unsigned int i_sha=0;
    //// assigning char to unsigned char, this has to be done to implement sha256
    int i_sha_a=0;
    while(i_sha_a<index+1)
    {

	   ch=file_content[i_sha_a];
      // printf("%c",ch);

       if(ch!='\n'){
        *(st+i_sha)=ch;
	//	printf("%c",*(st+i_sha));
		i_sha++;}
        i_sha_a++;

    }



   SHA256 sha;

	sha.update(st,i_sha);
   delete [] st;
   st=NULL;

	uint8_t * digest = sha.digest();
   // int it=0;//just an iterator

    char HEX_format_Digest[65];

    // code for making hex string
    int hex_count=0;
    for(int yui=0;yui<32;yui++)
    {
        // printf("%d ",*(digest+yui)>>4||0x0f);
         HEX_format_Digest[hex_count]=HEX_ele[*(digest+yui)>>4 & 0x0f];
         hex_count++;
         HEX_format_Digest[hex_count]=HEX_ele[*(digest+yui) & 0x0f];
         hex_count++;
    }
    HEX_format_Digest[hex_count]='\0';
    */
    char HEX_format_Digest[65];

    Sha256_implement(NULL,HEX_format_Digest,fptr,index+1);// implemented using file pointer
    fclose(fptr);

    
    static int pass=0;
    pass++;
    printf("\nhere is the     string %d: %s\n\n",pass,HEX_format_Digest);
   // delete []st;
    //st=NULL;

   char signature[500];//=(char*) malloc(sizeof(char)*500);
   fetch_tag(NULL,tag2,signature,file);
   
   // free(file_content);
    printf("signature\n%s\n",signature);
    // base64 decoder


   char hex[1000];//=(char*) malloc(sizeof(char)*(1000));

   base64decoder(&signature[0],hex);

   char decrypted_hex[1000];//=(char*) malloc(sizeof(char)*1000);

 //  check_debug_mp(hex,decrypted_hex, &key);

   // printf("\n%s\n",hex);
   //mp_expmod function input 1)hex  output 1)useful_decrypted_hex
    mp_int HEX1;
    mp_init(&HEX1);
    //mp_read_radix(&HEX1,"6676CB59FC89868EB6F2EF269CEF076E265C963779DE44B9E2E234A3391043B10E7667892400753214A9B1FD51AB7F48A429BD6AE73B0EC894785CCE3E0EFD735C4BBD54D2B9F7709629BC6C5A635F1AF52BBEBA1352D876154EDFA8EF4F1C58D4EFF9ADAB3EB81AF329B35595BA94B98505B67EBB814963B71C35312CA2904BA56CC2A4DDBD53D161BB900A74B5CF647531476A343293895433F70A0A35E7110EC220299A9F685BF6A98685925C3DA603BFC11EE0BF6E2216F47873DEF58EDB0CFB4CAE158F70E60E6233B09542CAA1F21722CAC5F24A8C09E4A32AFD34B879C8CA68E0DBD4CBE65F30A793333D5983B006EB91CC5FD86549939D526EBB3CE6",16);
    mp_read_radix(&HEX1,hex,16);
    mp_int modulus;
    mp_init(&modulus);
    mp_read_radix(&modulus,key.modulus,10);
    mp_int public_key;
    mp_init(&public_key);
    mp_read_radix(&public_key,"65537",10);
    mp_int decrypted;
    mp_init(&decrypted);
    mp_exptmod(&HEX1,&public_key,&modulus,&decrypted);
    mp_to_hex(&decrypted,decrypted_hex,sizeof(decrypted_hex));

    printf("\n%s\n",decrypted_hex);//this is the encrypted message

    char useful_decrypted_hex[100];//=(char*) malloc(sizeof(char)*400);

    int i=445;
    int i_counter=0;
    while(i<int(strlen(decrypted_hex))){
        useful_decrypted_hex[i_counter]=small_letter(decrypted_hex[i]);
        i++;
        i_counter++;
    }
    useful_decrypted_hex[i_counter]='\0';

  //  printf("heloo\n\n%s\n\n %s\n\n",useful_decrypted_hex,HEX_format_Digest);

  //  free(decrypted_hex);
   // free(useful_decrypted_hex);
  // mp_clear(&HEX1);
   //mp_clear(&modulus);
   //mp_clear(&public_key);
   //mp_clear(&decrypted);

   if (strcmp(useful_decrypted_hex,HEX_format_Digest)==0){
       return 1;// valid
   }else{

    return 0;// not valid
   }


}

// this function is for encrypting files that has to be remained inside the RFM
void encrypting_File(char *content, key key, char *fname)
{


int i=0;
char buf[50];
int aux;
mp_int public_key,private_key,modulus,aux_int,aux_int_result;
mp_init_multi(&public_key,&private_key,&modulus,&aux_int,&aux_int_result,NULL);
mp_read_radix(&public_key,"65537",10);//
//mp_int ;
//mp_init(&private_key);
mp_read_radix(&private_key,key.private_exponent,10);//
//mp_int modulus;
//mp_init(&modulus);
mp_read_radix(&modulus,key.modulus,10);//

char encrypted_content[30000];//=(char*) malloc(sizeof(char)*20000);
//mp_int aux_int;
//mp_init(&aux_int);
//mp_int aux_int_result;
//mp_init(&aux_int_result);
char snum[10];
while(content[i]!='\0'){
    aux=int(content[i]);
    sprintf(snum, "%d",aux );
    mp_read_radix(&aux_int, snum,10);
  // printf("\n %c  %d \n",content[i],aux);
    mp_exptmod(&aux_int,&private_key,&modulus,&aux_int_result);

    mp_to_decimal(&aux_int_result,buf,sizeof(buf));
  //  printf("\n modulus product==\n%s\n\n",buf);
  //  break;

    if(i==0){
    strcpy(encrypted_content,buf);
    strcat(encrypted_content,";");
    }else{
        strcat(encrypted_content,buf);
        strcat(encrypted_content,";");
    }



    i=i+1;


}
strcat(encrypted_content,"\0");
mp_clear_multi(&public_key,&private_key,&modulus,&aux_int,&aux_int_result,NULL);
FILE *fptr_encrypted;

fptr_encrypted=fopen(fname,"w");
fprintf(fptr_encrypted,"%s",encrypted_content);
fclose(fptr_encrypted);

//mp_clear(&aux_int_result);
//mp_clear(&modulus);
//mp_clear(&public_key);
//mp_clear(&private_key);

}

// this function is for decrypting files that has to be remained inside the RFM
void decrypting_File(char *file , key key, char *result){
    FILE *fptr;
    fptr=fopen(file,"r");

    char buf[50];

    mp_int public_key;
    mp_init(&public_key);
    mp_read_radix(&public_key,"65537",10);//
    mp_int private_key;
    mp_init(&private_key);
    mp_read_radix(&private_key,key.private_exponent,10);//
    mp_int modulus;
    mp_init(&modulus);
    mp_read_radix(&modulus,key.modulus,10);//
    mp_int aux_int;
    mp_init(&aux_int);
    mp_int aux_int_result;
    mp_init(&aux_int_result);

    signed char tao;
    char decrypt;
    char buff2[30];
    int k=0,di=0;

    while((tao=fgetc(fptr))!=EOF){

        if(tao!=';'){
            buff2[k]=tao;
            k++;
        }else{
            buff2[k]='\0';
            mp_read_radix(&aux_int, buff2,10);
            mp_exptmod(&aux_int,&public_key,&modulus,&aux_int_result);
            mp_to_decimal(&aux_int_result,buf,sizeof(buf));
      // printf("\n modulus product==\n%s\n\n",buf);
            decrypt = atoi(buf);
            result[di]=decrypt;
            di++;
            memset(buff2, 0, sizeof(buff2));
            k=0;
        }

    }
    fclose(fptr);
    result[di]='\0';
    FILE *fptr4;
    fptr4=fopen("/fs/microsd/keys","w");
    fprintf(fptr4,"\n\n%s\n\n",result);
    fclose(fptr4);
mp_clear(&aux_int);
mp_clear(&aux_int_result);
mp_clear(&modulus);
mp_clear(&public_key);
mp_clear(&private_key);

}



// this function is for fetching value of a particular tag in the given string
void fetch_tag(char *content, char *target, char *result, char *fptr){
    char TAG[30];
     int index;
     int i=0;
    strcpy(TAG,"<");
    strcat(TAG,target);
    if(fptr==NULL){

      index=isSubstring(TAG,content,NULL);
      index=index+int(strlen(TAG))+1;
      while(content[index]!='>')
      {
         result[i]=content[index];
         i++;
         index++;
      }
    result[i]='\0';
    }else{
      FILE *fptr_yu;
      fptr_yu=fopen("/fs/microsd/debug_fetch.txt","w");
      FILE *f_aux;
      f_aux=fopen(fptr,"r");
      index=isSubstring(TAG,NULL,f_aux);
      fclose(f_aux);
      f_aux=fopen(fptr,"r");

      index=index+int(strlen(TAG))+1;
      fprintf(fptr_yu,"\n\n%d",index);
      signed char ch;
      int check_count=0;
      while((ch=fgetc(f_aux))!=EOF){
         
         if(check_count>=index){
            if(ch=='>'){
               break;
            }
            result[i]=ch;
            fprintf(fptr_yu,"%c",result[i]);
            i++;
            
         }

         check_count++;
      }
      result[i]='\0';
      fclose(f_aux);
      fprintf(fptr_yu,"\n\n%s",result);
      fclose(fptr_yu);

   }
    
    



}


// this function is for fetching value of a particular tag in the given string
void fetch_tag(char *content, char *target, char *result){
    char TAG[30];
     int index;
    strcpy(TAG,"<");
    strcat(TAG,target);
    index=isSubstring(TAG,content);
    
    index=index+int(strlen(TAG))+1;
    int i=0;

    while(content[index]!='>'){
      result[i]=content[index];
      i++;
      index++;
    }
    result[i]='\0';


}


//function to generate DroneID.txt
void DroneIDcreation( ){

    char DroneID[]="ABBDDJEDNDJK";//any string
    char RFM_version[2]="0";
    char RPAS_category[10]="Small";


 // RFM_public_key: needs to be fetched from PublicPrivateInuse.txt(which is kept encrypted)
 // inside the rfm, decrypted using decrypting key.
    key RFM_private_key;

    char aux_fname[100];
    strcpy(aux_fname,dir);
    char fname[60]="/log/PublicPrivateInuse.txt";//"./log/PublicPrivateInuse.txt";
    strcat(aux_fname,fname);
    memset(fname,0,sizeof(fname));
    strcpy(fname,aux_fname);

    key Inside_RFM;
    strcpy(Inside_RFM.modulus,"180919775566931");
    strcpy(Inside_RFM.private_exponent,"32102716896161");

    char content[5000];//=(char*) malloc(sizeof(char)*5000);

    decrypting_File(fname , Inside_RFM, content);
   // printf("\n\nThe content of decrypted file is :\n\n%s\n",content);

    char target0[30]="Modulus";
    char target1[30]="PrivateKey";

    char value_modulus[2000];//=(char*) malloc(sizeof(char)*2000);
    char value_private[800];//=(char*) malloc(sizeof(char)*800);

    //fetch_tag(content,target0,value_modulus,NULL);
    fetch_tag(content,target0,value_modulus);
    //fetch_tag(content,target1,value_private,NULL);
    fetch_tag(content,target1,value_private);
    strcpy(RFM_private_key.private_exponent,value_private);
    strcpy(RFM_private_key.modulus,value_modulus);




    //printf("vadd %s",value_modulus);
   // free(content);

    char rfm_key0[30]="RFM_public_key_modulus";
    char rfm_key1[30]="RFM_public_key_exponent";



 //    this value is just for writing into the file
    char DigitalSky_public_key[1020]="MIIC8TCCAdmgAwIBAgIJAJRDnqfLydHvMA0GCSqGSIb3DQEBCwUAMA8xDTALBgNVBAMMBHRlc3QwHhcNMTkwMzI2MDcxMTQzWhcNMjkwMzIzMDcxMTQzWjAPMQ0wCwYDVQQDDAR0ZXN0MIIBIjANBgkqhkiG9w0BAQEFAAOCAQ8AMIIBCgKCAQEAq51cjR/mcgd0nWO33O3SM84yu3DRdaG8OMYSqzPixY5R+D8niOTVLZvOtaFROSneP1JmUAcaBn5sFhsFxgpJX8O6ee0m9PqLL+LKjexEs5dZ85IG8GqF+UJABaKfBeTPOgI5NAwoyZPBphzxsra1fH2OV2roaCf4ErMnYluuyey/VfFlHTVgC5+VX2wvO+o6pYUuzdNqCvgYwZrMEDCXm+08iZk/qpLgqgUCQTs8qGu/Y0d/EqwGmv9xN8tyxX+IbaeQM7uztN8PbMf8wY40OqdgNmgaVmMR4mfAO2XJiryR5Y8JACDGf3dhmcDrdtfmNjaHR109o2/wUPhSdWB/3QIDAQABo1AwTjAdBgNVHQ4EFgQUR4p2KJJXG5cZ8STI66RG6l2o7yowHwYDVR0jBBgwFoAUR4p2KJJXG5cZ8STI66RG6l2o7yowDAYDVR0TBAUwAwEB/zANBgkqhkiG9w0BAQsFAAOCAQEAV3uurlHMtyopefBpdGj59eLWCrpRYJLKbDtLFCj+tY1/uiwogUMNsEEHEBeEdwM+PIPuzWZ4tSYQ+SvdCCt4/6e9x+c2/1mZKhnRzL/s9o70RyWZXQO+Dz43B5aIIy/qARUhLxU2NVL42q90pInIh/ltT02IVkcibwDnsM4XJhsSyvQlRyYXdPzDeBjEOVYFpafLbC/7a5FBuNwfNKEMWhOj6AELnC8fWb3maNevhjSH5amGU2XrUp6yIdWUL2HuW7ReSer93Lg6iYujd/aaqk+pWE5bQsC+r2kHpNcpntHJLsd9E1cwzWCJiEM9zK4GXqKV/QDUdPC6FYfEf+ti9A==";

    char Firmware_version[4]="1.0";

    pair_set pairset[7];

    strcpy(pairset[0].tag,"DroneID");
    strcpy(pairset[0].value,DroneID);


    strcpy(pairset[1].tag,"RFM_version");
    strcpy(pairset[1].value,RFM_version);

    strcpy(pairset[2].tag,"RPAS_category");
    strcpy(pairset[2].value,RPAS_category);

    strcpy(pairset[3].tag,"Firmware_version");
    strcpy(pairset[3].value,Firmware_version);

    strcpy(pairset[4].tag,"DigitalSky_public_key");
    strcpy(pairset[4].value,DigitalSky_public_key);

    strcpy(pairset[5].tag,rfm_key0);
    strcpy(pairset[5].value,value_modulus);

    strcpy(pairset[6].tag,rfm_key1);
    strcpy(pairset[6].value,"65537");

    char fileName[100]="/log/DroneID.txt";

      //char aux_fname[100];
      strcpy(aux_fname,dir);
      strcat(aux_fname,fileName);
      memset(fileName,0,sizeof(fileName));
      strcpy(fileName,aux_fname);



  //  free(value_modulus);

    pair_file_write(pairset,7,fileName,RFM_private_key);





}


// this function  creates .txt file when an amendment takes
// place in the hardware(currently only gps)
void HardwareInuseCreation(int gps){
    char DroneID[]="ABBDDJEDNDJK";//any string
    char RFM_version[2]="0";
    char RPAS_category[10]="Small";
    char GPS_ID[20];
    sprintf(GPS_ID, "%d",gps);
 // RFM_public_key: needs to be fetched from PublicPrivateInuse.txt(which is kept encrypted)
 // inside the rfm, decrypted using decrypting key.
    key RFM_private_key;
    char fname[100]="/log/PublicPrivateInuse.txt";

    char aux_fname[100];
    strcpy(aux_fname,dir);
    strcat(aux_fname,fname);
    memset(fname,0,sizeof(fname));
    strcpy(fname,aux_fname);


    key Inside_RFM;
    strcpy(Inside_RFM.modulus,"180919775566931");
    strcpy(Inside_RFM.private_exponent,"32102716896161");

    char content[5000];//=(char*) malloc(sizeof(char)*5000);

    decrypting_File(fname , Inside_RFM, content);
    //printf("\n\nThe content of decrypted file is :\n\n%s\n",content);

    char target0[30]="Modulus";
    char target1[30]="PrivateKey";

    char value_modulus[2000];//=(char*) malloc(sizeof(char)*2000);
    char value_private[800];//=(char*) malloc(sizeof(char)*800);

   // fetch_tag(content,target0,value_modulus,NULL);
    fetch_tag(content,target0,value_modulus);
   // fetch_tag(content,target1,value_private,NULL);
    fetch_tag(content,target1,value_private);
    strcpy(RFM_private_key.private_exponent,value_private);
    strcpy(RFM_private_key.modulus,value_modulus);
   /// RFM private key fetched



   // printf("vadd %s",value_modulus);
  //  free(content);




    pair_set pairset[4];

    strcpy(pairset[0].tag,"DroneID");
    strcpy(pairset[0].value,DroneID);


    strcpy(pairset[1].tag,"RFM_version");
    strcpy(pairset[1].value,RFM_version);

    strcpy(pairset[2].tag,"RPAS_category");
    strcpy(pairset[2].value,RPAS_category);



    strcpy(pairset[3].tag,"GPS_ID");
    strcpy(pairset[3].value,GPS_ID);



    char fileName[100]="/log/HardwareInuse.txt";

    //char aux_fname[100];
    strcpy(aux_fname,dir);
    strcat(aux_fname,fileName);
    memset(fileName,0,sizeof(fileName));
    strcpy(fileName,aux_fname);



   // free(value_modulus);

    pair_file_write(pairset,4,fileName,RFM_private_key);
}



void get_RFM_Key(key *key1){
    char fname[100]="/fs/microsd/log/PublicPrivateInuse.txt";

    /*char aux_fname[100];
    strcpy(aux_fname,dir);
    strcat(aux_fname,fname);
    memset(fname,0,sizeof(fname));
    strcpy(fname,aux_fname);*/

    key Inside_RFM;
    strcpy(Inside_RFM.modulus,"180919775566931");
    strcpy(Inside_RFM.private_exponent,"32102716896161");

    char content[2000];//=(char*) malloc(sizeof(char)*5000);

    decrypting_File(fname , Inside_RFM, content);
   // printf("\n\nThe content of decrypted file is :\n\n%s\n",content);

    char target0[30]="Modulus";
    char target1[30]="PrivateKey";

  /*  char value_modulus[800];//=(char*) malloc(sizeof(char)*2000);
    char value_private[800];//=(char*) malloc(sizeof(char)*800);
*/
  //  fetch_tag(content,target0,value_modulus,NULL);
    fetch_tag(content,target0,key1->modulus);
   // fetch_tag(content,target1,value_private,NULL);
    fetch_tag(content,target1,key1->private_exponent);
   // strcpy(key1->private_exponent,value_private);
    //strcpy(key1->modulus,value_modulus);

}

/// this function combines validation and fetching of tag values of files
// coming from MC/MS  or staying inside RFM
int file_read(char *fname,char *tag,char *result,int file_type){
    // return types
    //if 1 : file is genuine
    //if 0 : file is not genuine
    // if 2: file is not present
    key *RFM_key;
    key RFM_KEY;
    RFM_key=&RFM_KEY;

    FILE *checking_presence;
    checking_presence=fopen(fname,"r");
    if(checking_presence==NULL){
      return 2;
    }
    fclose(checking_presence);

    get_RFM_Key(RFM_key);

    //file_type 0,1
    //0 = this is for file needing RFM key pair to get validated
    //1= this for file needing FMP key pair to get validated
    if (file_type==0){
        // RFM key pair needed to validate
        int validation_stats=Validating_File(fname,RFM_KEY);
        if (validation_stats==1){
            FILE *fptr;
            signed char ch;
            fptr=fopen(fname,"r");
            char *content=(char*) malloc(sizeof(char)*3000);
            int i=0;
            while((ch=fgetc(fptr))!=EOF){
                content[i]=ch;
                i++;
            }
            content[i]='\0';
            fclose(fptr);
            // Now fetch the tag value
            fetch_tag(content,tag,result);
            free(content);
            return 1;
        }
        else{
            return 0;
        }


    }

    else{
        //FMP_key defined here
        key FMP_key;
        strcpy(FMP_key.modulus,"26730675313584941186560749178137844398391258151106035922104564996357015594451047759599219270560135549941477603738118115764757775734508251718265058999430985835716076462248377819064664881403774398289832625702128624297432579628718638031709818597981265137629684514491781630727683626461964227525538709626953674033997874387729484132040534464656945770393501502839625253132645954405219070308408708085422661840555555025926057125687490593846474404214169083806088065122913843606434120163765097852620800693096123077396800135201599717234897499796697020907623367339541602816939734869450634288165688491010942289888232437040231818827");
        strcpy(FMP_key.private_exponent,"9531141503697955439637938672730292016747896920442587354131856429724746943118117766243718142643838557319261617775490854466802870185757493113087536791411008697514583743073684985212980285161716284529911482024006922693782207314400981636708958923930393575173268043008437266686673344571206800262703336536043051327018493439267232775893284980747319937113655022770270510978884647411012877253041625079190653496682971178532563034737589582352362760267741512014615969164270535601279273456871641667043609413678452331175640797803872310658739384463477457294633079275692433997764390282891429002147156439862257036670607929361882536897");
        int validation_stats=Validating_File(fname,FMP_key);
         if (validation_stats==1){
            FILE *fptr;
            signed char ch;
            fptr=fopen(fname,"r");
            char *content=(char*) malloc(sizeof(char)*3000);
            int i=0;
            while((ch=fgetc(fptr))!=EOF){
                content[i]=ch;
                i++;
            }
            content[i]='\0';
            fclose(fptr);
            // Now fetch the tag value
            fetch_tag(content,tag,result);
            free(content);
            return 1;
        }
        else{
            return 0;
        }

    }
}



//amendment in KeyLog.txt with fileID of a file
void KeyLog_Regen(char *fileID){
    FILE *fptr;

    char fileName[100]="/log/KeyLog.txt";

    char aux_fname[100];
    strcpy(aux_fname,dir);
    strcat(aux_fname,fileName);
    memset(fileName,0,sizeof(fileName));
    strcpy(fileName,aux_fname);

    fptr = fopen(fileName,"r");
    signed char ch;
    char *content=(char*) malloc(sizeof(char)*3000);
    int i=0;
    while((ch=fgetc(fptr))!=EOF){
        content[i]=ch;
        i++;
    }
    content[i]='\0';
    printf("\n\n%s\n\n",content);
    char Number_Keyrotation[4];
    char tag_changes[20]="No_Changes";
    char tag_changes1[20]="<No_Changes=";
    fetch_tag(content,tag_changes,Number_Keyrotation);

    int change_aux=atoi(Number_Keyrotation);
    printf("\n%d\n",change_aux);
    change_aux++;

    char changeAux[3];
    sprintf(changeAux,"%hu",(unsigned short)change_aux);////change1
    printf("\n\n%s\n\n",changeAux);
    // number of new lines char
    int index= isSubstring(tag_changes1,content);
    printf("\n\n\n%d\n\n\n",index);
    char* content2=(char*) malloc(3000*sizeof(char));

    int i_count=0;
    while(i_count<index+int(strlen(tag_changes1)))
    {
        content2[i_count]=content[i_count];
        i_count++;
    }
    int i_count_separate=i_count;


    strcat(content2,changeAux);

    i_count_separate=i_count_separate+int(strlen(changeAux));
    while(1){
        if(content[i_count]!='>'){
            printf("\n\nsign :%c\n\n",content[i_count]);
        i_count++;}
        else{
            break;
        }
    }

    content2[i_count_separate]=content[i_count];
    i_count++;
    i_count_separate++;
    content2[i_count_separate]=content[i_count];
    i_count++;
    i_count_separate++;
   // content2[i_count_separate]=content[i_count];
   //     i_count++;
    //    i_count_separate++;
    int new_line_count=0;
    while(1){
        if(content[i_count]=='\n'){
            new_line_count++;
        }
        if(new_line_count>=(change_aux-1)){
            break;
        }
        content2[i_count_separate]=content[i_count];
        i_count++;
        i_count_separate++;
    }
    free(content);
    if((change_aux-1)==0){
    printf("\n\n\n%s\n\n\n",content2);
    strcat(content2,"<");
    strcat(content2,fileID);
    strcat(content2,">\n");
    strcat(content2,"</content>\0");
    printf("\n\n\n%s\n\n\n",content2);}
    else{
        printf("\n\n\n%s\n\n\n",content2);
    strcat(content2,"\n<");
    strcat(content2,fileID);
    strcat(content2,">\n");
    strcat(content2,"</content>\0");
    printf("\n\n\n%s\n\n\n",content2);

    }

    char *signature=(char*) malloc(sizeof(char)*3000);
    key RFM_key;
   // char dd[]="1ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff003031300d060960864801650304020105000420";
  //  strcpy(RFM_key.modulus,"26495127767604337655831691918113833077956334480395285962585476150242592943334533493300410596845191511956515492348526447864429111984398076478393159921716146324008099192785429713391158778980470576914436564260777828713503499078118282954851831718742193757807790192236071369965695862117529898068098551948463482032984966019192036825055149445334134853954472696975028748893285534653718791714099153590644093991103487859923770318942338310035808291238607889065040628342091199251953687656515817829870399170890472577591915379199888970387732102895345522222635522973211224211091186930411354449072698100145010326153546937615803429429");
  //  strcpy(RFM_key.private_exponent,"24021354375322558781058142276736617999389802434596138079632937453577588041977071136990155131350955784632074057774459381406055342415566243621056270140966629267130220147961070276489248399064064585487617556117107847879807603464052857723291974564966716033712609329726458163489658005940146657360121454440603066594695330718634482791509166056105316975123379501455703860598165895993179884232993174716442550760635024519308123321406513256918808853580242186763289628649427009957647963922327454796638500088631677702805825426368000585849292006901367973110968487920147218164610608698409129440488923014780703396073978252084600942649");
    get_RFM_Key(&RFM_key);
    signing_support(RFM_key,content2,signature);
    //content2[int(strlen(content2))-1]='\n';
    strcat(content2,"\n<Sign>\n");
    strcat(content2,"<Signature=");
    strcat(content2,signature);
    strcat(content2,">\n");
    strcat(content2,"</Sign>");
    strcat(content2,"\0");

    fclose(fptr);

   char fileName1[100]="/log/KeyLog.txt";
    //char aux_fname[100];
    strcpy(aux_fname,dir);
    strcat(aux_fname,fileName1);
    memset(fileName1,0,sizeof(fileName));
    strcpy(fileName1,aux_fname);
    remove(fileName1);

    printf("\n\n%s\n\n",content2);
    char aux_copy[3000];
    strcpy(aux_copy,content2);
    free(content2);
    FILE *fptr2;

    char filename2[100]="/log/KeyLog.txt";
    strcpy(aux_fname,dir);
    strcat(aux_fname,filename2);
    memset(filename2,0,sizeof(filename2));
    strcpy(filename2,aux_fname);

    fptr2=fopen(filename2,"w");
    fprintf(fptr2,"%s",aux_copy);

    fclose(fptr2);




}



// this function will tell if the current time lies inside the start and end
bool In_Time(Date_time current, Date_time start, Date_time end)
{

	// full check
   // printf("/n/n Inside In_Time function /n/n");
	if(start.Year<=current.Year && current.Year<=end.Year)
    {
        printf("\n year is ok \n ");
		if(start.Month<=current.Month && current.Month<=end.Month)
        {
            printf("\n Month is ok \n");
            if(start.date<=current.date && current.date<=end.date)
            {
                printf("\n Date is ok \n");
            }
		    else{
                printf("%\nThis is when date is out of limits d \n",current.Year);
			    return 0;
		    }
		}else{
			return 0;
		}
	}
	else{
        	return 0;
	}
	//time check
	if(start.Hours<=current.Hours || current.Hours<=end.Hours)//||
    {
        printf(" \n Hours is ok \n");
        if(start.Hours==current.Hours||current.Hours==end.Hours){
            printf("both the hours are equal");
            if(start.Hours==current.Hours)
            {
                printf("\nstart hours are equal");
                if(start.Minutes<=current.Minutes)
                {
                    printf(" minutes are ok");
                    if(start.Minutes==current.Minutes)
                    {
                        printf("miutes are equal");
                        if(start.Seconds<=current.Seconds){
                            return 1;
                        }else
                        {   return 0;
                        }
                    }else
                    {
                        return 1;}
                }else
                {
                    return 0;
                }
            }
            else{
                if(end.Minutes>=current.Minutes)
                {
                    if(end.Minutes==current.Minutes){
                        if(end.Seconds>current.Seconds){
                            return 1;
                        }else{
                            return 0;
                        }
                    }else{return 1;}
                }else
                {
                    return 0;
                }
            }
        }
        return 1;
    }else
    {
        return 0;
    }
return 0;
}



// this function will convert string format of time to Date_time struct
void conversion_to_dateTime(Date_time *dt,char *str){
    char aux_buf1[60];
    int length_str=0;
    int i_iter=0;
        while(length_str<(int)strlen(str))
        {
                if (str[i_iter]=='.')
                {
                    break;
                }

                if ( str[i_iter]!='T' && str[i_iter]!='-' && str[i_iter]!=':')
                {
                    aux_buf1[length_str]=str[i_iter];

                    if(i_iter==3)
                    {
                        aux_buf1[length_str+1]='\0';
                       // printf(" %s ",aux_buf1);
                        dt->Year=strtol(aux_buf1,NULL,10);
                      //  printf("year: %d ",dt->Year);
                        memset(aux_buf1,'0',4);

                    }
                    if(i_iter==6){
                        aux_buf1[length_str+1]='\0';
                        printf(" %s ",aux_buf1);
                        dt->Month=strtol(aux_buf1,NULL,10);
                       // printf("month: %d ",dt->Month);
                        memset(aux_buf1,'0',7);
                    }
                    if(i_iter==9){
                        aux_buf1[length_str+1]='\0';
                        printf(" %s ",aux_buf1);
                        dt->date=strtol(aux_buf1,NULL,10);
                       // printf("date: %d ",dt->date);
                        memset(aux_buf1,'0',10);
                    }
                    if(i_iter==12){
                        aux_buf1[length_str+1]='\0';
                        printf(" %s ",aux_buf1);
                        dt->Hours=strtol(aux_buf1,NULL,10);
                      //  printf("Hours: %d ",Rv.start.Hours);
                        memset(aux_buf1,'0',10);
                    }
                    if(i_iter==15){
                        aux_buf1[length_str+1]='\0';
                        printf(" %s ",aux_buf1);
                        dt->Minutes=strtol(aux_buf1,NULL,10);
                      //  printf("Minutes: %d ",Rv.start.Minutes);
                        memset(aux_buf1,'0',12);
                    }
                    if(i_iter==18)
                    {
                        aux_buf1[length_str+1]='\0';
                        printf(" %s ",aux_buf1);
                        dt->Seconds=strtol(aux_buf1,NULL,10);
                     //   printf("Seconds: %d ",Rv.start.Seconds);
                        memset(aux_buf1,'0',12);
                    }

                        length_str++;
                }
                i_iter++;

        }
        //printf("\n\n");
        memset(aux_buf1,'0',60);


}

//this function will tell if current time lies inside the start and end time
int check_time(char *start,char *end){
    // getting current time
    int vehi_gps_pos=orb_subscribe(ORB_ID(vehicle_gps_position));
    time_t timestamp;
    struct vehicle_gps_position_s raw;
    orb_copy(ORB_ID(vehicle_gps_position),vehi_gps_pos,&raw);
    timestamp=(raw.time_utc_usec)/1000000;//microsecond to seconds
	printf("\ntimestamp is seconds: %u\n",timestamp);
    struct tm  ts;
    Date_time DT;//current time
    ts = *localtime(&timestamp);
    DT.Year=ts.tm_year+1900;
    DT.Month=ts.tm_mon+1;
    DT.date=ts.tm_mday;
	DT.Hours=ts.tm_hour;
	DT.Minutes=ts.tm_min;
	DT.Seconds=ts.tm_sec;
    // converting start and end to the struct format of Date_time
    Date_time startTime;
    Date_time endTime;

    conversion_to_dateTime(&startTime,start);
    conversion_to_dateTime(&endTime,end);

    int intime_status=In_Time(DT,startTime,endTime);
    if(intime_status==1){
        return 1;// in time
    }

    return 0;// not in time

}


// Check_recentPA() function
/*
check the validity of the recentPA.txt
if the current time lies outside the start time and end time, then start the bundling
if it lies inside then check for number of flights done and number of permissible flights
when landing or free fall instance is noted down, at the same time regenerate recentPA.txt
with the number of updated flights done

return : 0 :: (not in time)start bundling and wait for fetched file from MC.
return : 1 :: (frequency complete)start bundling and wait for fetched file from MC
return : 2 :: no need to start bundling
return : 3 :: invalid recentPA.txt file, someone is trying to hack


retunr
*/
int check_recentPA(char *paID){
    char filename[100]="/log/recentPA.txt";

    char aux_fname[100];
    strcpy(aux_fname,dir);
    strcat(aux_fname,filename);
    memset(filename,0,sizeof(filename));
    strcpy(filename,aux_fname);

    key RFM_key;
    get_RFM_Key(&RFM_key);
    int valid_status= Validating_File(filename , RFM_key);
    if(valid_status!=1){
        //file is not valid
        return 3;
    }
    // file is valid and can be further used
    FILE *fptr;
    fptr=fopen(filename,"r");
    int i=0;
    signed char aux;
    char content[5000];
    while((aux=fgetc(fptr))!=EOF){
        content[i]=aux;
        i++;
    }
    content[i]='\0';
    fclose(fptr);
    //tags to be fetched are:
    //1) start time and end time
    //2) allowable frequency
    //3) frequency of done flights

    //fetching start time and end time
    char start_time_tag[50]="Start_time";
    char start_time_content[70];
    char end_time_tag[40]="End_time";
    char end_time_content[70];
    char PA_ID_tag[50]="permissionArtifactId";

    fetch_tag(content,start_time_tag,start_time_content);
    fetch_tag(content,end_time_tag,end_time_content);
    fetch_tag(content,PA_ID_tag,paID);

    // now, checking if current time lies inside the limits
    int check_status=check_time(start_time_content,end_time_content);
    if(check_status!=1){
         // not in time :: start bundling
        return 0;
    }
    // at this point it is sure that current time is in limits

    // now check for frequency
    // two frequencies need to be fetched :1) allowable 2)done_flights
    char allow_frequency_tag[50]="allowable_frequency";
    char done_frequency_tag[50]="done_frequency";

    char allow_freq_content[60];
    char done_freq_content[60];

    fetch_tag(content,allow_frequency_tag,allow_freq_content);
    fetch_tag(content,done_frequency_tag,done_freq_content);

    if(strcmp(allow_freq_content,done_freq_content)==0){
        return 1;
    }



    //in time and in_limit_frequency :: no need to start bundle

    return 2;
}



////////////////////////////  PA_EXTRACT.cpp

ParsedData parse_artifact()
{


    ParsedData result{};

    char file_name[100]="/log/permission_artifact_breach.xml"; // read mode
    char aux_fname[100];
    strcpy(aux_fname,dir);
    strcat(aux_fname,file_name);
    memset(file_name,0,sizeof(file_name));
    strcpy(file_name,aux_fname);

    char  start_time[200], end_time[200];//,long_lat_coords[20][20];

    char Reference_canonilized[5000];// variable to hold canonicalized reference section

    char xmlExeclusive[5000];
    char outputih[5000];
        // This function will extract the reference section out of the xml file (PA)
	    //and willl perform canonicalization step over it
    Reference_canon(file_name,Reference_canonilized);
    xmlExeclusiveCanon(Reference_canonilized,outputih);
    cleanerXML(outputih,xmlExeclusive);

    char tag_flight_start[100]="flightStartTime=";

    char tag_flight_end[100]="flightEndTime=";
    printf("\ncanon part %s\n",xmlExeclusive);
    int index_start_flight=isSubstring(tag_flight_start,xmlExeclusive);
    printf("\nindex start %d\n",index_start_flight);
    index_start_flight=index_start_flight+int(strlen(tag_flight_start))+1;
    printf("\nstart char %c\n",xmlExeclusive[index_start_flight]);
    int aux_fetch1=0;
    while(xmlExeclusive[index_start_flight]!='\"'){


        start_time[aux_fetch1]=xmlExeclusive[index_start_flight];
        index_start_flight++;
        aux_fetch1++;

    }
    start_time[aux_fetch1]='\0';
    printf("\n\nstart time is :%s\n\n",start_time);


    int index_end_flight=isSubstring(tag_flight_end,xmlExeclusive);
       printf("\nindex start %d\n",index_end_flight);
    index_end_flight=index_end_flight+int(strlen(tag_flight_end))+1;
    aux_fetch1=0;
    while(xmlExeclusive[index_end_flight]!='\"'){


        end_time[aux_fetch1]=xmlExeclusive[index_end_flight];
        aux_fetch1++;
        index_end_flight++;
    }
    end_time[aux_fetch1]='\0';
    printf("\n\nend time is :%s\n\n",end_time);


    strcpy(result.start_time,start_time);
    printf("\nstart time inside parseartifact:%s\n",start_time);
    strcpy(result.end_time , end_time);

    char tag_coordinates[60]="<Coordinates>";

    char tag_coordinates_end[60]="</Coordinates>";

    int index_coordinates_start=isSubstring(tag_coordinates,xmlExeclusive);

    int index_coordinates_end=isSubstring(tag_coordinates_end,xmlExeclusive);

    index_coordinates_start=index_coordinates_start+int(strlen(tag_coordinates));

    int co_ind= index_coordinates_start;
    tag_value geo_coordinate[10];
    int number_of_tags=0;
    int start1=0;
    int start2=0;
    int end1=0;
    int end2=0;
    char aux_tag[100];
    char aux_value[100];
     printf("\nhelllloooo\n");
    while(co_ind<index_coordinates_end){

        if(xmlExeclusive[co_ind]==' ')
        {
            //space arrived
            while(xmlExeclusive[co_ind]!='>')
            {
                start1=co_ind+1;
                while(xmlExeclusive[co_ind]!='='){
                    co_ind++;
                }
                end1=co_ind-1;
                start2=co_ind+2;
                co_ind=co_ind+2;

                while(xmlExeclusive[co_ind]!='\"'){
                    co_ind++;
                }
                end2=co_ind-1;
                co_ind++;


                formPairstrings(aux_tag,aux_value,start1,start2,end1,end2,xmlExeclusive);


                strcpy(geo_coordinate[number_of_tags].tag,aux_tag);
                strcpy(geo_coordinate[number_of_tags].value,aux_value);
                memset(aux_tag, ' ', (end1-start1)*sizeof(char));
                memset(aux_value, ' ', (end2-start2)*sizeof(char));



                //printf("%s %s\n",value_aux,tag_aux);
                // add to tag value array


                number_of_tags++;
            }

        }
        co_ind++;
    }
/*
    for(int i=0;i<number_of_tags;i++){
        printf("\n%s\n",geo_coordinate[i].tag);
        printf("\n%s\n",geo_coordinate[i].value);
          printf("\n\n");
    }

*/
    for(int i = 1; i<=4;i++)
        {
            result.long_lat_coords[i-1][0] = atof(geo_coordinate[2*i-2].value);
            result.long_lat_coords[i-1][1] = atof(geo_coordinate[2*i-1].value);
        }




    return result;
}



bool In_Place(ParsedData geo,int latti, int longi){

    double max_lat=-700;
    double min_lat=700;
    for(int i=0;i<4;i++)
    {
        if(geo.long_lat_coords[i][0]<=min_lat){
            min_lat=geo.long_lat_coords[i][0];
        }
        if(geo.long_lat_coords[i][0]>=max_lat){
            max_lat=geo.long_lat_coords[i][0];
        }
    }
    double max_lon=-700;
    double min_lon=700;
    for(int i=0;i<4;i++)
    {
        if(geo.long_lat_coords[i][0]<=min_lon){
            min_lon=geo.long_lat_coords[i][1];
        }
        if(geo.long_lat_coords[i][0]>=max_lon){
            max_lon=geo.long_lat_coords[i][1];
        }
    }

    if(double(latti)>=min_lat && double(latti)<=max_lat ){
        if(double(longi)>=min_lon && double(longi)<=max_lon ){
            return 1;
        }else{
            return 0;
        }

    }else{
        return 0;
    }

}



bool In_Time(Date_time current, GEO_DATE_TIME_XML Xml)
{

	// full check
   // printf("/n/n Inside In_Time function /n/n");
	if(Xml.start.Year<=current.Year && current.Year<=Xml.end.Year)
    {
        printf("\n year is ok \n ");
		if(Xml.start.Month<=current.Month && current.Month<=Xml.end.Month)
        {
            printf("\n Month is ok \n");
            if(Xml.start.date<=current.date && current.date<=Xml.end.date)
            {
                printf("\n Date is ok \n");
            }
		    else{
                printf("%\nThis is when date is out of limits d \n",current.Year);
			    return 0;
		    }
		}else{
			return 0;
		}
	}
	else{
        	return 0;
	}
	//time check
	if(Xml.start.Hours<=current.Hours || current.Hours<=Xml.end.Hours)//||
    {
        printf(" \n Hours is ok \n");
        if(Xml.start.Hours==current.Hours||current.Hours==Xml.end.Hours){
            printf("both the hours are equal");
            if(Xml.start.Hours==current.Hours)
            {
                printf("\nstart hours are equal");
                if(Xml.start.Minutes<=current.Minutes)
                {
                    printf(" minutes are ok");
                    if(Xml.start.Minutes==current.Minutes)
                    {
                        printf("miutes are equal");
                        if(Xml.start.Seconds<=current.Seconds){
                            return 1;
                        }else
                        {   return 0;
                        }
                    }else
                    {
                        return 1;}
                }else
                {
                    return 0;
                }
            }
            else{
                if(Xml.end.Minutes>=current.Minutes)
                {
                    if(Xml.end.Minutes==current.Minutes){
                        if(Xml.end.Seconds>current.Seconds){
                            return 1;
                        }else{
                            return 0;
                        }
                    }else{return 1;}
                }else
                {
                    return 0;
                }
            }
        }
    }else
    {
        return 0;
    }
return 0;
}


void data_fetch_delete(){
    // this function is called when pa has been validated for time, place, sign and drone ID
    // data that needs to be entered in recentPA.txt
    // 1)time : start time and end time
    // 2)location: gps coordinates (4)
    // 3)frequency upto which drone can make flights
    // 4)pa_id
    // 5)maximum altitude
    // 6) number of flghts done
    // 7)previous log hash

    ParsedData data= parse_artifact();
    FILE *fptr;
    char fileName[100]="/log/permission_artifact_breach.xml";
    char aux_fname[100];
    strcpy(aux_fname,dir);
    strcat(aux_fname,fileName);
    memset(fileName,0,sizeof(fileName));
    strcpy(fileName,aux_fname);

    fptr=fopen(fileName,"r");
    char content[5000];
    signed char ch;
    int i=0;
    while((ch=fgetc(fptr))!=EOF){
        content[i]=ch;
        i++;
    }
    content[i]='\0';
    fclose(fptr);
   // printf("\ndata fetch delete:::%s\n",content);
    // fetching frequency
    char tag_frequency[30]="frequency=";
    char tag_freq[30]="frequency_allowable";
    int index_frequency=isSubstring(tag_frequency,content);

    index_frequency=index_frequency+int(strlen(tag_frequency))+1;

    char tag_content_freq[30];

    i=0;
    while(content[index_frequency]!='\"'){
        tag_content_freq[i]=content[index_frequency];
        index_frequency++;
        i++;
    }
    tag_content_freq[i]='\0';//frequencies fetched.



    // fetching maximum altitude
    char tag_maxalt[30]="maxAltitude=";
    char tag_max[30]="maxAltitude";
    int index_maxalt=isSubstring(tag_maxalt,content);

    index_maxalt=index_maxalt+int(strlen(tag_maxalt))+1;

    char tag_content_maxalt[30];

    i=0;
    while(content[index_maxalt]!='\"'){
        tag_content_maxalt[i]=content[index_maxalt];
        index_maxalt++;
        i++;
    }
    tag_content_maxalt[i]='\0';//frequencies fetched.


    // fetching PA ID
    char tag_paID[70]="permissionArtifactId=";
    char tag_PAid[70]="permissionArtifactId";
    int index_paID=isSubstring(tag_paID,content);

    index_paID=index_paID+int(strlen(tag_paID))+1;

    char tag_content_paID[70];

    i=0;
    while(content[index_paID]!='\"'){
        tag_content_paID[i]=content[index_paID];
        index_paID++;
        i++;
    }
    tag_content_paID[i]='\0';//frequencies fetched.
    printf("\n\ntag_content_paID::::%s\n\n",tag_content_paID);
    pair_set pairset[14];

    strcpy(pairset[0].tag,tag_PAid);
    strcpy(pairset[0].value,tag_content_paID);


    strcpy(pairset[1].tag,tag_freq);
    strcpy(pairset[1].value,tag_content_freq);

    strcpy(pairset[2].tag,tag_max);
    strcpy(pairset[2].value,tag_content_maxalt);


   // int fetch_i=0;
   // char start_time[60];
    //char end_time[60];
    char aux_float[40];
   // sprintf(aux_float,"%f",data.start_time);


    //// start time

    strcpy(pairset[3].tag,"Start_time");
    strcpy(pairset[3].value,data.start_time);

    //// end time
    strcpy(pairset[4].tag,"End_time");
    strcpy(pairset[4].value,data.end_time);

    sprintf(aux_float,"%f",data.long_lat_coords[0][0]);
    //char *f1=gcvt(data.long_lat_coords[0][0], 8, aux_float);
    //printf("\n%s\n",f1);
    printf("\n mainnf %s\n",aux_float);
    //lattitude1
    strcpy(pairset[5].tag,"lattitude1");
    strcpy(pairset[5].value,aux_float);
    memset(aux_float,' ',sizeof(aux_float));

    //longitude1
    sprintf(aux_float,"%f",data.long_lat_coords[0][1]);
   //f1=gcvt(data.long_lat_coords[0][1], 8, aux_float);
    //printf("\n%s\n",f1);
    printf("\n mainnf %s\n",aux_float);
    strcpy(pairset[6].tag,"longitude1");
    strcpy(pairset[6].value,aux_float);
    memset(aux_float,' ',sizeof(aux_float));

    //lattitude2
    sprintf(aux_float,"%f",data.long_lat_coords[1][0]);
   // f1=gcvt(data.long_lat_coords[1][0], 8, aux_float);
   //printf("\n%s\n",f1);
    printf("\n mainnf %s\n",aux_float);
    strcpy(pairset[7].tag,"lattitude2");
    strcpy(pairset[7].value,aux_float);
    memset(aux_float,' ',12*sizeof(char));

    //longitude2
   // printf("printitng %d",int((10^7)*data.long_lat_coord)s[2][1]));
    sprintf(aux_float,"%f",data.long_lat_coords[2][1]);
   // f1=gcvt(data.long_lat_coords[2][1], 8, aux_float);
    printf("\n mainnf %f %s\n",data.long_lat_coords[2][1],aux_float);
   // printf("\n%s\n",f1);
    strcpy(pairset[8].tag,"longitude2");
    strcpy(pairset[8].value,aux_float);
    memset(aux_float,' ',sizeof(aux_float));

    strcpy(pairset[9].tag,"frequencies_done");
    strcpy(pairset[9].value,"0");


    strcpy(pairset[10].tag,"fetch_required");
    strcpy(pairset[10].value,"0");


   strcpy(pairset[11].tag,"previous_log_hash");
   strcpy(pairset[11].value,"None");


   //publishing Permission Artefact information to pa_data.msg

	struct pa_data_s pa_data_raw;
	memset(&pa_data_raw, 0, sizeof(pa_data_raw));
	orb_advert_t pa_data_raw_pub = orb_advertise(ORB_ID(pa_data), &pa_data_raw);

   pa_data_raw.lattitude1=data.long_lat_coords[0][0];
   pa_data_raw.longitude1=data.long_lat_coords[0][1];
   pa_data_raw.lattitude2=data.long_lat_coords[1][0];
   pa_data_raw.longitude2=data.long_lat_coords[1][1];
   pa_data_raw.allowable_height=strtol(tag_content_maxalt,NULL,10);

   Date_time startTime;
   Date_time endTime;
   conversion_to_dateTime(&startTime,data.start_time);
   conversion_to_dateTime(&endTime,data.end_time);

   pa_data_raw.start_hours=startTime.Hours;
   pa_data_raw.start_date=startTime.date;
   pa_data_raw.start_minutes=startTime.Minutes;
   pa_data_raw.start_month=startTime.Month;
   pa_data_raw.start_year=startTime.Year;
   pa_data_raw.start_seconds=startTime.Seconds;


   pa_data_raw.end_hours=endTime.Hours;
   pa_data_raw.end_date=endTime.date;
   pa_data_raw.end_minutes=endTime.Minutes;
   pa_data_raw.end_month=endTime.Month;
   pa_data_raw.end_year=endTime.Year;
   pa_data_raw.end_seconds=endTime.Seconds;

   int vehi_gps_pos=orb_subscribe(ORB_ID(vehicle_gps_position));

   struct vehicle_gps_position_s raw;
   orb_copy(ORB_ID(vehicle_gps_position),vehi_gps_pos,&raw);

   pa_data_raw.home_altitude=pow(10,-3)*(raw.alt);// in meters

   pa_data_raw.updated=1;


   orb_publish(ORB_ID(pa_data), pa_data_raw_pub, &pa_data_raw);




   char fileName1[100]="/fs/microsd/log/recentPA.txt";

    //char aux_fname[100];
    /*strcpy(aux_fname,dir);
    strcat(aux_fname,fileName1);
    memset(fileName1,0,sizeof(fileName1));
    strcpy(fileName1,aux_fname);*/





   key RFM_private_key;

   get_RFM_Key(&RFM_private_key);

  //  free(value_modulus);

   pair_file_write(pairset,12,fileName1,RFM_private_key);




   // remove("./log/permission_artifact_breach.xml");

}

int date_time_extract_and_check()
{

	FILE * file;
	int vehi_gps_pos=orb_subscribe(ORB_ID(vehicle_gps_position));
    time_t timestamp;
    struct vehicle_gps_position_s raw;
    orb_copy(ORB_ID(vehicle_gps_position),vehi_gps_pos,&raw);
    int lattitude=raw.lat;
    int longitude=raw.lon;
    int alttitude=raw.alt_ellipsoid;
    printf("\n\ncurrent lattitude: %d \ncurrent longitude: %d \ncurrent height: %d ",lattitude,longitude,alttitude);
    timestamp=(raw.time_utc_usec)/1000000;//microsecond to seconds
	printf("\ntimestamp is seconds: %u\n",timestamp);
    struct tm  ts;
    Date_time DT;//current time
    ts = *localtime(&timestamp);
    DT.Year=ts.tm_year+1900;
    DT.Month=ts.tm_mon+1;
    DT.date=ts.tm_mday;
	DT.Hours=ts.tm_hour;
	DT.Minutes=ts.tm_min;
	DT.Seconds=ts.tm_sec;
	printf("\n%d/%d/%d current date and time from GPS %d:%d:%d \n",DT.Year,DT.Month,DT.date,DT.Hours,DT.Minutes,DT.Seconds);

    char fileName[100]="/log/permission_artifact_breach.xml";
    char aux_fname[100];
    strcpy(aux_fname,dir);
    strcat(aux_fname,fileName);
    memset(fileName,0,sizeof(fileName));
    strcpy(fileName,aux_fname);

	file = fopen(fileName, "r");// ./log/ for posix /log/ for nuttx
	if (file)
    {
        fclose(file);
        GEO_DATE_TIME_XML Rv;
        ParsedData data= parse_artifact();
            printf("\n  %s    from parse_artifact function      %s \n ",data.start_time,data.end_time);
            char aux_buf1[60];
            // extracting start year
        printf("\n");
        int length_str=0;
        int i_iter=0;
        while(length_str<(int)strlen(data.start_time))
        {
                if (data.start_time[i_iter]=='.')
                {
                    break;
                }

                if ( data.start_time[i_iter]!='T' && data.start_time[i_iter]!='-' && data.start_time[i_iter]!=':')
                {
                    aux_buf1[length_str]=data.start_time[i_iter];

                    if(i_iter==3)
                    {
                        aux_buf1[length_str+1]='\0';
                        printf(" %s ",aux_buf1);
                        Rv.start.Year=strtol(aux_buf1,NULL,10);
                        printf("year: %d ",Rv.start.Year);
                        memset(aux_buf1,'0',4);

                    }
                    if(i_iter==6){
                        aux_buf1[length_str+1]='\0';
                        printf(" %s ",aux_buf1);
                        Rv.start.Month=strtol(aux_buf1,NULL,10);
                        printf("month: %d ",Rv.start.Month);
                        memset(aux_buf1,'0',7);
                    }
                    if(i_iter==9){
                        aux_buf1[length_str+1]='\0';
                        printf(" %s ",aux_buf1);
                        Rv.start.date=strtol(aux_buf1,NULL,10);
                        printf("date: %d ",Rv.start.date);
                        memset(aux_buf1,'0',10);
                    }
                    if(i_iter==12){
                        aux_buf1[length_str+1]='\0';
                        printf(" %s ",aux_buf1);
                        Rv.start.Hours=strtol(aux_buf1,NULL,10);
                        printf("Hours: %d ",Rv.start.Hours);
                        memset(aux_buf1,'0',10);
                    }
                    if(i_iter==15){
                        aux_buf1[length_str+1]='\0';
                        printf(" %s ",aux_buf1);
                        Rv.start.Minutes=strtol(aux_buf1,NULL,10);
                        printf("Minutes: %d ",Rv.start.Minutes);
                        memset(aux_buf1,'0',12);
                    }
                    if(i_iter==18)
                    {
                        aux_buf1[length_str+1]='\0';
                        printf(" %s ",aux_buf1);
                        Rv.start.Seconds=strtol(aux_buf1,NULL,10);
                        printf("Seconds: %d ",Rv.start.Seconds);
                        memset(aux_buf1,'0',12);
                    }

                        length_str++;
                }
                i_iter++;

        }
        printf("\n\n");
        memset(aux_buf1,'0',60);
        // extracting end year
        printf("\n");
        length_str=0;
        i_iter=0;
        while(length_str<(int)strlen(data.end_time)){
            if (data.end_time[i_iter]=='.'){
                break;
            }

            if ( data.end_time[i_iter]!='T' && data.end_time[i_iter]!='-' && data.end_time[i_iter]!=':'){


            aux_buf1[length_str]=data.end_time[i_iter];

            if(i_iter==3){
                aux_buf1[length_str+1]='\0';
                printf(" %s ",aux_buf1);
                Rv.end.Year=strtol(aux_buf1,NULL,10);
                printf("year: %d ",Rv.end.Year);
                memset(aux_buf1,'0',4);

            }
            if(i_iter==6){
                aux_buf1[length_str+1]='\0';
                printf(" %s ",aux_buf1);
                Rv.end.Month=strtol(aux_buf1,NULL,10);
                printf("month: %d ",Rv.end.Month);
                memset(aux_buf1,'0',7);
            }
            if(i_iter==9){
                aux_buf1[length_str+1]='\0';
                printf(" %s ",aux_buf1);
                Rv.end.date=strtol(aux_buf1,NULL,10);
                printf("date: %d ",Rv.start.date);
                memset(aux_buf1,'0',10);
            }
            if(i_iter==12){
                aux_buf1[length_str+1]='\0';
                printf(" %s ",aux_buf1);
                Rv.end.Hours=strtol(aux_buf1,NULL,10);
                printf("Hours: %d ",Rv.end.Hours);
                memset(aux_buf1,'0',10);
            }
            if(i_iter==15){
                aux_buf1[length_str+1]='\0';
                printf(" %s ",aux_buf1);
                Rv.end.Minutes=strtol(aux_buf1,NULL,10);
                printf("Minutes: %d ",Rv.end.Minutes);
                memset(aux_buf1,'0',12);
            }
            if(i_iter==18){
                aux_buf1[length_str+1]='\0';
                printf(" %s ",aux_buf1);
                Rv.end.Seconds=strtol(aux_buf1,NULL,10);
                printf("Seconds: %d ",Rv.end.Seconds);
                memset(aux_buf1,'0',12);
            }

                length_str++;
            }
            i_iter++;

        }


        printf("\n\n");
        printf("\nfrom xml end date and time:: %d/%d/%d and  %d:%d:%d",Rv.end.Year,Rv.end.Month,Rv.end.date,Rv.end.Hours,Rv.end.Minutes,Rv.end.Seconds);
        printf("\n\n");
        printf("\nfrom xml start date and time:: %d/%d/%d and  %d:%d:%d",Rv.start.Year,Rv.start.Month,Rv.start.date,Rv.start.Hours,Rv.start.Minutes,Rv.start.Seconds);


            //int has_artifact = 1;
            //param_set(param_find("PERM_ARTIFACT"),&has_artifact);
            //printf("Permission Artifact Found \n");

        bool in_place=In_Place(data,lattitude,longitude);

        bool in_time=In_Time(DT,Rv);
        if (in_time!=1){
            printf("\n permission denied\n");
            return 0;
        }else{
            if(in_place!=1){
              return 2;
            }
            return 1;

        }

	}
	else
    {
   	//file doesn't exists or cannot be opened (es. you don't have access permission)
	//int has_artifact = 0;
	//param_set(param_find("PERM_ARTIFACT"),&has_artifact);
	return 0;
	}
}

int DroneIDverification(char *paID){
    char file_name[100]="/log/permission_artifact_breach.xml";

    char aux_fname[100];
    strcpy(aux_fname,dir);
    strcat(aux_fname,file_name);
    memset(file_name,0,sizeof(file_name));
    strcpy(file_name,aux_fname);

    char Reference_canonilized[5000];// variable to hold canonicalized reference section

    char content[5000];
    char outputih[5000];
        // This function will extract the reference section out of the xml file (PA)
	    //and willl perform canonicalization step over it
    Reference_canon(file_name,Reference_canonilized);
    xmlExeclusiveCanon(Reference_canonilized,outputih);
    cleanerXML(outputih,content);

    char paID_tag[70]="permissionArtifactId";
    char paID_contains[100];
    fetch_tag(content,paID_tag,paID_contains);
    if(strcmp(paID_contains,paID)!=0){
        return 2;
    }
    char tag[30]="UUID=";
    char tag_content[50];
    int index=isSubstring(tag,content);
    int i=index;
    while(content[i]!='='){
        i++;
    }
    i=i+2;
    int k=0;
    while(content[i]!='\"'){
        tag_content[k]=content[i];
        i++;
        k++;
    }
    tag_content[k]='\0';
    printf("\n%s\n",content);
    printf("\n%s\n",tag_content);
    char dro[40]="ABBDDJEDNDJK";
    if(strcmp(tag_content,dro)==0){
        return 1;
    }else{
        return 0;
    }


}


int DroneIDverification(){
    char file_name[100]="/log/permission_artifact_breach.xml";

    char aux_fname[100];
    strcpy(aux_fname,dir);
    strcat(aux_fname,file_name);
    memset(file_name,0,sizeof(file_name));
    strcpy(file_name,aux_fname);
    char Reference_canonilized[5000];// variable to hold canonicalized reference section

    char content[5000];
    char outputih[5000];
        // This function will extract the reference section out of the xml file (PA)
	    //and willl perform canonicalization step over it
    Reference_canon(file_name,Reference_canonilized);
    xmlExeclusiveCanon(Reference_canonilized,outputih);
    cleanerXML(outputih,content);

    /*char paID_tag[70]="permissionArtifactId";
    char paID_contains[100];
    fetch_tag(content,paID_tag,paID_contains);
    if(strcmp(paID_contains,paID)!=0){
        return 2;
    }*/
    char tag[30]="UUID=";
    char tag_content[50];
    int index=isSubstring(tag,content);
    int i=index;
    while(content[i]!='='){
        i++;
    }
    i=i+2;
    int k=0;
    while(content[i]!='\"'){
        tag_content[k]=content[i];
        i++;
        k++;
    }
    tag_content[k]='\0';
    printf("\n%s\n",content);
    printf("\n%s\n",tag_content);
    char dro[40]="ABBDDJEDNDJK";
    if(strcmp(tag_content,dro)==0){
        return 1;
    }else{
        return 0;
    }


}

void formPairstrings(char *tag,char *value,int st1,int st2,int end1,int end2,char *sample){
    int count=0;
    for(int i=st1;i<=end1;i++){
        tag[count]=sample[i];
        count++;
    }
    tag[count]='\0';
    count =0;
    for(int i=st2;i<=end2;i++){
        value[count]=sample[i];
        count++;
    }
    value[count]='\0';


}

int min(int a,int b){
    if(a<b){
        return a;
    }
    else{
        return b;
    }

}

void sortAndmake(tag_value *aux,int number_of_tags,char *result){
   // tag_value res[10];
   //for(int i=0;i<number_of_tags;i++){
        //printf("\n%s\n",aux[i].tag);
       // printf("\n%s\n",aux[i].value);
     //   printf("\n\n");
   // }

    char aux_tag[100];
    char aux_value[100];

   // int mini=10000;
    int k=0;

    for(int i=0;i<number_of_tags;i++)
    {
       for(int j=i+1;j<number_of_tags;j++)
       {

           k=min(int(strlen(aux[i].tag)),int(strlen(aux[j].tag)));

           for(int v=0;v<k;v++)
           {

               if(int(aux[i].tag[v])> int(aux[j].tag[v])){
                   //swapping needed
                   strcpy(aux_tag,aux[i].tag);
                   strcpy(aux_value,aux[i].value);

                   strcpy(aux[i].tag,aux[j].tag);
                   strcpy(aux[i].value,aux[j].value);

                   strcpy(aux[j].tag,aux_tag);
                   strcpy(aux[j].value,aux_value);



               }else if(int(aux[i].tag[v])< int(aux[j].tag[v])){
                   break;
               }else {
                   continue;
               }
           }


       }

    }
    for(int i=0;i<number_of_tags;i++){
        printf("\n%s\n",aux[i].tag);
        printf("\n%s\n",aux[i].value);
        printf("\n\n");
       if(i==0){
        strcpy(result,aux[i].tag);
       }else{
           strcat(result," ");
           strcat(result,aux[i].tag);
       }
       strcat(result,"=");
        strcat(result,"\"");
        strcat(result,aux[i].value);

        strcat(result,"\"");
    }
    strcat(result,">\0");
    printf("\n%s\n",result);


}


void xmlExeclusiveCanon(char *sample,char *result){
    int i=0;
    int r=0;
    int start1=0,start2=0;
    int end1=0,end2=0;
    int element_head_set=0;
    int check_element=0;

    char *aux_buff=(char*) malloc(1000*sizeof(char));
     char *tag_aux=(char*) malloc(100*sizeof(char));
    char *value_aux=(char*) malloc(100*sizeof(char));


    while(sample[i]!='\0')
    {
        if(sample[i]=='<'){
            check_element=1;
           // element_head_set=1;
           result[r]=sample[i];
            r++;
            i++;
            continue;


        }

        if(check_element==1 && sample[i]=='>' && element_head_set==0){
            //type <xxxx> </xxx>
            //just passby
            check_element=0;
            result[r]=sample[i];
            r++;
            i++;


            continue;

        }
        else if (check_element==1 && sample[i]==' ')
        {
            //type <xxx > this has tag and values
           // element_head_set==1;
            int number_of_tags=0;
            result[r]=sample[i];
            r++;
            tag_value *aux=(tag_value*) malloc(10*sizeof(tag_value));

            //create an element
            while(sample[i]!='>')
            {

                start1=i+1;
                while(sample[i]!='=')
                {
                    i++;
               //     printf("%c",sample[i]);
                }
              //  printf("%c\n",sample[i]);
                end1=i-1;
                start2=i+2;

                i=i+2;
               // printf("\n%c",sample[i]);
                while(sample[i]!='\"')
                {
                    i++;
                }
                end2=i-1;
                i++;//coming at space or >

                // form tag and value strings
             /*   printf(" %c ",sample[start2]);

                printf(" %c \n",sample[end2]);

                printf(" %c ",sample[start1]);
                printf(" %c \n",sample[end1]);*/
                formPairstrings(tag_aux,value_aux,start1,start2,end1,end2,sample);

                strcpy(aux[number_of_tags].tag,tag_aux);
                strcpy(aux[number_of_tags].value,value_aux);
                memset(tag_aux, ' ', (end1-start1)*sizeof(char));
                if(sample[start2]!='\"'){
                    //no value for tag;
                    memset(value_aux, ' ', (end2-start2)*sizeof(char));
                }


                //printf("%s %s\n",value_aux,tag_aux);
                // add to tag value array


                number_of_tags++;

            }
             // send array to the sort and make function
            sortAndmake(aux,number_of_tags,aux_buff);
            // function sends a sorted string which is catenated at the end
            printf("\n on return \n %s \n",aux_buff);
           // strcat(result,aux_buff);
           // strcat(result,"\0");
           // r=r+int(strlen(aux_buff));
           int y=0;
           while(aux_buff[y]!='\0'){
               result[r]=aux_buff[y];
               r++;
               y++;
           }

             memset(aux_buff, ' ', (int(strlen(aux_buff)))*sizeof(char));

                // of the result array
           // element_head_set=0;
            check_element=0;
            free(aux);
            i++;
            continue;
        }

        result[r]=sample[i];
        r++;
        i++;

    }
   free(aux_buff);
   free(value_aux);
   free(tag_aux);
    result[r]='\0';


}

void Sha256_implement(char *Input, char *output, FILE *fptr, int index_count){
   int aux0;
   int count_digest=0;
   int it=0;
   if (fptr==NULL){
      char ch;
      //  unsigned char st_arr[3000];
	   unsigned char *st=(unsigned char*) malloc(2000*sizeof(unsigned char));
      unsigned int iff=0;

        // assigning char to unsigned char (this is done to have asci value in 8 bits rather than 7 bits)
      while(1)
      {
	      ch=Input[iff];
	      if(ch=='\0')
         {
		      break;
	      }
         *(st+iff)=ch;
	      //	printf("%c",*(st+iff));
		   iff++;
      }


      //printf("\n");
      //printf("\n");

      SHA256 sha;

      sha.update(st,iff);
       uint8_t * digest = sha.digest();
      free(st);
        while(it<32)
          {

            aux0=(*(digest+it)<<24)|(*(digest+it+1)<<16)|(*(digest+it+2)<<8)|*(digest+it+3);

            for(int h=28;h>=0;h=h-4){
               output[count_digest] =HEX_DIGITS[(aux0>>h)&0xf];
               count_digest++;
            }

		    it=it+4;
	   }

      output[count_digest]='\0';
   }else
   {
      SHA256 sha;
      sha.update(fptr,index_count);
       uint8_t * digest  = sha.digest();
         while(it<32)
          {

            aux0=(*(digest+it)<<24)|(*(digest+it+1)<<16)|(*(digest+it+2)<<8)|*(digest+it+3);

            for(int h=28;h>=0;h=h-4){
               output[count_digest] =HEX_DIGITS[(aux0>>h)&0xf];
               count_digest++;
            }

		    it=it+4;
	   }

      output[count_digest]='\0';
   }


   /*int it=0;//just an iterator
	//printing digest
	printf("\n");

   while(it<32){
		printf("%02x",*(digest+it));
		it++;
	}*/


//   printf("\nSha of reference section that has been canonicalized :%s\n",Sha_of_Reference);
}


int Is_PA_VAlid(){
   return 1;
	char Tag_Digest_value0[12]="DigestValue";
 //  char Digest_Value0[100];// in base 64
   char *Digest_Value0=(char*) malloc(100*sizeof(char));
  // char Digest_Value_in_hex[100];// in hex
 char *Digest_Value_in_hex=(char*) malloc(100*sizeof(char));

   char Tag_Signed_Value0[17]="SignatureValue";

 //  char Signature_Value0[550];// in base 64
   char *Signature_Value0=(char*) malloc(550*sizeof(char));

  // char Signature_Value_in_hex[500];// in hex
   char *Signature_Value_in_hex=(char*) malloc(500*sizeof(char));

   //char  Extracted_sha_from_Signature_value[100];
   char *Extracted_sha_from_Signature_value=(char*) malloc(100*sizeof(char));

   char file_name[100]="/fs/microsd/log/permission_artifact_breach.xml";//"permission_artifact_1.xml";

   /*char aux_fname[100];
   strcpy(aux_fname,dir);
   strcat(aux_fname,file_name);
   memset(file_name,0,sizeof(file_name));
   strcpy(file_name,aux_fname);*/

   //char Sha_of_Reference[300];
   char *Sha_of_Reference=(char*) malloc(300*sizeof(char));

  // char Sha_of_SignedInfo[300];
   char *Sha_of_SignedInfo=(char*) malloc(300*sizeof(char));



  // char Reference_canonilized[5000];// variable to hold canonicalized reference section
   char *Reference_canonilized=(char*) malloc(5000*sizeof(char));

 //  char xmlExeclusive[5000];
   char *xmlExeclusive=(char*) malloc(5000*sizeof(char));

  // char output[5000];
    char *output=(char*) malloc(5000*sizeof(char));

        // This function will extract the reference section out of the xml file (PA)
	    //and willl perform canonicalization step over it
        Reference_canon(file_name,Reference_canonilized);
        cleanerXML(Reference_canonilized,output);
     //   printf("\noutput:::::::::::::: %s\n",output);
        xmlExeclusiveCanon(output,xmlExeclusive);
        free(Reference_canonilized);
        free(output);

        printf("\n\n execlusive xml  :%s\n\n",xmlExeclusive);

      //  char SignedInfo_canonilized[4000];
      char *SignedInfo_canonilized=(char*) malloc(1000*sizeof(char));
      //char outputS[3000];
    char *outputS=(char*) malloc(3000*sizeof(char));
	    //This section will extract the signed info section out of the xml file (PA)
	    //and will perform canonicalization step over it
	    SignedInfo_canon(file_name,SignedInfo_canonilized);
        cleanerXML(SignedInfo_canonilized,outputS);
        printf("\nSigned Info :%s\n",SignedInfo_canonilized);
        free(SignedInfo_canonilized);

        printf("\n");


	    // SHA value calculated for extracted canonicalized SIgned info and reference section
	    //start
        char ch;
      //  unsigned char st_arr[3000];
	    unsigned char *st=(unsigned char*) malloc(3000*sizeof(unsigned char));
        unsigned int iff=0;
        // assigning char to unsigned char (this is done to have asci value in 8 bits rather than 7 bits)
      while(1){
	      ch=xmlExeclusive[iff];
	      if(ch=='\0'){
		      break;
	      }
         *(st+iff)=ch;
	      //	printf("%c",*(st+iff));
		   iff++;
      }
      free(xmlExeclusive);

	//printf("\n");
	//printf("\n");

        SHA256 sha;

	sha.update(st,iff);

	uint8_t * digest = sha.digest();
   free(st);

   int it=0;//just an iterator
	//printing digest
	printf("\n");

        while(it<32){
		printf("%02x",*(digest+it));
		it++;
	}
       int aux0;
       int count_digest=0;
       it=0;
       while(it<32)
          {

            aux0=(*(digest+it)<<24)|(*(digest+it+1)<<16)|(*(digest+it+2)<<8)|*(digest+it+3);

            for(int h=28;h>=0;h=h-4){
               Sha_of_Reference[count_digest] =HEX_DIGITS[(aux0>>h)&0xf];
               count_digest++;
            }

		    it=it+4;
	   }

      Sha_of_Reference[count_digest]='\0';
   printf("\nSha of reference section that has been canonicalized :%s\n",Sha_of_Reference);

   //   printf("\n");


   char ch1;
  // unsigned char st1_arr[3000];
   unsigned char *st1=(unsigned char*) malloc(3000*sizeof(unsigned char));
   unsigned int i1=0;
    //// assigning char to unsigned char
    while(1){
	   ch1=outputS[i1];
	   if(ch1=='\0'){
		   break;
	   }
        *(st1+i1)=ch1;
	//	printf("%c",*(st1+i1));
		i1++;

    }
   free(outputS);

	printf("\n");
	printf("\n");

   SHA256 sha1;

	sha1.update(st1,i1);

	uint8_t * digest1 = sha1.digest();
   free(st1);
   it=0;//just an iterator
	//printing digest
	printf("\n");

   while(it<32){
		printf("%02x",*(digest1+it));
		it++;
	}
   int aux1;
   count_digest=0;
   it=0;
   while(it<32)
      {

            aux1=(*(digest1+it)<<24)|(*(digest1+it+1)<<16)|(*(digest1+it+2)<<8)|*(digest1+it+3);

            for(int h=28;h>=0;h=h-4){
               Sha_of_SignedInfo[count_digest] =HEX_DIGITS[(aux1>>h)&0xf];
               count_digest++;
            }

		    it=it+4;
	   }

   Sha_of_SignedInfo[count_digest]='\0';
   printf("\nsha of SignedInfo section that has been canonicalized :%s\n",Sha_of_SignedInfo);

//	printf("\n");


   // SHA value calculated for extracted canonicalized SIgned info and reference section
   // end

   //In coming functions :
   //1) extracting digest value
   //2)converting it to hex from base64
   //same above two steps for Signature value.

    // start
    // storing Digest value
   getTagvalue(Tag_Digest_value0,Digest_Value0,file_name);
   //  printf("\nhello  %d\n",Digestv.len);


  // printf("\n Digest value in the permission artefact is %s\n",Digest_Value0);

  // printf("\n");
   base64decoder(Digest_Value0, Digest_Value_in_hex);

   free(Digest_Value0);
  // printf("\n");
  // printf("\n Digest value in hex format %s \n",Digest_Value_in_hex);


   // storing Signed value
   getTagvalue(Tag_Signed_Value0,Signature_Value0, file_name);



   printf("\n Signature value in the permission artefact is%d\n %s\n",int(strlen(Signature_Value0)),Signature_Value0);
   printf("\n");

   base64decoder(Signature_Value0, Signature_Value_in_hex);
   free(Signature_Value0);
   printf("\n");
   printf("\n Signature value in hex format %s \n",Signature_Value_in_hex);

   //end

   // message: bignum object created for storing value of Signature value
   // now Performing modulus operation on bignumber integers

   mp_int message,modulus,public_key,Decrypted;
   mp_init_multi(&message,&modulus,&public_key,&Decrypted,NULL);
   mp_read_radix(&message,Signature_Value_in_hex,16);
   free(Signature_Value_in_hex);

   /// Public key of dgca: This has to be taken from a reserved file in directory./ firmware update

 //  mp_int ;
    char Modulus[513]="d6912a773335d3a193c742071762794e26bcac49d78b6b65a784e1c18d18f2f88b9f13d7fc41a5b9e11e3d75905b0f8373dbac5658962940d104711467814b7319774014644df82e9764b4e852cb7547d56de3224aada000db231b1356ffe86391b3eca2ddf67d44f40ea6e84092cc67387db3a2042487b8dacfe5b588738973a59f5fa813a33e4bb8fbafa407794db990a307201b0a0fbd92ddd181a868a6cfa64625353714e9b1de6f81f3addbf5b202030dbd9e51385db9314591f01af01f07cc92f18ebe60545cea53eb54438e251fd88e7380b5a0612c29ccb7dfd88f7a7d7f0dd07a9b596923602798ad9f1bbb89fffd3cdb8fb94c48640d27fa681e47";
  // char Modulus[513]="ab9d5c8d1fe67207749d63b7dcedd233ce32bb70d175a1bc38c612ab33e2c58e51f83f2788e4d52d9bceb5a1513929de3f526650071a067e6c161b05c60a495fc3ba79ed26f4fa8b2fe2ca8dec44b39759f39206f06a85f9424005a29f05e4cf3a0239340c28c993c1a61cf1b2b6b57c7d8e576ae86827f812b327625baec9ecbf55f1651d35600b9f955f6c2f3bea3aa5852ecdd36a0af818c19acc1030979bed3c89993faa92e0aa0502413b3ca86bbf63477f12ac069aff7137cb72c57f886da79033bbb3b4df0f6cc7fcc18e343aa76036681a566311e267c03b65c98abc91e58f090020c67f776199c0eb76d7e6363687475d3da36ff050f85275607fdd";
   mp_read_radix(&modulus,Modulus,16);

  // mp_int ;
   mp_read_radix(&public_key,"65537",10);

   //mp_int ;
   mp_exptmod(&message, &public_key,&modulus,&Decrypted);   //                       this part for decrypting encrypted text

   char *message_string_hex=(char*) malloc(513*sizeof(char));
   mp_to_hex(&Decrypted,message_string_hex,sizeof(message_string_hex));
   printf("\n Decrypted signature value : %s\n",message_string_hex);//this is the decrypted message
   // At this point message has been decrypted
   // Now having decrypted message in small letters
   mp_clear_multi(&message,&modulus,&public_key,&Decrypted,NULL);
   char padding_SHA256[]="1ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff003031300d060960864801650304020105000420";
   for(int o=0;o<444;o++){
       if(padding_SHA256[o]!=small_letter(message_string_hex[o])){
           return 2;

       }
   }
  // printf("\n\n");
   int k_cout=0;
   for(int i2=445;message_string_hex[i2]!='\0';i2++){
      Extracted_sha_from_Signature_value[k_cout]=small_letter(message_string_hex[i2]);
      k_cout++;
     // printf("%c",message_string_hex[i]);
   }
   free(message_string_hex);
   //with Extracted SHA we are reffering to the expected SHA of SIgned info section
   Extracted_sha_from_Signature_value[k_cout]='\0';//
   printf("\n%s\n",Extracted_sha_from_Signature_value);

   if(strcmp(Extracted_sha_from_Signature_value,Sha_of_SignedInfo)==0){
        printf("SEE if next two string strings are same");
       // printf("\n%s\n ",Sha_of_Reference);
       // printf("%s\n",Digest_Value_in_hex);
        if(strcmp(Sha_of_Reference,Digest_Value_in_hex)==0){
            printf("\n 100000000000PA is valid and non tampered \n");
            int pa_done=1;
            param_set(param_find("PA_CHECK"),&pa_done);
            free(Sha_of_SignedInfo);
            free(Sha_of_Reference);
            free(Extracted_sha_from_Signature_value);
            free(Digest_Value_in_hex);
            return 1;
        }else{
            printf("\n PA is not valid\n");
            free(Sha_of_SignedInfo);
            free(Sha_of_Reference);
            free(Digest_Value_in_hex);
            free(Extracted_sha_from_Signature_value);
            return 0;
        }

   }else{
      printf("\n PA is not valid \n");
    //  PX4_INFO("PA is valid???????????");
      free(Sha_of_SignedInfo);
      free(Sha_of_Reference);
      free(Digest_Value_in_hex);
      free(Extracted_sha_from_Signature_value);
      return 0;
   }
/// at this point validation is complete
}


int Is_PA_VAlid1(){

	char Tag_Digest_value0[12]="DigestValue";
 //  char Digest_Value0[100];// in base 64
   char *Digest_Value0=(char*) malloc(100*sizeof(char));
  // char Digest_Value_in_hex[100];// in hex
 char *Digest_Value_in_hex=(char*) malloc(100*sizeof(char));

   char Tag_Signed_Value0[17]="SignatureValue";

 //  char Signature_Value0[550];// in base 64
   char *Signature_Value0=(char*) malloc(550*sizeof(char));

  // char Signature_Value_in_hex[500];// in hex
   char *Signature_Value_in_hex=(char*) malloc(500*sizeof(char));

   //char  Extracted_sha_from_Signature_value[100];
   char *Extracted_sha_from_Signature_value=(char*) malloc(100*sizeof(char));

   char file_name[100]="/fs/microsd/log/permission_artifact_breach.xml";//"permission_artifact_1.xml";

   /*char aux_fname[100];
   strcpy(aux_fname,dir);
   strcat(aux_fname,file_name);
   memset(file_name,0,sizeof(file_name));
   strcpy(file_name,aux_fname);*/

   //char Sha_of_Reference[300];
   char *Sha_of_Reference=(char*) malloc(300*sizeof(char));

  // char Sha_of_SignedInfo[300];
   char *Sha_of_SignedInfo=(char*) malloc(300*sizeof(char));



  // char Reference_canonilized[5000];// variable to hold canonicalized reference section
   char *Reference_canonilized=(char*) malloc(5000*sizeof(char));

 //  char xmlExeclusive[5000];
   char *xmlExeclusive=(char*) malloc(5000*sizeof(char));

  // char output[5000];
    char *output=(char*) malloc(5000*sizeof(char));

        // This function will extract the reference section out of the xml file (PA)
	    //and willl perform canonicalization step over it
        Reference_canon(file_name,Reference_canonilized);
        cleanerXML(Reference_canonilized,output);
     //   printf("\noutput:::::::::::::: %s\n",output);
        xmlExeclusiveCanon(output,xmlExeclusive);
        free(Reference_canonilized);
        free(output);

        printf("\n\n execlusive xml  :%s\n\n",xmlExeclusive);

      //  char SignedInfo_canonilized[4000];
      char *SignedInfo_canonilized=(char*) malloc(1000*sizeof(char));
      //char outputS[3000];
    char *outputS=(char*) malloc(3000*sizeof(char));
	    //This section will extract the signed info section out of the xml file (PA)
	    //and will perform canonicalization step over it
	    SignedInfo_canon(file_name,SignedInfo_canonilized);
        cleanerXML(SignedInfo_canonilized,outputS);
        printf("\nSigned Info :%s\n",SignedInfo_canonilized);
        free(SignedInfo_canonilized);

        printf("\n");


	    // SHA value calculated for extracted canonicalized SIgned info and reference section
	    //start
        char ch;
      //  unsigned char st_arr[3000];
	    unsigned char *st=(unsigned char*) malloc(3000*sizeof(unsigned char));
        unsigned int iff=0;
        // assigning char to unsigned char (this is done to have asci value in 8 bits rather than 7 bits)
      while(1){
	      ch=xmlExeclusive[iff];
	      if(ch=='\0'){
		      break;
	      }
         *(st+iff)=ch;
	      //	printf("%c",*(st+iff));
		   iff++;
      }
      free(xmlExeclusive);

	//printf("\n");
	//printf("\n");

        SHA256 sha;

	sha.update(st,iff);

	uint8_t * digest = sha.digest();
   free(st);

   int it=0;//just an iterator
	//printing digest
	printf("\n");

        while(it<32){
		printf("%02x",*(digest+it));
		it++;
	}
       int aux0;
       int count_digest=0;
       it=0;
       while(it<32)
          {

            aux0=(*(digest+it)<<24)|(*(digest+it+1)<<16)|(*(digest+it+2)<<8)|*(digest+it+3);

            for(int h=28;h>=0;h=h-4){
               Sha_of_Reference[count_digest] =HEX_DIGITS[(aux0>>h)&0xf];
               count_digest++;
            }

		    it=it+4;
	   }

      Sha_of_Reference[count_digest]='\0';
   printf("\nSha of reference section that has been canonicalized :%s\n",Sha_of_Reference);

   //   printf("\n");


   char ch1;
  // unsigned char st1_arr[3000];
   unsigned char *st1=(unsigned char*) malloc(3000*sizeof(unsigned char));
   unsigned int i1=0;
    //// assigning char to unsigned char
    while(1){
	   ch1=outputS[i1];
	   if(ch1=='\0'){
		   break;
	   }
        *(st1+i1)=ch1;
	//	printf("%c",*(st1+i1));
		i1++;

    }
   free(outputS);

	printf("\n");
	printf("\n");

   SHA256 sha1;

	sha1.update(st1,i1);

	uint8_t * digest1 = sha1.digest();
   free(st1);
   it=0;//just an iterator
	//printing digest
	printf("\n");

   while(it<32){
		printf("%02x",*(digest1+it));
		it++;
	}
   int aux1;
   count_digest=0;
   it=0;
   while(it<32)
      {

            aux1=(*(digest1+it)<<24)|(*(digest1+it+1)<<16)|(*(digest1+it+2)<<8)|*(digest1+it+3);

            for(int h=28;h>=0;h=h-4){
               Sha_of_SignedInfo[count_digest] =HEX_DIGITS[(aux1>>h)&0xf];
               count_digest++;
            }

		    it=it+4;
	   }

   Sha_of_SignedInfo[count_digest]='\0';
   printf("\nsha of SignedInfo section that has been canonicalized :%s\n",Sha_of_SignedInfo);

//	printf("\n");


   // SHA value calculated for extracted canonicalized SIgned info and reference section
   // end

   //In coming functions :
   //1) extracting digest value
   //2)converting it to hex from base64
   //same above two steps for Signature value.

    // start
    // storing Digest value
   getTagvalue(Tag_Digest_value0,Digest_Value0,file_name);
   //  printf("\nhello  %d\n",Digestv.len);


  // printf("\n Digest value in the permission artefact is %s\n",Digest_Value0);

  // printf("\n");
   base64decoder(Digest_Value0, Digest_Value_in_hex);

   free(Digest_Value0);
  // printf("\n");
  // printf("\n Digest value in hex format %s \n",Digest_Value_in_hex);


   // storing Signed value
   getTagvalue(Tag_Signed_Value0,Signature_Value0, file_name);



   printf("\n Signature value in the permission artefact is%d\n %s\n",int(strlen(Signature_Value0)),Signature_Value0);
   printf("\n");

   base64decoder(Signature_Value0, Signature_Value_in_hex);
   free(Signature_Value0);
   printf("\n");
   printf("\n Signature value in hex format %s \n",Signature_Value_in_hex);

   //end

   // message: bignum object created for storing value of Signature value
   // now Performing modulus operation on bignumber integers

   mp_int message,modulus,public_key,Decrypted;
   mp_init_multi(&message,&modulus,&public_key,&Decrypted,NULL);
   mp_read_radix(&message,Signature_Value_in_hex,16);
   free(Signature_Value_in_hex);

   /// Public key of dgca: This has to be taken from a reserved file in directory./ firmware update

 //  mp_int ;
    char Modulus[513]="d6912a773335d3a193c742071762794e26bcac49d78b6b65a784e1c18d18f2f88b9f13d7fc41a5b9e11e3d75905b0f8373dbac5658962940d104711467814b7319774014644df82e9764b4e852cb7547d56de3224aada000db231b1356ffe86391b3eca2ddf67d44f40ea6e84092cc67387db3a2042487b8dacfe5b588738973a59f5fa813a33e4bb8fbafa407794db990a307201b0a0fbd92ddd181a868a6cfa64625353714e9b1de6f81f3addbf5b202030dbd9e51385db9314591f01af01f07cc92f18ebe60545cea53eb54438e251fd88e7380b5a0612c29ccb7dfd88f7a7d7f0dd07a9b596923602798ad9f1bbb89fffd3cdb8fb94c48640d27fa681e47";
  // char Modulus[513]="ab9d5c8d1fe67207749d63b7dcedd233ce32bb70d175a1bc38c612ab33e2c58e51f83f2788e4d52d9bceb5a1513929de3f526650071a067e6c161b05c60a495fc3ba79ed26f4fa8b2fe2ca8dec44b39759f39206f06a85f9424005a29f05e4cf3a0239340c28c993c1a61cf1b2b6b57c7d8e576ae86827f812b327625baec9ecbf55f1651d35600b9f955f6c2f3bea3aa5852ecdd36a0af818c19acc1030979bed3c89993faa92e0aa0502413b3ca86bbf63477f12ac069aff7137cb72c57f886da79033bbb3b4df0f6cc7fcc18e343aa76036681a566311e267c03b65c98abc91e58f090020c67f776199c0eb76d7e6363687475d3da36ff050f85275607fdd";
   mp_read_radix(&modulus,Modulus,16);

  // mp_int ;
   mp_read_radix(&public_key,"65537",10);

   //mp_int ;
   mp_exptmod(&message, &public_key,&modulus,&Decrypted);   //                       this part for decrypting encrypted text

   char *message_string_hex=(char*) malloc(513*sizeof(char));
   mp_to_hex(&Decrypted,message_string_hex,sizeof(message_string_hex));
   printf("\n Decrypted signature value : %s\n",message_string_hex);//this is the decrypted message
   // At this point message has been decrypted
   // Now having decrypted message in small letters
   mp_clear_multi(&message,&modulus,&public_key,&Decrypted,NULL);
   char padding_SHA256[]="1ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff003031300d060960864801650304020105000420";
   for(int o=0;o<444;o++){
       if(padding_SHA256[o]!=small_letter(message_string_hex[o])){
           return 2;

       }
   }
  // printf("\n\n");
   int k_cout=0;
   for(int i2=445;message_string_hex[i2]!='\0';i2++){
      Extracted_sha_from_Signature_value[k_cout]=small_letter(message_string_hex[i2]);
      k_cout++;
     // printf("%c",message_string_hex[i]);
   }
   free(message_string_hex);
   //with Extracted SHA we are reffering to the expected SHA of SIgned info section
   Extracted_sha_from_Signature_value[k_cout]='\0';//
   printf("\n%s\n",Extracted_sha_from_Signature_value);

   if(strcmp(Extracted_sha_from_Signature_value,Sha_of_SignedInfo)==0){
        printf("SEE if next two string strings are same");
       // printf("\n%s\n ",Sha_of_Reference);
       // printf("%s\n",Digest_Value_in_hex);
        if(strcmp(Sha_of_Reference,Digest_Value_in_hex)==0){
            printf("\n 100000000000PA is valid and non tampered \n");
            int pa_done=1;
            param_set(param_find("PA_CHECK"),&pa_done);
            free(Sha_of_SignedInfo);
            free(Sha_of_Reference);
            free(Extracted_sha_from_Signature_value);
            free(Digest_Value_in_hex);
            return 1;
        }else{
            printf("\n PA is not valid\n");
            free(Sha_of_SignedInfo);
            free(Sha_of_Reference);
            free(Digest_Value_in_hex);
            free(Extracted_sha_from_Signature_value);
            return 0;
        }

   }else{
      printf("\n PA is not valid \n");
    //  PX4_INFO("PA is valid???????????");
      free(Sha_of_SignedInfo);
      free(Sha_of_Reference);
      free(Digest_Value_in_hex);
      free(Extracted_sha_from_Signature_value);
      return 0;
   }
/// at this point validation is complete
}

int Is_PA_VAlid_0(){
	char Tag_Digest_value0[12]="DigestValue";
        char Digest_Value0[100];// in base 64
        char Digest_Value_in_hex[100];// in hex

        char Tag_Signed_Value0[17]="SignatureValue";
        char Signature_Value0[550];// in base 64
        char Signature_Value_in_hex[500];// in hex
        char  Extracted_sha_from_Signature_value[100];

        char file_name[48]="./log/permission_artifact_breach.xml";//"permission_artifact_1.xml";



        char Sha_of_Reference[300];

        char Sha_of_SignedInfo[300];

        char SignedInfo_canonilized[4000];

        char Reference_canonilized[5000];// variable to hold canonicalized reference section

        char xmlExeclusive[5000];
        char output[5000];
        char outputS[3000];
        // This function will extract the reference section out of the xml file (PA)
	    //and willl perform canonicalization step over it
        Reference_canon(file_name,Reference_canonilized);
        cleanerXML(Reference_canonilized,output);
        printf("\noutput:::::::::::::: %s\n",output);
        xmlExeclusiveCanon(output,xmlExeclusive);

        printf("\n\n execlusive xml  :%s\n\n",xmlExeclusive);


	    //This section will extract the signed info section out of the xml file (PA)
	    //and will perform canonicalization step over it
	    SignedInfo_canon(file_name,SignedInfo_canonilized);
        cleanerXML(SignedInfo_canonilized,outputS);
        printf("\nSigned Info :%s\n",SignedInfo_canonilized);

        printf("\n");

      Sha256_implement(xmlExeclusive,Sha_of_Reference,NULL,0);

       Sha256_implement(outputS,Sha_of_SignedInfo,NULL,0);

/*
	    // SHA value calculated for extracted canonicalized SIgned info and reference section
	    //start
        char ch;
        unsigned char st_arr[3000];
	    unsigned char *st=&st_arr[0];
        unsigned int iff=0;
        // assigning char to unsigned char (this is done to have asci value in 8 bits rather than 7 bits)
        while(1){
	      ch=xmlExeclusive[iff];
	     if(ch=='\0'){
		   break;
	     }
         *(st+iff)=ch;
	//	printf("%c",*(st+iff));
		iff++;
        }


	printf("\n");
	printf("\n");

        SHA256 sha;

	sha.update(st,iff);

	uint8_t * digest = sha.digest();
        int it=0;//just an iterator
	//printing digest
	printf("\n");

        while(it<32){
		printf("%02x",*(digest+it));
		it++;
	}
       int aux0;
       int count_digest=0;
       it=0;
       while(it<32)
          {

            aux0=(*(digest+it)<<24)|(*(digest+it+1)<<16)|(*(digest+it+2)<<8)|*(digest+it+3);

            for(int h=28;h>=0;h=h-4){
               Sha_of_Reference[count_digest] =HEX_DIGITS[(aux0>>h)&0xf];
               count_digest++;
            }

		    it=it+4;
	   }

      Sha_of_Reference[count_digest]='\0';
   printf("\nSha of reference section that has been canonicalized :%s\n",Sha_of_Reference);
*/
   //   printf("\n");

/*
   char ch1;
   unsigned char st1_arr[3000];
   unsigned char *st1=&st1_arr[0];
   unsigned int i1=0;
    //// assigning char to unsigned char
    while(1){
	   ch1=outputS[i1];
	   if(ch1=='\0'){
		   break;
	   }
        *(st1+i1)=ch1;
	//	printf("%c",*(st1+i1));
		i1++;

    }


	printf("\n");
	printf("\n");

   SHA256 sha1;

	sha1.update(st1,i1);

	uint8_t * digest1 = sha1.digest();
   it=0;//just an iterator
	//printing digest
	printf("\n");

   while(it<32){
		printf("%02x",*(digest1+it));
		it++;
	}
   int aux1;
   count_digest=0;
   it=0;
   while(it<32)
      {

            aux1=(*(digest1+it)<<24)|(*(digest1+it+1)<<16)|(*(digest1+it+2)<<8)|*(digest1+it+3);

            for(int h=28;h>=0;h=h-4){
               Sha_of_SignedInfo[count_digest] =HEX_DIGITS[(aux1>>h)&0xf];
               count_digest++;
            }

		    it=it+4;
	   }

   Sha_of_SignedInfo[count_digest]='\0';
   printf("\nsha of SignedInfo section that has been canonicalized :%s\n",Sha_of_SignedInfo);
*/
//	printf("\n");


   // SHA value calculated for extracted canonicalized SIgned info and reference section
   // end

   //In coming functions :
   //1) extracting digest value
   //2)converting it to hex from base64
   //same above two steps for Signature value.

    // start
    // storing Digest value
   getTagvalue(Tag_Digest_value0,Digest_Value0,file_name);
   //  printf("\nhello  %d\n",Digestv.len);


  // printf("\n Digest value in the permission artefact is %s\n",Digest_Value0);

  // printf("\n");
   base64decoder(Digest_Value0, Digest_Value_in_hex);
  // printf("\n");
  // printf("\n Digest value in hex format %s \n",Digest_Value_in_hex);


   // storing Signed value
   getTagvalue(Tag_Signed_Value0,Signature_Value0, file_name);



   printf("\n Signature value in the permission artefact is%d\n %s\n",int(strlen(Signature_Value0)),Signature_Value0);
   printf("\n");

   base64decoder(Signature_Value0, Signature_Value_in_hex);
   printf("\n");
   printf("\n Signature value in hex format %s \n",Signature_Value_in_hex);

   //end

   // message: bignum object created for storing value of Signature value
   // now Performing modulus operation on bignumber integers

   mp_int message,modulus,public_key,Decrypted;
   mp_init_multi(&message,&modulus,&public_key,&Decrypted,NULL);
   mp_read_radix(&message,Signature_Value_in_hex,16);

   /// Public key of dgca: This has to be taken from a reserved file in directory./ firmware update

 //  mp_int ;
    char Modulus[513]="d6912a773335d3a193c742071762794e26bcac49d78b6b65a784e1c18d18f2f88b9f13d7fc41a5b9e11e3d75905b0f8373dbac5658962940d104711467814b7319774014644df82e9764b4e852cb7547d56de3224aada000db231b1356ffe86391b3eca2ddf67d44f40ea6e84092cc67387db3a2042487b8dacfe5b588738973a59f5fa813a33e4bb8fbafa407794db990a307201b0a0fbd92ddd181a868a6cfa64625353714e9b1de6f81f3addbf5b202030dbd9e51385db9314591f01af01f07cc92f18ebe60545cea53eb54438e251fd88e7380b5a0612c29ccb7dfd88f7a7d7f0dd07a9b596923602798ad9f1bbb89fffd3cdb8fb94c48640d27fa681e47";
  // char Modulus[513]="ab9d5c8d1fe67207749d63b7dcedd233ce32bb70d175a1bc38c612ab33e2c58e51f83f2788e4d52d9bceb5a1513929de3f526650071a067e6c161b05c60a495fc3ba79ed26f4fa8b2fe2ca8dec44b39759f39206f06a85f9424005a29f05e4cf3a0239340c28c993c1a61cf1b2b6b57c7d8e576ae86827f812b327625baec9ecbf55f1651d35600b9f955f6c2f3bea3aa5852ecdd36a0af818c19acc1030979bed3c89993faa92e0aa0502413b3ca86bbf63477f12ac069aff7137cb72c57f886da79033bbb3b4df0f6cc7fcc18e343aa76036681a566311e267c03b65c98abc91e58f090020c67f776199c0eb76d7e6363687475d3da36ff050f85275607fdd";
   mp_read_radix(&modulus,Modulus,16);

  // mp_int ;
   mp_read_radix(&public_key,"65537",10);

   //mp_int ;
   mp_exptmod(&message, &public_key,&modulus,&Decrypted);   //                       this part for decrypting encrypted text

   char message_string_hex[513];
   mp_to_hex(&Decrypted,message_string_hex,sizeof(message_string_hex));
   printf("\n Decrypted signature value : %s\n",message_string_hex);//this is the decrypted message
   // At this point message has been decrypted
   // Now having decrypted message in small letters
   mp_clear_multi(&message,&modulus,&public_key,&Decrypted,NULL);
   char padding_SHA256[]="1ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff003031300d060960864801650304020105000420";
   for(int o=0;o<444;o++){
       if(padding_SHA256[o]!=small_letter(message_string_hex[o])){
           return 2;

       }
   }
  // printf("\n\n");
   int k_cout=0;
   for(int i2=445;message_string_hex[i2]!='\0';i2++){
      Extracted_sha_from_Signature_value[k_cout]=small_letter(message_string_hex[i2]);
      k_cout++;
     // printf("%c",message_string_hex[i]);
   }
   //with Extracted SHA we are reffering to the expected SHA of SIgned info section
   Extracted_sha_from_Signature_value[k_cout]='\0';//
   printf("\n%s\n",Extracted_sha_from_Signature_value);

   if(strcmp(Extracted_sha_from_Signature_value,Sha_of_SignedInfo)==0){
        printf("SEE if next two string strings are same");
        printf("\n%s\n ",Sha_of_Reference);
        printf("%s\n",Digest_Value_in_hex);
        if(strcmp(Sha_of_Reference,Digest_Value_in_hex)==0){
            printf("\n 100000000000PA is valid and non tampered \n");
            int pa_done=1;
            param_set(param_find("PA_CHECK"),&pa_done);
            return 1;
        }else{
            printf("\n PA is not valid\n");
            return 0;
        }

   }else{
      printf("\n PA is not valid \n");
    //  PX4_INFO("PA is valid???????????");
      return 0;
   }

/// at this point validation is complete

}


char* strip(char *str)
      {
        size_t len = strlen(str);
        memmove(str, str+1, len-2);
        str[len-2] = 0;
        return str;
      }


void updateStack_remove(char *ptr, stack *ss){
    int i=0;
    ss->counter=ss->counter-1;
    while(1){
      char t=ss->list_att[ss->counter][i];
      if (t=='\0'){
         break;
      }else{
         t=' ';
      }
      i=i+1;
    }

}

void updateStack_add(char *ptr,int set, stack *ss)
{
    int i=0;

    while(1){

        if(*(ptr+set+i+1)==' ' || *(ptr+set+i+1)=='>'){

            break;
        }

        ss->list_att[ss->counter][i]=ptr[set+i+1];
        i++;
    }
    ss->list_att[ss->counter][i]='\0';
   // printf("ha %s ha",ss->list_att[ss->counter]);
    ss->counter=ss->counter+1;


}

void replaceTag(char *ptr, stack *ss,char *result,int *ptr_count){
    char aux[40];
    int i=0;
    ss->counter=ss->counter-1;
    while(1){
        char aux_v=ss->list_att[ss->counter][i];
        if (aux_v=='\0' || aux_v==' '){
            if(aux_v=='\0'){
               // printf("endarray");
            }
            if(aux_v==' '){
               // printf("space");
            }
            break;
        }
        aux[i]=aux_v;
        i++;
    }

  //  printf("%c",'>');
    result[*(ptr_count)]='>';
    *(ptr_count)=*(ptr_count)+1;
  //  printf("%c",'<');
    result[*(ptr_count)]='<';
    *(ptr_count)=*(ptr_count)+1;

   // printf("%c",'/');
    result[*(ptr_count)]='/';
    *(ptr_count)=*(ptr_count)+1;


    for(int f=0;f<i;f++){

       // printf("%c",aux[f]);
        result[*(ptr_count)]=aux[f];
        *(ptr_count)=*(ptr_count)+1;

    }
}

int isSignature(char *ptr){
    /* if signature start ::0
       if signature end ::1
       if none  ::2
    */
    char check0[18]="<Signature xmlns=";

    char check1[14]="</Signature>";

    char check2[14]="</SignedInfo>";

    char check3[13]="<SignedInfo>";

    int check0_status=1;

    for(int i=0;i<(int)strlen(check0);i++)
    {
        if(*(ptr+i)!=check0[i])
        {
            check0_status=0;
            break;
        }
    }

    int check1_status=1;

    for(int i=0;i<(int)strlen(check1);i++)
    {
        if(*(ptr+i)!=check1[i])
        {
            check1_status=0;
            break;
        }
    }

    int check2_status=1;
    for(int i=0;i<(int)strlen(check2);i++){
        if(*(ptr+i)!=check2[i]){
            check2_status=0;
            break;
        }
    }


    int check3_status=1;
    for(int i=0;i<(int)strlen(check3);i++){
        if(*(ptr+i)!=check3[i]){
            check3_status=0;
            break;
        }
    }

    if(check0_status==1){
        return 0;
    }
    if(check1_status==1){
        return 1;
    }
    if(check2_status==1){
        return 2;
    }
    if(check3_status==1){
        return 3;
    }
    return 4;

}


int aux_space_count;
void Reference_canon(char *file_name, char *res){


    FILE *fp;
    //char res[5000];
    fp=fopen(file_name,"r");
    if(fp)
    {

	    //exit(1);

    signed char ch;
    int count_result=0;
    int *ptr_count_result;
    ptr_count_result=&count_result;
    int i=0;



    char *content=(char*) malloc(sizeof(char)*5000);
    //char content[5000];
    int stop_variable=0;
   // int stop_variabletype2=0;
    while(1)
    {

        ch=fgetc(fp);
        if(ch==EOF){
            break;
        }
        if(content[i-1]=='<' && ch=='!'){
            stop_variable=1;
            i--;
        }
        if(content[i-1]=='<' && ch=='?'){
            stop_variable=1;
            i--;
        }
        /*if(content[i-1]=='>' && ch==' '){
            stop_variabletype2=1;
            i--;
        }*/
        if(stop_variable==1 && ch=='>'){
            stop_variable=0;
            continue;
        }
      /*  if(stop_variabletype2==1 && ch=='<'){
            stop_variabletype2=0;
            continue;
        }*/
        if(stop_variable!=1 && ch!='\n'){
            content[i]=ch;
            i++;
        }


    }
    stack *dd;
    stack dds;
    dds.counter=0;
    dd=&dds;

    //int comment_start=0;
    int flag0=0,flag1=0;
    int aux_check_flag;
    int space_count=0;
    for(int u=0;u<i;u++)
    {



        if(content[u]=='<')
        {
            aux_check_flag=isSignature(content+u);
           /* if(content[u+1]=='!'){
                comment_start=1;
            }*/
            if(aux_check_flag==0)
                {   aux_space_count=space_count;
                    flag0=1;
                }
            if(aux_check_flag==1){
                flag1=1;
                aux_space_count=space_count;
                u=u+12;

            }
            if(aux_check_flag==2){

            }
        }

        if((flag0==0 && flag1==0) || (flag0==1 && flag1==1) )
        {

            if(content[u]=='<')
            {
                if(content[u+1]!='/')
                {
                    updateStack_add(content,u,dd);
                    space_count=space_count+2;
                }

            }

            if(content[u]=='/' && content[u+1]=='>')
            {

                replaceTag(content+u,dd,res,ptr_count_result);
                space_count=space_count-2;

            }else{
                res[count_result]=content[u];
                count_result=count_result+1;
            }

            if(content[u]=='>' && content[u+1]=='<')
            {


                if(content[u+2]=='/')
                {
                    updateStack_remove(content+u,dd);

                     space_count=space_count-2;
                }

            }

        }


    }
    res[count_result]='\0';
  //  strcpy(res1,res);
    printf("\ncacccccccccccc:::::: %s\n",res);
    free(content);
    }else{
        printf("PA File does not exist");
    }
    fclose(fp);
}
void cleanerXML(char *input,char *output){
    int i=0;
    int o=0;
    int checkspace=0;
    while(input[i]!='\0')
    {
        if(input[i]=='>'){
            checkspace=1;
        }
        if(checkspace==1 && input[i]==' '){
            i++;
            continue;

        }
        if(checkspace==1 && input[i]=='<'){
            checkspace=0;

        }
        output[o]=input[i];
        o++;
        i++;

    }
    output[o]='\0';


}

void SignedInfo_canon(char *file_name, char *res)
{

    FILE *fp;

    fp=fopen(file_name,"r");
    if(fp)
    {


    char *content=(char*) malloc(sizeof(char)*10000);

    signed char ch;
    char aux_copy[300];
    int count=0;
    int flag_ini=0;
    int flag_end=0;
     int *ptr_count_result;

    int res_count=0;

    ptr_count_result=&res_count;

    while(1)
    {

        ch=fgetc(fp);
        if(ch==EOF){
            break;
        }
        if(ch!='\n'){
        content[count]=ch;
        count++;}
    }

    stack dde;
    stack *de;
    de=&dde;
    dde.counter=0;
    int space_count;
    int first_pass_flag=0;
    int in_count=0;//aux length
    int first_time_flag=0;
    for(int i=0;i<count;i++)

    {
        if(content[i]=='<')
        {

        int check_Signature;

        check_Signature=isSignature(content+i);
        if(check_Signature==0)
        {  // printf("\nfound start of signed info\n");
            int k=11;
            while(content[k+i]!='>')
            {
                aux_copy[k-11]=content[k+i];
                k++;
                in_count++;

            }
            aux_copy[k]='\0';
          //  printf("%s\n",aux_copy);
            flag_ini=1;
        }

        if(check_Signature==2)
        {
            //end()
           // printf("\njhghggg\n");
            flag_end=1;
            for(int f=0;f<13;f++)
            {
            res[res_count]=content[i+f];
            res_count++;
            }
        }

        if(check_Signature==3)
        {
           // content[i+12]=' ';
            space_count=aux_space_count;
            for(int o=0;o<11;o++)
            {
                res[res_count]=content[i+o];
              //  printf("%c",res[res_count]);
                res_count++;
            }
            i=i+12;
            res[res_count]=' ';
            res_count++;

            for(int j=0;j<in_count;j++)
            {
                res[res_count]=aux_copy[j];
              //  printf("%c",res[res_count]);
                res_count++;
            }
            res[res_count]='>';
            res_count++;
          //  res[res_count]='\n';
          //  res_count++;

            first_pass_flag=1;
        }
        }

        if(flag_ini==1 && first_pass_flag==1 && flag_end!=1)
        {
            if(first_time_flag==0){
                first_time_flag=1;
                space_count=space_count+4;
                for(int y=0;y<space_count;y++)
                {
                   // res[res_count]=space;
                   // res_count++;
                }
            }

            if(content[i]=='<')
            {
                if(content[i+1]!='/')
                {
                    updateStack_add(content,i,de);
                    space_count=space_count+2;
                }
                if(content[i+1]=='/' && content[i-1]!='>')
                {

                    updateStack_remove(content+i,de);
                   /*  printf("\nprinting stack::\n");
                    for(int i=0;i<de->counter;i++){
                        printf("%s ",de->list_att[i]);
                    }
                    printf("\n\n");*/
                    space_count=space_count-2;
                }
              //  printf("\npassed updatestack\n");

            }

            if(content[i]=='/' && content[i+1]=='>')
                {

                    replaceTag(content+i,de,res,ptr_count_result);
                    space_count=space_count-2;

                }else{

              //  printf("%c",content[i]);
                res[res_count]=content[i];
                res_count=res_count+1;
                }

            if(content[i]=='>' && content[i+1]=='<')
            {
               // printf("\n");//new line

               // res[res_count]=0x0a;//'\n';

               // res_count=res_count+1;


                if(content[i+2]=='/')
                {
                    updateStack_remove(content+i,de);
                  /*  printf("\nprinting stack::\n");
                    for(int i=0;i<de->counter;i++){
                        printf("%s ",de->list_att[i]);
                    }
                    printf("\n\n");*/
                     space_count=space_count-2;
                }





            }


        }

    }
    res[res_count]='\0';
    free(content);
    fclose(fp);
    }
    else{
        printf("PA file does not exist");

    }



}

//function get tag value
void getTagvalue(char *Tag, char *certificate,char *file_nam_perm){
    //certificate is the value of tag.
   FILE *fp;
   fp=fopen(file_nam_perm,"r");
   // Digest DD;

    char *buf=(char*) malloc(1000*sizeof(char));// earlier it was static defined [3000]

    char tagi[30];
    char tage[30];

    memset(tagi, 0, sizeof(tagi));
    strcpy(tagi, "<");
    strcat(tagi,Tag );
    strcat(tagi, ">");

    memset(tage, 0, sizeof(tage));
    strcpy(tage, "</");
    strcat(tage,Tag );
    strcat(tage, ">");

    char *ptr1,*ptr2;

    int c_flag_i=0,c_flag_e=0,line=0;

    int index[2][2]={{0,0},{0,0}};

    int  aux=0;

    int c=0;

   while(fscanf(fp, "%s", buf) != EOF )
    {
    line=line+1;
    if(isSubstring(tagi,buf)!=-1 || isSubstring(tage,buf)!=-1 ){
      //  printf(" %d ",line);
        if (isSubstring(tagi,buf)!=-1){       //for "<X509Certificate>"
          // printf("start");
           index[0][0]=line;
           index[0][1]=isSubstring(tagi,buf);
           c_flag_i=1;
           aux=1;

        }
        if(isSubstring(tage,buf)!=-1){                          //for "</X509Certificate>"
          // printf("  end");
           index[1][0]=line;
           index[1][1]=isSubstring(tage,buf);
          // printf("\n%d   %d\n",index[1][0],index[1][1]);
           c_flag_e=1;


        }
        }

        if(c_flag_i==1 && c_flag_e==0 ){
        if (aux==1){

        for(int u=index[0][1]+strlen(tagi);u<(int)strlen(buf);u++){
            certificate[c]=buf[u];

            c++;
        }

        aux=0;
        }
        else{ptr1=certificate+c;
            ptr2=&(buf[0]);
            memcpy(ptr1,ptr2,strlen(buf));
            c=c+strlen(buf);
        }
        }
        else if(c_flag_i==1 && c_flag_e==1 )
        {
            if(index[0][0]==index[1][0])//start line and end line is same
        {
            for(int u=index[0][1]+strlen(tagi);u<index[1][1];u++)
            {
            certificate[c]=buf[u];
            c++;
            }

        }else
        {  //end line is different than start line
           for(int u=0;u<index[1][1];u++)
            {
            certificate[c]=buf[u];
            c++;
            }

        }
        break;
        }else{continue;}

    }
    certificate[c]='\0';
    free(buf);
    fclose(fp);

}



////////////////////////

//function for maintaining hardware integrity
//return types
//0:not genuine harware used
//1 : genuine hardware used
int RPAS_identifier(FILE *fptr1){

    // taking gps id from the attached gps sensor
	 int vehi_gps=orb_subscribe(ORB_ID(sensor_gps));
    struct sensor_gps_s raw;
    orb_copy(ORB_ID(sensor_gps),vehi_gps,&raw);
    int deviceID=raw.device_id;

	fprintf(fptr1,"Device ID of GPS is: %d ",deviceID);

    deviceID=2348919;
    // reading device id of gps from HardwareInuse.txt

    char fname_0[100]="/fs/microsd/log/HardwareInuse.txt";
    /*char aux_fname[100];
    strcpy(aux_fname,dir);
    strcat(aux_fname,fname_0);
    memset(fname_0,0,sizeof(fname_0));
    strcpy(fname_0,aux_fname);*/

    char tag_0[30]="GPS_ID";
    char GPS_ID_file[30];

    // this function first validates the file (if its written inside rfm) and then
    // fetches the required tag
    int status_gps=file_read(fname_0,tag_0,GPS_ID_file,0);// file name

    if (status_gps==1){
        // file is valid
        // now compare the gps id of attached hardware to the one in the file
        char aux_gps_id[30];

        sprintf(aux_gps_id,"%d",deviceID);

        int gps_id_status=strcmp(aux_gps_id,GPS_ID_file);

        if (gps_id_status==0){

            // both strings matches and hardware attached is genuine
            return 1;// all okay: hardware is genuine

        }else{
            // possiblity of new authentic hardware or unauthentic hardware
            //check for HardwareChange.txt file
            char Hardware_fname[130]="/fs/microsd/log/HardwareChange.txt" ;
            FILE *fptr;
            fptr=fopen("/fs/microsd/log/HardwareChange.txt","r");
            if (fptr==NULL){
               // printf("\n\nhello inside RPAS identifier\n\n");
                // no such file exists , it is an attempt to use unautherised hardware
            //    write_log("UH");//UH =attemp to use unauthorized hardware
                return 0;
            }else{
                //file has been found, now validate it using FMP public key
                fclose(fptr);
                char id_gps_change[20];
                int file_status = file_read(Hardware_fname,tag_0,id_gps_change,1);
                if (file_status==1){
                   //file is valid
                   //compare strings in id_gps_change and aux_gps_id
                   if(!strcmp(aux_gps_id,id_gps_change)){
                       //the strings are same
                       //start the modification of HardwareInuse.txt
                       remove("/fs/microsd/log/HardwareInuse.txt");
                       HardwareInuseCreation(deviceID);//new HardwareInuse.txt file created
                       // and return 1
                       return 1;
                   }else{
                       // id written in the file (HArdwareCHange.txt ) is different from the one
                       // retrieved from attached hardware.
                       // different hardware is attached to the RPAS
                       return 0;

                   }

                }else
                {
                   //file is not valid
                   // case where file exists but is not valid(not sgned by MS)
                   return 0;

                }

            }

        }

    }
    else{
        //file is not valid
        // firmware update is required
        // No genuine HardwareInuse.txt found, needs factory reset
        //write_log for firmware update.
        return 2;
    }




}



void   call_DroneIDcreation(){
    DroneIDcreation( );
}



int call_file_read(char *fname,char *tag,char *result,int file_type){
    printf("\ninside call file read \n");
    int return_value=file_read(fname,tag,result,file_type);
    return return_value;

}



int CHECK_REUSAGE(char *file_id){
    //check the given id in file KeyLog.txt
    //return types 1=file_id is beng reused, dont do key rotation
    //             0=file_id is fresh , start key rotation and amend KeyLog.txt
    char fname[100]="/fs/microsd/log/KeyLog.txt";
    FILE *fptr;
    fptr=fopen(fname,"r");
    if(fptr!=NULL){
        fclose(fptr);
        key RFM_key;
       // strcpy(RFM_key.modulus,"26495127767604337655831691918113833077956334480395285962585476150242592943334533493300410596845191511956515492348526447864429111984398076478393159921716146324008099192785429713391158778980470576914436564260777828713503499078118282954851831718742193757807790192236071369965695862117529898068098551948463482032984966019192036825055149445334134853954472696975028748893285534653718791714099153590644093991103487859923770318942338310035808291238607889065040628342091199251953687656515817829870399170890472577591915379199888970387732102895345522222635522973211224211091186930411354449072698100145010326153546937615803429429");
       // strcpy(RFM_key.private_exponent,"24021354375322558781058142276736617999389802434596138079632937453577588041977071136990155131350955784632074057774459381406055342415566243621056270140966629267130220147961070276489248399064064585487617556117107847879807603464052857723291974564966716033712609329726458163489658005940146657360121454440603066594695330718634482791509166056105316975123379501455703860598165895993179884232993174716442550760635024519308123321406513256918808853580242186763289628649427009957647963922327454796638500088631677702805825426368000585849292006901367973110968487920147218164610608698409129440488923014780703396073978252084600942649");
        get_RFM_Key(&RFM_key);
       int valid_status= Validating_File(fname,RFM_key);

       if (valid_status==1){
           //file is valid it can searched for current file_id
            FILE *fptr1;
            fptr1=fopen(fname,"r");
            signed char aux_c;
            char *content=(char*) malloc(sizeof(char)*3000);
            int i=0;
            while((aux_c=fgetc(fptr1))!=EOF){
                content[i]=aux_c;
                i++;
            }
            content[i]='\0';
            fclose(fptr1);
           int prescence=isSubstring(file_id,content);
           free(content);
           if(prescence!=-1){
               return 0;
           }else{
               return 1;
           }
       }else{
           return 1;
       }
    }else{
        return 2;
    }


}



//for calling KeyLog_Regen
void call_KeyLog_Regen(char *fileID){
    KeyLog_Regen(fileID);
}


// Key rotation function
//two files are generated
//1)PublicKeyNew.txt: later sent to MS
//2)PublicPrivateNew.txt:kept inside RFM encrypted
void Key_rotation_start(char *file_id){
  /// NOW prime number generation begins for KEY pair creation
    mp_int  p1, q1;
    char Prime1[1000];
    char Prime2[1000];
    char    buf[4096];
    char PrivateKey[1000];
    char Modulus[1000];
    char PublicKey[20];

    //mp_int z,r;
    mp_init(&p1);
    mp_init(&q1);

    // mp_prime_rand(&p1,1,90,MP_PRIME_SAFE);
    mp_prime_rand(&p1,1,1024,MP_PRIME_2MSB_ON);
    mp_to_decimal(&p1,Prime1,sizeof(Prime1));
    //printf("\nprime1 == %s \n",buf);

    printf("\n Key Generation...... \n");
    //mp_prime_rand(&q1,1,90,MP_PRIME_SAFE);
    mp_prime_rand(&q1,1,1024,MP_PRIME_2MSB_ON);
    mp_to_decimal(&q1,Prime2,sizeof(Prime2));
    //mp_clear(&p1);
   // mp_clear(&q1);

   // printf("\n prime2 == %s \n",buf);
    printf("\n Key Generation...... \n");
    //phi(n)=(p1-1)*(q1-1)
    mp_int s1,s2,product_p1_q1;
    mp_init_multi(&s1,&s2,&product_p1_q1,NULL);



    mp_mul(&p1,&q1,&product_p1_q1);//modulus n=p1*q1

    mp_sub_d(&p1, 1uL, &s1);
  //  mp_to_decimal(&s1, buf, sizeof(buf));
  //  printf("\n\ns1 == %s\n", buf);//p1-1
    printf("\n Key Generation...... \n");
    mp_sub_d(&q1, 1uL, &s2);
    mp_to_decimal(&s2, buf, sizeof(buf));
   // printf("\n\ns2 == %s\n", buf);//q1-1
    printf("\n Key Generation...... \n");
    mp_int product_s1_s2,e,d;
    mp_init_multi(&e,&d,&product_s1_s2,NULL);
    mp_mul(&s1,&s2,&product_s1_s2);// phi(n)=(p1-1)*(q1-1)
    //mp_int e;
    mp_read_radix(&e,"65537",10);//
    //  mp_int d;
    mp_invmod(&e,&product_s1_s2,&d);

    //writing public key
    mp_to_decimal(&e,PublicKey,sizeof(PublicKey));
   // printf("\n Public key : e===\n%s\n\n",buf);
    printf("\n Key Generation...... \n");



    ///writing private key
    mp_to_decimal(&d,PrivateKey,sizeof(PrivateKey));
  //  printf("\n Private key : d==\n%s\n\n",buf);
    printf("\n Key Generation...... \n");





    mp_to_decimal(&p1,buf,sizeof(buf));
  //  printf("\n prime number 1: p1==\n%s\n\n",buf);
    mp_to_decimal(&q1,buf,sizeof(buf));
  //  printf("\n prime number 2: q1==\n%s\n\n",buf);



    //writing modulus in file
    mp_to_decimal(&product_p1_q1,Modulus,sizeof(Modulus));
  //  printf("\n modulus product==\n%s\n\n",buf);
    printf("\n Key Generation...... \n");

    mp_clear(&d);
    mp_clear(&p1);
    mp_clear(&q1);
    mp_clear(&e);
    mp_clear(&product_s1_s2);
    mp_clear(&s1);
    mp_clear(&s2);
    mp_clear(&product_p1_q1);



    char DroneID[]="ABBDDJEDNDJK";//any string
    char RFM_version[2]="0";
    char RPAS_category[10]="Small";
    pair_set pairset0[7];
    pair_set pairset1[6];

    strcpy(pairset0[0].tag,"DroneID");
    strcpy(pairset0[0].value,DroneID);

    strcpy(pairset1[0].tag,"DroneID");
    strcpy(pairset1[0].value,DroneID);


    strcpy(pairset0[1].tag,"RFM_version");
    strcpy(pairset0[1].value,RFM_version);

    strcpy(pairset1[1].tag,"RFM_version");
    strcpy(pairset1[1].value,RFM_version);

    strcpy(pairset0[2].tag,"RPAS_category");
    strcpy(pairset0[2].value,RPAS_category);

    strcpy(pairset1[2].tag,"RPAS_category");
    strcpy(pairset1[2].value,RPAS_category);


    strcpy(pairset0[3].tag,"PublicKey");
    strcpy(pairset0[3].value,PublicKey);

    strcpy(pairset1[3].tag,"PublicKey");
    strcpy(pairset1[3].value,PublicKey);

    strcpy(pairset0[4].tag,"PrivateKey");
    strcpy(pairset0[4].value,PrivateKey);

    strcpy(pairset0[5].tag,"Modulus");
    strcpy(pairset0[5].value,Modulus);

    strcpy(pairset1[4].tag,"Modulus");
    strcpy(pairset1[4].value,Modulus);


    strcpy(pairset0[6].tag,"KEY_ID");
    strcpy(pairset0[6].value,file_id);

    strcpy(pairset1[5].tag,"KEY_ID");
    strcpy(pairset1[5].value,file_id);


    char fileName[50]="/fs/microsd/log/PublicPrivateNew_aux.txt";

    char filename1[40]="/fs/microsd/log/PublicKeyNew.txt";

    key RFM_key;
    get_RFM_Key(&RFM_key);


    pair_file_write(pairset1,6,filename1,RFM_key);
    pair_file_write(pairset0,7,fileName,RFM_key);//this needs to be encrypted

    // Now signed file has been made
    FILE *fptr;
    fptr=fopen("/fs/microsd/log/PublicPrivateNew_aux.txt","r");

    signed char ch;
    char content[3000];//=(char*) malloc(sizeof(char)*3000);
    int i=0;
    while((ch=fgetc(fptr))!=EOF){
        content[i]=ch;
        i++;
    }
    content[i]='\0';
    fclose(fptr);
    remove("/fs/microsd/log/PublicPrivateNew_aux.txt");
    key Inside_RFM;

    char fname[42]="/fs/microsd/log/PublicPrivateNew.txt";
    strcpy(Inside_RFM.modulus,"180919775566931");
    strcpy(Inside_RFM.private_exponent,"32102716896161");
    encrypting_File(content,Inside_RFM,fname);
    char decrypting[3000];
    decrypting_File(fname,Inside_RFM,decrypting);
    printf("\n\n%s\n\n",decrypting);
  //  free(content);
    // new key is encrypted and saved inside the file





/// KEY pair generation ends here

}

/*In this function contents of PublicPrivateNew.txt are fetched. The corresponding tags and values
inside the PublicPrivateInuse.txt are replaced by the ones that are fetched earlier.
The process invloves 1)decrypting of PublicPrivateNew.txt
                     2)fetching Modulus, Privatekey, PublicKey
                     3)making new PublicPrivateInuse.txt
*/
void KEY_CHANGE_INITIATION(){

    char fname[60]="/fs/microsd/log/PublicPrivateNew.txt";
    key Inside_RFM;
    strcpy(Inside_RFM.modulus,"180919775566931");
    strcpy(Inside_RFM.private_exponent,"32102716896161");

    char content[5000];//=(char*) malloc(sizeof(char)*5000);

    decrypting_File(fname , Inside_RFM, content);
   // printf("\n\nThe content of decrypted file is :\n\n%s\n",content);

    char target0[30]="Modulus";
    char target1[30]="PrivateKey";

    char value_modulus[2000];//=(char*) malloc(sizeof(char)*2000);
    char value_private[800];//=(char*) malloc(sizeof(char)*800);

    fetch_tag(content,target0,value_modulus);
    fetch_tag(content,target1,value_private);

    /*pair_set pair_key[2];


    strcpy(pair_key[0].tag,"Modulus");
    strcpy(pair_key[0].value,value_modulus);

    strcpy(pair_key[1].tag,"Privatekey");
    strcpy(pair_key[1].value,value_private);*/


    char content1[5000];

    char Content_Start[10]="<content>";
    strcpy(content1,Content_Start);
    strcat(content1,"\n");

    char tag_private[20]="<PrivateKey=";
    strcat(content1,tag_private);
    strcat(content1,value_private);
    strcat(content,">\n");

    char tag_modulus[20]="<Modulus=";
    strcat(content1,tag_modulus);
    strcat(content1,value_modulus);
    strcat(content1,">\n");


    char Content_End[12]="</content>";
    strcat(content1,Content_End);
    strcat(content1,"\0");
    //remove("./log/PublicPrivateInuse.txt");
    char key_file_name[40]="/fs/microsd/log/PublicPrivateInuse.txt";
    encrypting_File(content1,Inside_RFM,key_file_name);
    // at this point new PublicPrivateInuse.txt is formed
    remove("/fs/microsd/log/PublicPrivateNew.txt");
    remove("/fs/microsd/log/KeyChangePerm.txt");
    remove("/fs/microsd/log/PublicKeyNew.txt");



}


/*This function comes in handy when ParamChangePerm.txt is present and one
has to modify already existig ParamInuse.txt file
read ParamChangePerm.txt
if any similar then change the value
if any new then add the value
*/
void ParamInuseModify(){

     // first store the information inside ParamInuse.txt
    FILE *fptr_INuse;
    fptr_INuse=fopen("/fs/microsd/log/ParamInuse.txt","r");

    signed char ch;
    char ParamInuseContent[5000];
    int i=0;
    while((ch=fgetc(fptr_INuse))!=EOF){
        ParamInuseContent[i]=ch;
        i++;
    }
    ParamInuseContent[i]='\0';
    fclose(fptr_INuse);


    //checking in ParamChangePerm
    FILE *fptrParamChange;
    fptrParamChange=fopen("/fs/microsd/log/ParamChangePerm.txt","r");

    char ParamChangeContent[4000];
    i=0;
    while((ch=fgetc(fptrParamChange))!=EOF){
        ParamChangeContent[i]=ch;
        i++;
    }
    ParamChangeContent[i]='\0';
    fclose(fptrParamChange);



    char no_params_tag[20]="No_Params";
    char No_params[10];
    fetch_tag(ParamInuseContent,no_params_tag,No_params);

    int index=isSubstring(no_params_tag,ParamInuseContent);

    while(ParamInuseContent[index]!='\n'){
        index++;
    }
    index++;

   // printf("\ncharacter kfdlkfsldk %c\n",ParamInuseContent[index]);


  //  char no_params_tag[20]="No_Params";
    char No_params_change[10];
    fetch_tag(ParamChangeContent,no_params_tag,No_params_change);
    int no_params_change=atoi(No_params_change);

    int no_params=atoi(No_params);
    printf("dff %d",no_params);



    int newline_count=0;
    int start_param=0;
    int start_value=0;
    char buff_param[100];
    char buff_value[100];
    int buff_it=0;
    int jump_value=0;
    //pair *param_value=(pair*) malloc(no)
    pair_set param_value[400];
    int param_value_count=1;
    printf("\n%s\n",ParamInuseContent);
    while(newline_count<no_params)
    {

       // printf("\ncharacter inside 000 %c\n",ParamInuseContent[index]);

        if(ParamInuseContent[index]=='>'){
          //  printf("\ninsiid\n");
            newline_count++;
            index=index+2;
            start_value=0;
            buff_value[buff_it]='\0';
            buff_it=0;
            if(jump_value!=1){
                strcpy(param_value[param_value_count].tag,buff_param);
                strcpy(param_value[param_value_count].value,buff_value);
                param_value_count++;


            }else{
                jump_value=0;
            }
        }

        if(start_value==1 && jump_value==0){
            buff_value[buff_it]=ParamInuseContent[index];
            buff_it++;
        }


        if(ParamInuseContent[index]=='='){
            start_param=0;
            start_value=1;
            buff_param[buff_it]='\0';
            buff_it=0;
            int check_availability=isSubstring(buff_param,ParamChangeContent);
            if(check_availability!=-1){
                // tag is available in ParamChange also
                // fetch the value from ParamChangeContent

                fetch_tag(ParamChangeContent,buff_param,buff_value);
                strcpy(param_value[param_value_count].tag,buff_param);
                strcpy(param_value[param_value_count].value,buff_value);
                param_value_count++;

                jump_value=1;
            }
        }



        if(start_param==1){
            buff_param[buff_it]=ParamInuseContent[index];
           // printf("\n adcdsnfdfds this is char %c\n",buff_param[buff_it]);
            buff_it++;
        }


        if(ParamInuseContent[index]=='<'){
            start_param=1;
        }

      //  printf("\ncharacter inside %c\n",ParamInuseContent[index]);
        index++;


    }



    start_param=0;
    start_value=0;
    buff_it=0;
    jump_value=0;

    index=isSubstring(no_params_tag,ParamChangeContent);

    while(ParamChangeContent[index]!='\n'){
        index++;
    }
    index++;
    newline_count=0;
    while(newline_count<no_params_change){

            if(ParamChangeContent[index]=='>'){
                newline_count++;
                index=index+2;
                start_value=0;
                buff_value[buff_it]='\0';
                buff_it=0;
            if(jump_value!=1){
                strcpy(param_value[param_value_count].tag,buff_param);
                strcpy(param_value[param_value_count].value,buff_value);
                param_value_count++;

            }else{
                jump_value=0;
            }
            }
            if(start_value==1 && jump_value==0){
                buff_value[buff_it]=ParamChangeContent[index];
                buff_it++;
            }
        if(ParamChangeContent[index]=='='){
            start_param=0;
            start_value=1;
            buff_param[buff_it]='\0';
            buff_it=0;
            printf("\n helloooooooo %s  \n",buff_param);
            int check_availability=isSubstring(buff_param,ParamInuseContent);
            if(check_availability==-1){
                // tag is available in ParamInuse also
                // no need to add
                //fetch_tag(ParamChangeContent,buff_param,buff_value);
                //strcpy(param_value[param_value_count].tag,buff_param);
                //strcpy(param_value[param_value_count].value,buff_value);
                //param_value_count++;

            }else{
                jump_value=1;
            }
        }

        if(start_param==1){
            buff_param[buff_it]=ParamChangeContent[index];
            buff_it++;
        }

        if(ParamChangeContent[index]=='<'){
            start_param=1;
        }


        index++;



    }
    char aux_int_no[20];

    sprintf(aux_int_no,"%d",param_value_count);
    strcpy(param_value[0].tag,"No_Params");
    strcpy(param_value[0].value,aux_int_no);

    char file_name_inuse[40]="/fs/microsd/log/ParamInuseNew.txt";
    key RFM_key;

    get_RFM_Key(&RFM_key);

    pair_file_write(param_value,param_value_count,file_name_inuse,RFM_key);

    remove("/fs/microsd/log/ParamChangePerm.txt");

}


/*
function to start setting the parameters as specified in ParamInuse.txt
*/
void ParamSetfile(){
     // first store the information inside ParamInuse.txt
    FILE *fptr_INuse;
    fptr_INuse=fopen("/fs/microsd/log/ParamInuse.txt","r");

    signed char ch;
    char ParamInuseContent[5000];
    int i=0;
    while((ch=fgetc(fptr_INuse))!=EOF){
        ParamInuseContent[i]=ch;
        i++;
    }
    ParamInuseContent[i]='\0';
    fclose(fptr_INuse);

    char no_params_tag[20]="No_Params";
    char No_params[10];
    fetch_tag(ParamInuseContent,no_params_tag,No_params);

    int index=isSubstring(no_params_tag,ParamInuseContent);

    while(ParamInuseContent[index]!='\n'){
        index++;
    }
    index++;

    int no_params=atoi(No_params);
   // printf("dff %d",no_params);



    int newline_count=0;
    int start_param=0;
    int start_value=0;
    char buff_param[100];
    char buff_value[100];
    int buff_it=0;
   // int jump_value=0;
    //pair *param_value=(pair*) malloc(no)

    printf("\n%s\n",ParamInuseContent);
    while(newline_count<no_params)
    {

      //  printf("\ncharacter inside 000 %c\n",ParamInuseContent[index]);

        if(ParamInuseContent[index]=='>'){
            //printf("\ninsiid\n");
            newline_count++;
            index=index+2;
            start_value=0;
            buff_value[buff_it]='\0';
            buff_it=0;


            int gfint=atoi(buff_value);

            param_set(param_find(buff_param),&gfint);

        }

        if(start_value==1 ){
            buff_value[buff_it]=ParamInuseContent[index];
            buff_it++;
        }


        if(ParamInuseContent[index]=='='){
            start_param=0;
            start_value=1;
            buff_param[buff_it]='\0';
            buff_it=0;

        }



        if(start_param==1){
            buff_param[buff_it]=ParamInuseContent[index];
          //  printf("\n adcdsnfdfds this is char %c\n",buff_param[buff_it]);
            buff_it++;
        }


        if(ParamInuseContent[index]=='<'){
            start_param=1;
        }

       // printf("\ncharacter inside %c\n",ParamInuseContent[index]);
        index++;


    }


}




//this function will make the amendments in recentPA.txt file
/* type == 0,1,2;
fetch_required    0
freq_done         1
previous_log_hash 2
*/
void update_recentPA(int type,char *value){
   FILE *fptr;

    fptr=fopen("/fs/microsd/log/recentPA.txt","r");
    if(fptr==NULL){
       printf("\nproblem in opening file.\n");
    }
    signed char ch;
    int i=0;

    char content[10000];

    while((ch=fgetc(fptr))!=EOF){
        content[i]=ch;
        i++;
    }
    content[i]='\0';
   printf("\ndata_fetched\n");
   fclose(fptr);
   // type == 0,1,2;
   // first fetching all the constant tags
 /*permissionArtifactId
   frequency_allowable
   maxAltitude
   Start_time
   End_time
   lattitude1
   longitude1
   lattitude2
   longitude2*/
   pair_set pairset[14];

   //permissionArtifactId
   char tag_PAid[70]="permissionArtifactId";
   char tag_content_paID[100];
   fetch_tag(content,tag_PAid,tag_content_paID);
   strcpy(pairset[0].tag,tag_PAid);
   strcpy(pairset[0].value,tag_content_paID);

   // frequency allowable
   char tag_freq[70]="frequency_allowable";
   char tag_content_freq[100];
   fetch_tag(content,tag_freq,tag_content_freq);
   strcpy(pairset[1].tag,tag_freq);
   strcpy(pairset[1].value,tag_content_freq);

   //    maxAltitude
   char tag_max[70]="maxAltitude";
   char tag_content_max[100];
   fetch_tag(content,tag_max,tag_content_max);
   strcpy(pairset[2].tag,tag_max);
   strcpy(pairset[2].value,tag_content_max);



    //// start time
   char tag_st[70]="Start_time";
   char tag_content_st[100];
   fetch_tag(content,tag_st,tag_content_st);
   strcpy(pairset[3].tag,tag_st);
   strcpy(pairset[3].value,tag_content_st);


    //// end time
   char tag_et[70]="End_time";
   char tag_content_et[100];
   fetch_tag(content,tag_et,tag_content_et);
   strcpy(pairset[4].tag,tag_et);
   strcpy(pairset[4].value,tag_content_et);

   printf("\nproblem at 9950\n");
   //lattitude1
   char tag_lat1[70]="lattitude1";
   char tag_content_lat1[100];
   fetch_tag(content,tag_lat1,tag_content_lat1);
   strcpy(pairset[5].tag,tag_lat1);
   strcpy(pairset[5].value,tag_content_lat1);


   //longitude1
   char tag_lon1[70]="longitude1";
   char tag_content_lon1[100];
   fetch_tag(content,tag_lon1,tag_content_lon1);
   strcpy(pairset[6].tag,tag_lon1);
   strcpy(pairset[6].value,tag_content_lon1);

   //lattitude2
   char tag_lat2[70]="lattitude2";
   char tag_content_lat2[100];
   fetch_tag(content,tag_lat2,tag_content_lat2);
   strcpy(pairset[7].tag,tag_lat2);
   strcpy(pairset[7].value,tag_content_lat2);


   //longitude2
   char tag_lon2[70]="longitude2";
   char tag_content_lon2[100];
   fetch_tag(content,tag_lon2,tag_content_lon2);
   strcpy(pairset[8].tag,tag_lon2);
   strcpy(pairset[8].value,tag_content_lon2);


   printf("\nproblem at 9982\n");
   if(type==0){
      // changing fetch_required

      char tag_freq_done[70]="frequencies_done";
      char tag_content_freq_done[100];
      fetch_tag(content,tag_freq_done,tag_content_freq_done);
      strcpy(pairset[9].tag,tag_freq_done);
      strcpy(pairset[9].value,tag_content_freq_done);

printf("\nproblem at 9992\n");
      char tag_fetch[70]="fetch_required";
      //char tag_content_fetch[100];
      //fetch_tag(content,tag_fetch,tag_content_fetch);
      strcpy(pairset[10].tag,tag_fetch);
      strcpy(pairset[10].value,value);

      char tag_prev_log_hash[70]="previous_log_hash";
      char tag_content_prev_log_hash[100];
      fetch_tag(content,tag_prev_log_hash,tag_content_prev_log_hash);
      strcpy(pairset[11].tag,tag_prev_log_hash);
      strcpy(pairset[11].value,tag_content_prev_log_hash);


   }else if(type==1){
      //changing freq_done

      char tag_freq_done[70]="frequencies_done";
      char tag_content_freq_done[100];
      fetch_tag(content,tag_freq_done,tag_content_freq_done);
      int aux;
      char aux1[20];
      aux=strtol(tag_content_freq_done,NULL,10);
      aux=aux+1;
      sprintf(aux1,"%d",aux);
      printf("\n\n%s %s %d\n\n",tag_content_freq_done,aux1,aux);
      strcpy(pairset[9].tag,tag_freq_done);
      strcpy(pairset[9].value,aux1);

      char tag_fetch[70]="fetch_required";
      char tag_content_fetch[100];
      fetch_tag(content,tag_fetch,tag_content_fetch);
      strcpy(pairset[10].tag,tag_fetch);
      strcpy(pairset[10].value,tag_content_fetch);

      char tag_prev_log_hash[70]="previous_log_hash";
      char tag_content_prev_log_hash[100];
      fetch_tag(content,tag_prev_log_hash,tag_content_prev_log_hash);
      strcpy(pairset[11].tag,tag_prev_log_hash);
      strcpy(pairset[11].value,tag_content_prev_log_hash);

   }else{
      // changing previous_log_hash
      char tag_freq_done[70]="frequencies_done";
      char tag_content_freq_done[100];
      fetch_tag(content,tag_freq_done,tag_content_freq_done);
      strcpy(pairset[9].tag,tag_freq_done);
      strcpy(pairset[9].value,tag_content_freq_done);


      char tag_fetch[70]="fetch_required";
      char tag_content_fetch[100];
      fetch_tag(content,tag_fetch,tag_content_fetch);
      strcpy(pairset[10].tag,tag_fetch);
      strcpy(pairset[10].value,tag_content_fetch);

      char tag_prev_log_hash[70]="previous_log_hash";
      //char tag_content_fetch[100];
      //fetch_tag(content,tag_fetch,tag_content_fetch);
      strcpy(pairset[11].tag,tag_prev_log_hash);
      strcpy(pairset[11].value,value);
   }


   char fileName[100]="/fs/microsd/log/recentPA.txt";
   key RFM_private_key;

   get_RFM_Key(&RFM_private_key);

  //  free(value_modulus);

   pair_file_write(pairset,12,fileName,RFM_private_key);
}


// functions remaining to define.

void Bundling_begins(){
   // files are already present, this function will generate a combined hash
   // file name:: (pa_Id_first_term)_log_(<=frequnecy_done)
   // generate bundledLog.txt containing:1)combine hash 2)signed combined hash 3)sign
   FILE *fptr;
   signed char ch;
   int i=0;
   char *content=(char*) malloc(4000*sizeof(char));
   fptr=fopen("/fs/microsd/log/log_of_logs.txt","r");

   while((ch=fgetc(fptr))!=EOF){
      content[i]=ch;
      i++;
   }

   content[i]='\0';
   fclose(fptr);

   // preparing for bundled.txt creation
   pair_set *new_content=(pair_set*) malloc(20*sizeof(pair_set));
   int new_cont=0;


   char tag_pa_id[50]="PA_ID";
   char tag_pa_id_contain[70];
   fetch_tag(content,tag_pa_id,tag_pa_id_contain);
   strcpy(new_content[new_cont].tag,tag_pa_id);
   strcpy(new_content[new_cont].value,tag_pa_id_contain);
   new_cont++;


   char tag_num_logs[70]="number_of_logs";
   char tag_logs_contain[40];
   fetch_tag(content,tag_num_logs,tag_logs_contain);
   strcpy(new_content[new_cont].tag,tag_num_logs);
   strcpy(new_content[new_cont].value,tag_logs_contain);
   new_cont++;



   char file_ini[30];

   i=0;
   while(tag_pa_id_contain[i]!='-'){
      file_ini[i]=tag_pa_id_contain[i];
      i++;
   }

   file_ini[i]='\0';

   int log_count=strtol(tag_logs_contain,NULL,10);
   char aux0[100];
   strcpy(aux0,"/fs/microsd/log/");
   strcat(aux0,file_ini);
   strcat(aux0,"_log_");
   char aux1[100];
   char aux_count[50];

   char *aux_hold=(char*) malloc(150*sizeof(char));

   char *combined_hash=(char*) malloc(1500*sizeof(char));


   //fetching hashes of each log and concating it into aux_hold.
   for (int u=1;u<=log_count;u++){

      strcpy(aux1,aux0);
      sprintf(aux_count,"%d",u);
      strcat(aux1,aux_count);
      strcat(aux1,".json");

      //fetching the hash
      fetch_tag(content,aux1,aux_hold);

      strcpy(new_content[new_cont].tag,aux1);
      strcpy(new_content[new_cont].value,aux_hold);
      new_cont++;
      if(u==1){
         strcpy(combined_hash,aux_hold);
      }else{
         strcat(combined_hash,aux_hold);
      }



   }
   free(aux_hold);


   free(content);

   strcat(combined_hash,"\0");

   printf("combined_hash ::: %s",combined_hash);

   key RFM_private_key;

   get_RFM_Key(&RFM_private_key);

   char hash_combined_hash[180];
   char signed_combined_hash[1080];

   signing_support_0(RFM_private_key,combined_hash,signed_combined_hash,hash_combined_hash);


   free(combined_hash);


   strcpy(new_content[new_cont].tag,"combined_hash");
   strcpy(new_content[new_cont].value,hash_combined_hash);
   new_cont++;


   strcpy(new_content[new_cont].tag,"signed_combined_hash");
   strcpy(new_content[new_cont].value,signed_combined_hash);
   new_cont++;


   char filename[30]="/fs/microsd/log/bundled.txt";
   pair_file_write(new_content,new_cont,filename,RFM_private_key);
   free(new_content);
}


int read_for_fetch(){
   char tag_fetch[30]="fetch_required";
   char tag_fetch_contain[20];
   char fname[40]="/fs/microsd/log/recentPA.txt";
   int status = call_file_read(fname,tag_fetch,tag_fetch_contain,0);
   if(status==1){
      // file is genuine
      int aux = strtol(tag_fetch_contain,NULL,10);
      if(aux==1){
         // fetching is required
         return 1;
      }else{
         // fetching is not required
         return 2;
      }

   }else{
      // file is not genuine
      return 3;
   }

}


/*
int check_fetch(){
   char fname[40]="./log/fetched.txt";

   FILE *fptr;
   fptr=fopen(fname,"r");
   if(fptr!=NULL){
      fclose(fptr);
      // fetched.txt is present

      key MC_key_public;// for validating the source to be no other than management client

      get_MC_public(&MC_key_public);


      int valid_status=Validating_File(fname,MC_key_public);
      // match the pa_id
      //if matches then delete the recentPA.txt, log of logs and log files (leaving space)

   }else{
      // fetched.txt is not present
   }

}
*/

void update_log_of_logs(char *tag,char *value,char *done_freq,char *pa_id,char *start_file){
   // this file will keep track of the log file names with the their hash value
   // then when bundling_begins function is evoked the file is used to form a bundled.txt
   // file with signed combined hash
   FILE *fptr;
   fptr=fopen("/fs/microsd/log/log_of_logs.txt","r");
   if(fptr!=NULL){

     signed  char ch;
      int a=0;
      char content_file[5000];
      while((ch=fgetc(fptr))!=EOF){
         content_file[a]=ch;
         a++;
      }
      content_file[a]='\0';
      fclose(fptr);
      // this is not the first time log is generated
      //first fetch the information from the file and then regenerate the file with additional tag and value.
      // note down the done_freq.
      //  done_freq-1 log file names would be there in the log_of_logs.txt.

      int count=atoi(done_freq);
      printf("\n%d\n",count);
 printf("\nproblem at 10126\n");
      pair_set *content=(pair_set*) malloc((count+4)*sizeof(pair_set));
      int content_count=0;
 printf("\nproblem at 10129\n");
      strcpy(content[content_count].tag,"PA_ID");
      strcpy(content[content_count].value,pa_id);
      content_count++;
 printf("\nproblem at 10133\n");
      strcpy(content[content_count].tag,"number_of_logs");
      strcpy(content[content_count].value,done_freq);
      content_count++;
 printf("\nproblem at 10137\n");
      char fetching_tg[100];
      strcpy(fetching_tg,start_file);
      char aux_st[100];
      char fetch_tag_contain[100];
      char aux0[100];
      printf("\nproblem at 10143\n");
      for(int i=1;i<count;i++){
         printf("\ni value   ;;;;;;;  %d\n",i);
         strcpy(aux_st,fetching_tg);
         sprintf(aux0,"%d",i);
         strcat(aux_st,aux0);
         strcat(aux_st,".json");

         fetch_tag(content_file,aux_st,fetch_tag_contain);


         strcpy(content[content_count].tag,aux_st);
         strcpy(content[content_count].value,fetch_tag_contain);
         content_count++;

         memset(aux_st,0,sizeof(aux_st));
         memset(fetch_tag_contain,0,sizeof(fetch_tag_contain));



      }

      strcpy(content[content_count].tag,tag);
      strcpy(content[content_count].value,value);
      content_count++;
      printf("\nproblem at 10165\n");
      char wr_filename[70]="/fs/microsd/log/log_of_logs.txt";

      key RFM_private_key;

      get_RFM_Key(&RFM_private_key);
       printf("\nproblem at 10171\n");

      pair_file_write(content,content_count,wr_filename,RFM_private_key);
       printf("\nproblem at 10174\n");
      free(content);


   }else{
      // this is the first time log is generated

      //char buff1[5000];
      pair_set content_holder[3];

      strcpy(content_holder[0].tag,"PA_ID");
      strcpy(content_holder[0].value,pa_id);

      strcpy(content_holder[1].tag,"number_of_logs");
      strcpy(content_holder[1].value,"1");


      strcpy(content_holder[2].tag,tag);
      strcpy(content_holder[2].value,value);


      char filename[70]="/fs/microsd/log/log_of_logs.txt";

      key private_rfm_key;

      get_RFM_Key(&private_rfm_key);

      pair_file_write(content_holder,3,filename,private_rfm_key);

   }






}


void fetching_publish_padata(){

   FILE *fptr;
   fptr=fopen("/fs/microsd/log/recentPA.txt","r");
   if(fptr!=NULL){
      signed char ch;
      int i=0;
      char content[5000];
      while((ch=fgetc(fptr))!=EOF){

         content[i]=ch;
         i++;

      }
      content[i]='\0';
     fclose(fptr);

   //publishing Permission Artefact information to pa_data.msg
   struct pa_data_s pa_data_raw;
   memset(&pa_data_raw, 0, sizeof(pa_data_raw));
   orb_advert_t pa_data_raw_pub = orb_advertise(ORB_ID(pa_data), &pa_data_raw);
      char lat1_tag[40]="lattitude1";
      char lattitude1[40];
      char lat2_tag[40]="lattitude2";
      char lattitude2[40];
      char long1_tag[40]="longitude1";
      char longitude1[40];
      char long2_tag[40]="longitude2";
      char longitude2[40];
      char st_tag[40]="Start_time";
      char start_time[80];
      char et_tag[40]="End_time";
      char end_time[80];
      char max_alt_tag[80]="maxAltitude";
      char maxalt[70];

      //fetching lattitude and longitude
      fetch_tag(content,lat1_tag,lattitude1);
      pa_data_raw.lattitude1=atof(lattitude1);

      fetch_tag(content,long1_tag,longitude1);
      pa_data_raw.longitude1=atof(longitude1);

      fetch_tag(content,lat2_tag,lattitude2);
      pa_data_raw.lattitude2=atof(lattitude2);



      fetch_tag(content,long2_tag,longitude2);
      pa_data_raw.longitude2=atof(longitude2);

      fetch_tag(content,max_alt_tag,maxalt);
      pa_data_raw.allowable_height=strtol(maxalt,NULL,10);

      //fetching start time and end time.
      fetch_tag(content,st_tag,start_time);
      fetch_tag(content,et_tag,end_time);

      Date_time startTime;
      Date_time endTime;
      conversion_to_dateTime(&startTime,start_time);
      conversion_to_dateTime(&endTime,end_time);

      pa_data_raw.start_hours=startTime.Hours;
      pa_data_raw.start_date=startTime.date;
      pa_data_raw.start_minutes=startTime.Minutes;
      pa_data_raw.start_month=startTime.Month;
      pa_data_raw.start_year=startTime.Year;
      pa_data_raw.start_seconds=startTime.Seconds;


      pa_data_raw.end_hours=endTime.Hours;
      pa_data_raw.end_date=endTime.date;
      pa_data_raw.end_minutes=endTime.Minutes;
      pa_data_raw.end_month=endTime.Month;
      pa_data_raw.end_year=endTime.Year;
      pa_data_raw.end_seconds=endTime.Seconds;

      int vehi_gps_pos=orb_subscribe(ORB_ID(vehicle_gps_position));

      struct vehicle_gps_position_s raw;
      orb_copy(ORB_ID(vehicle_gps_position),vehi_gps_pos,&raw);

      pa_data_raw.home_altitude=pow(10,-3)*(raw.alt);// in meters

      pa_data_raw.updated=1;


      orb_publish(ORB_ID(pa_data), pa_data_raw_pub, &pa_data_raw);



   }else{
      printf("\nfile is not opening\n");
   }

}




/////////////////////  logging json content


double max(double a,double b){
    if(a<b){
        return b;
    }
    else{
        return a;
    }

}


int fetch_previous_log_hash(char *result){
    // this function will fetch hash value from recentPA.txt
    FILE *fptr;
    fptr=fopen("/fs/microsd/log/recentPA.txt","r");
    signed char ch;
    int i=0;
    char *content=(char*)malloc(sizeof(char)*5000);
    while((ch=fgetc(fptr))!=EOF){
        content[i]=ch;
        i++;
    }
    content[i]='\0';
     fclose(fptr);
    char pr_log_hash_tag[50]="previous_log_hash";
    char pr_log_hash_contain[100];

    fetch_tag(content,pr_log_hash_tag,pr_log_hash_contain);
    free(content);
    if(strcmp( pr_log_hash_contain,"None")==0){
        // this is the first log for a particular PA.
        return 0;

    }else{
        strcpy(result,pr_log_hash_contain);
        return 1;
    }


}



void signing_log_content(char *HEX,char *cipher_string_hex){
      printf("\nproblem at line 668 possible mp iniialization %s\n",HEX);

    // modulus in public key
    mp_int modulus;
    mp_init(&modulus);
    char Modulus[513] ="a7f5781267c153741d140afab399e1d1052bcfcc075f146c09ac7ba0eef6912c1a8e89f4c8c2440d8123f1b6f8b0545650260135a65d93fda39d540935e7462bae44124e6dd2a1eb62c570b4f276d863870ba9e5d3feb6e1923a5bda8ef28b1355cada3e710f8c0227035848303c0966a577d567d1b4dec1651f244cfe7b06a5e853b83517d00939bfa6450c3dec551bf0a673d6fac7563b737a4c80f11cf6c57aeb5997e951a1521d36d881e4ca0e9b83726c78ac6f6b5eedc6e1ffc9a247e20fb9ac306fa92e8749f689d1353a32fd9f5eeb2c71a839975c8c7e0fdc25bbda1c678d6fa71962be337854cd8e6d0043ab4ce1e5de035798451d73ee25c487f3";
   // char Modulus[513] ="c8f770862ad1bc391310984fde3bb0879590e44e25dfbf9177431b27af55bb21bb6e93358d22d6ea60a38ce9e08cb30077e3d8e6908a2cc46bde10e767c0cc8977d66322be0fc196e5f399675a22a666563ce7de7af613c9dc23f79df8abd8e99391d957de97e9950bd2e1a8671404e8a9d19c32cae9dd5f6d29474cf45b613a31596bead1db45c7a929994fc99bdf68426ab3e9857b61fee826c01ee81c9fc030754252e26822ee26d5931b8835d221cc7889cdefec20ccc4312934d52b4de3cdfcb258dfefbf73ad365fb04a18edd3b8d3f68196fb910b791f6dbcd60b5c25b7f646176f9b009305c47d1b6c2cfe1ebddb58c08ad2da9dc932851ba9e258ef";
	//char Modulus[513] = "bc07d529450214ef63a8d61966987e8ca0594d9a7ec4f1881117b4f8ecbdc74b8769f6c98bfe931c9474116be8bd36527acfd95f6633d12cc8a960ab3d3e7a0b4b3e4990b594ee61af3b56315337501225525fb997b65c38118d614601dcb8bd631673a510498f2c3dab44d723d8b6daa697d0108e7fcb4d27525f386e7fcd9ce29c4ab12c4258aa77872259a25804791a1eaef54b65226ec84765442ac839db30467d86910e700d802807de1f4fef5235738d66359cb0a2707cb9cd90e90bb1f2d0d807aafbd048b1ddbb156d4984cfbbaa9a435b9230d213140dd5be64b5e594945d474665eaf5267fc598a5f75b99f83b029971b80c4149891d43abe62b95";
	//char Modulus[513]="cfb7e3ed5fb094ca81a2da9b07db403fc3882a33acf45f34c57d8fe677b3c787e3f034b4cfeba70f8a7beb95a53520780a59b1494ce00393aa812edf96abdafab551dd6728cd5d4ac4371c6049acecc6b63972737e60773129fa86dd8c5277ed64d13febb808776f630f7afb895789f609b28ba37caf0b14381f48ab123c093dc2621c34c570144d912de27fbf1d5669d427a07c9cf027013c93167d7faab8e8b4c95cfd0a9be0e4a92ae132713891a63af36315fce8d6f376cd2be410ef2155be9182caba4772275c93b47d5b8eed6cfe3fda8be884f79d49de09d02e42d93f350a25486a90653b0c72b6fb4d7bf9d40a2bc1559785b4bd0dd24a70aca7ef93";
    mp_read_radix(&modulus,Modulus,16);

        //private key exponent
    mp_int Private_key;
    mp_init(&Private_key);
    char private_hex[513]="011a1d3591b4b50578035fa711729b06b20ffd870b2d5686f6f148c65f8b029cd577c5f3f33543190b95deca228b95a213588c7d7b9ff58e9e7a33c8f3af96c846966fc04ffc27cfd190161339dc09c36d69682df7dc1dfb10e88d1dbbfe5f673b12dfa7b53a32e2f8ba8ee3ba5d4a7a7fef6f59050938c4f012ea3c0f17638dc459cde1898e47cabcb1d3c92e9ea09775738a84a4cf951cbe697a6ae9f332bd0351a5b7be91a8cb939f7d5b76b2753fbb903721cd181f705f724303bfdc56927fbb44aa19ed9a76b4c5165e356d5cf06fc543803e55431757d0355a48e90c6b73c5210f4f886964a2ae8432af68223bf40fe1a0d68a32b1b80a0595747ae801";
    //char private_hex[513]="184885e9405d4d801c04a252ec488c2125fa770bd659bdfd26cb0e09f28eca68de0c136fa219369ce5867dad78fba75984231cff6731bb0d14f7a55540dd3419dc48247c7b38ce2c9ca69dbfb64d7f8bd819cdeebd2ee4df3c6180372f681c72c4e917b91d657fcd09bbb696b1b5e28df68f246fa2c33583a55e1a867af45bc004583335a5d5a79d046b2d9c74bb5b7282a4ffb0167c3f0c982051fe6e253a8c20d063e1eccf0dfcdcaf4c34a3f0ca2cf8ebbb37d03c8176ed157b55f738c5f47fe4d0d42b05c6c0e9a58cc3766ec87926ee8e38843b25171886ec413c644b931acbb92895c4665274a17355d1d72401339d59ee0ec69a65d613e9b5c215a7c9";
    // char private_hex[513] = "9b78ea83264133684182400d5eaca6aec68330cc97176712f7f71f3758210f61df44f9beead78372753987922f2e0c75a480aa1edc95e9d65ad0da529ce044ef83b6ac03507125ae75c2dd61098ac9d54730d65fd21702278633dd8392549c18548f22ee100a92aca50d316da68131a897691dac22f77df57c96fa8ee1a7212db313396410a5c9c8a31f6f940724cff2b2db5eb078eedad92b6ff29a8636fcd370e99773e96168f34839693f84b7a083597bfbe0f674c79b2348b038ca730cada30bcf2dd9cde27dd555891d3cc10b7831b23e7cda163570635727f11d569492a201f55c56d9a92d46f71b6ecea30f28f8c040f834a2da43f72a1ec927df9441";
    //char private_hex[513]="6beefbcaae7c4cf46524403f6a77ad0cf5075e1677fa8b361aa0c2135983db5c6b3eb7c4747dd8d3247c7bcfc886b0966f9a679ad50d5a0e72fca964992037ab2a689d892b147b338c7dae8b01fd8f133a40e38dcbcf48600d96165a2cbdf57f2f71e3ab1277a3c8074b55f63a497870965d665dcf3e0d9db603db78b902e5317897c3538a24057c29f25f98284c7262bd965f4923da6a117f37c88e78ec6ce34da38eae946c06cc4ac166df70a7dd24b1be2254c8a25797bde3c58918b9f8ae0b5bf8019f8f60fb9e7e087d20bcb935b0b9c3a43c0b7f64ca9027c8bfa9e83525781dfd7b529a88a305ae8b4ec56e2cadff1ded53f2ec461ae2bd56e26231f1";
    mp_read_radix(&Private_key,private_hex,16);//
    //char yu[]="b59cab400e9c64525f566e85c97c2d65a2d0620f0479cb8ef8e6eb1d3c1cf093be5eb827a4d01d9a585b0043f148d9e3c213676d819818be05bf39a7657cd953bd4a7b2e361d03e1d64d99c755036eb1f97223832b1d596ee131f0adb84879cbd6763f1fbe61c214ff6514551900ef6034feb20cf0daf995ebaf9d7e5c3566cf012e1aa32e3733747be37955ec0d127d32fbddfdd8451fc5254dfe09a4b1b26bd9ec3487cd04121fd4b965456c024e21086da38c81670b8cb57fdad93f8f90eb6bfa109d57ffda156c105fe15c4c78d588de116ff06abc1ea03546ef36b1ac711cfe4feb64250264ac5ff2704a76a6b1fe04480af906e9ee71836190c345bf67";
        /*const char *content_to_hash =(const char*) malloc(sizeof(char)*4000);

        int rsa_encryption_count_aux=0;
        while(HEX_format_Digest[rsa_encryption_count_aux]!='\0'){
            content_to_hash[rsa_encryption_count_aux]=HEX_format_Digest[rsa_encryption_count_aux];
            rsa_encryption_count_aux++;
        }*/
    //making the signing process pkcs.1.15 compatible
    char padding_SHA256[810]="1ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff003031300d060960864801650304020105000420";
    strcat(padding_SHA256,HEX);
    printf("\npadded sha digest :%s\n",padding_SHA256);
    mp_int hash_to_sign;
    mp_init(&hash_to_sign);

    mp_read_radix(&hash_to_sign,padding_SHA256,16);
    char *aux_hash_ex=(char*) malloc(1014*sizeof(char));
    memset(aux_hash_ex,0,1014);
    // char aux_hash_ex[1014];
    mp_to_hex(&hash_to_sign,aux_hash_ex,1014);
    printf("\n   hash to sign ===\n%s\n\n",aux_hash_ex);
    free(aux_hash_ex);

    //Signing begins

    mp_int cipher;
    mp_init(&cipher);

    mp_exptmod(&hash_to_sign,&Private_key,&modulus,&cipher);
    //
    //char cipher_string[1013];
    char *cipher_string=(char*) malloc(1013*sizeof(char));
    mp_to_hex(&cipher,cipher_string,1013);
    printf("\ncipher string in hex %s\n",cipher_string);//this is the encrypted message
        strcpy(cipher_string_hex,cipher_string);
        free(cipher_string);
        mp_clear(&modulus);
        mp_clear(&Private_key);
        mp_clear(&hash_to_sign);
        mp_clear(&cipher);


}


void get_PA_ID(char *res){
    FILE *fptr;
    fptr=fopen("/fs/microsd/log/recentPA.txt","r");
    if(fptr==NULL){
        printf("\nfile is not opening\n");
    }
    //char content[5000];
    char *content=(char*)malloc(sizeof(char)*5000);
    signed char ch;
    int i=0;

    while((ch=fgetc(fptr))!=EOF){
        content[i]=ch;
        i++;
    }
    content[i]='\0';
    fclose(fptr);
    char tag[50]="permissionArtifactId";

    fetch_tag(content,tag,res);
    free(content);


}

void writingKey_Value(char *str_ptr,char *Key,int *count,char *value,int last){

    int i=*count;

    str_ptr[i]='"';
    i++;

    for(int y=0;y<int(strlen(Key));y++){
        str_ptr[i]=*(Key+y);
        i++;
    }
    str_ptr[i]='"';
    i++;
    str_ptr[i]=':';
    i++;
    str_ptr[i]=' ';
    i++;


    //this is for string

    str_ptr[i]='"';
    i++;
    for(int y=0;y<int(strlen((char*)value));y++){
        str_ptr[i]=*((char*)value+y);
        i++;
    }
    str_ptr[i]='"';
    i++;


    if(last==1){
        str_ptr[i]='\n';
        i++;
    }else{
        str_ptr[i]=',';
        i++;
         str_ptr[i]='\n';
        i++;
    }
    *count=i;

}


void writingKey_Value(char *str_ptr,char *Key,int *count,int value,int last){

    int i=*count;

    str_ptr[i]='"';
    i++;

    for(int y=0;y<int(strlen(Key));y++){
        str_ptr[i]=*(Key+y);
        i++;
    }
    str_ptr[i]='"';
    i++;
    str_ptr[i]=':';
    i++;
    str_ptr[i]=' ';
    i++;

    //this is for integer
    char str[30];
    sprintf(str, "%d", value);

    //putting char from str to str_ptr
    int i_count=0;
    while(str[i_count]!='\0'){
    str_ptr[i]=str[i_count];
    i_count++;
    i++;
    }



    if(last==1){
        str_ptr[i]='\n';
        i++;
    }else{
        str_ptr[i]=',';
        i++;
         str_ptr[i]='\n';
        i++;
    }
    *count=i;

}


void writingKey_Value(char *str_ptr,char *Key,int *count,long value,int last){

    int i=*count;

    str_ptr[i]='"';
    i++;

    for(int y=0;y<int(strlen(Key));y++){
        str_ptr[i]=*(Key+y);
        i++;
    }
    str_ptr[i]='"';
    i++;
    str_ptr[i]=':';
    i++;
    str_ptr[i]=' ';
    i++;

    //this is for integer
    char str[30];
    sprintf(str, "%ld", value);
    printf("\n\n\n string unsigned int  ::::::::::  %s\n\n\n",str);

    //putting char from str to str_ptr
    int i_count=0;
    while(str[i_count]!='\0'){
    str_ptr[i]=str[i_count];
    i_count++;
    i++;
    }



    if(last==1){
        str_ptr[i]='\n';
        i++;
    }else{
        str_ptr[i]=',';
        i++;
         str_ptr[i]='\n';
        i++;
    }
    *count=i;

}


void writingKey_Value(char *str_ptr,char *Key,int *count,double value,int last){

    int i=*count;

    str_ptr[i]='"';
    i++;

    for(int y=0;y<int(strlen(Key));y++){
        str_ptr[i]=*(Key+y);
        i++;
    }
    str_ptr[i]='"';
    i++;
    str_ptr[i]=':';
    i++;
    str_ptr[i]=' ';
    i++;


    // this is for float
    char str[30];
    sprintf(str, "%f", value);

    //putting char from str to str_ptr
    int i_count=0;
    while(str[i_count]!='\0'){
    str_ptr[i]=str[i_count];
    i_count++;
    i++;
    }




    if(last==1){
        str_ptr[i]='\n';
        i++;
    }else{
        str_ptr[i]=',';
        i++;
         str_ptr[i]='\n';
        i++;
    }
    *count=i;

}


void writingKey(char *str_ptr,char *Key,int *count){

    int i=*count;

    str_ptr[i]='"';
    i++;

    for(int y=0;y<int(strlen(Key));y++){
        str_ptr[i]=*(Key+y);
        i++;
    }
    str_ptr[i]='"';
    i++;
    str_ptr[i]=':';
    i++;
    str_ptr[i]=' ';
    i++;
    *count=i;

}



void putSpace(int count_Space,char *ptr_string ,int *i_ptr){
    int y=*i_ptr;
    for(int u=0;u<count_Space;u++){
        ptr_string[y]=' ';
        y++;
    }
    *i_ptr=y;
}

void objectCreator(char *ptr,Geo_tag geo,int last,int space_count,int *count){

    /// this function will make object entries inside LogEntries
    /*
    {
        "Entry_type": "TAKEOFF/ARM",
        "TimeStamp": 41,
        "Longitude": 380.000000,
        "Latitude": 443.5000000,
        "Altitude": 704.3400,
    }
    */


    char Geolat[]="Latitude";// float  Latitude in Degrees East
    char Geolon[]="Longitude";//float  Longitude in Degrees North
    char GeoTime[]="TimeStamp";//milliseconds, type : integer
    char  GeoAlt[]="Altitude";//Ellipsoidal Height in Meters  type: integer
    char GeoEntryType[]="Entry_type";//entry_type

    int iterator=*count;

    ptr[iterator]='{';
    iterator++;
    ptr[iterator]='\n';
    iterator++;
    space_count=space_count+2;
    // putting spaces
    putSpace(space_count,ptr,&iterator);
      /// Lattitude
    writingKey_Value(ptr,Geolat,&iterator,geo.Lattitude,0);
    putSpace(space_count,ptr,&iterator);

    // starting with entrytype
    char entry[20];
    if(geo.Entrytype==0){
         strcpy(entry,"GEOFENCE_BREACH");
    }else if(geo.Entrytype==1){
         strcpy(entry,"TAKEOFF/ARM");
    }else if(geo.Entrytype==2){
         strcpy(entry,"TIME_BREACH");
    }else if(geo.Entrytype==3){
         strcpy(entry,"LAND");
    }else{
        strcpy(entry,"CRASH");
    }
    writingKey_Value(ptr,GeoEntryType,&iterator,entry,0);
    putSpace(space_count,ptr,&iterator);

    /// Altitude
    writingKey_Value(ptr,GeoAlt,&iterator,geo.Altitude,0);

    putSpace(space_count,ptr,&iterator);
    /// Longitude
    writingKey_Value(ptr,Geolon,&iterator,geo.Longitude,0);
    putSpace(space_count,ptr,&iterator);

     /// TimeStamp
    writingKey_Value(ptr,GeoTime,&iterator,geo.Timestamp,1);
    space_count=space_count-2;
    putSpace(space_count,ptr,&iterator);


    //closing

    ptr[iterator]='}';
    iterator++;

    *count=iterator;


}



void  log_naming_support(char *paID_firstTerm,char *done_freq){

    FILE *fptr;

    fptr=fopen("/fs/microsd/log/recentPA.txt","r");
    signed char ch;
    int i=0;
    char *content=(char*)malloc(sizeof(char)*5000);
   // char content[5000];

    while((ch=fgetc(fptr))!=EOF){
        content[i]=ch;
        i++;
    }
    content[i]='\0';
    fclose(fptr);
    char tag_paID[50]="permissionArtifactId";
    char tag_paID_contain[80];
    char tag_done_freq[50]="frequencies_done";
    char tag_done_freq_contain[100];

    fetch_tag(content,tag_paID,tag_paID_contain);
    fetch_tag(content,tag_done_freq,tag_done_freq_contain);
    free(content);
    i=0;
    while(tag_paID_contain[i]!='-'){
        paID_firstTerm[i]=tag_paID_contain[i];
        i++;
    }
    paID_firstTerm[i]='\0';

    strcpy(done_freq,tag_done_freq_contain);

}

void main_json_file_writing(Geo_tag *data_array, int length_dat){
    char *buff=(char*) malloc(sizeof(char)*10000);
    //char buff[10000];// this is where json would be written, later copied to a file

    int space_count=0;
    //object {}
    char object_start='{';
  //  char object_end='}';

    //array []
    char array_start='[';
    //char array_end=']';
    // ""
   // char quote_start='"';
   // char quote_end='"';
    //colon :
  //  char colon=':';
    //comma
   // char comma=',';


    // writting JSON flight log
    int i=0;

    buff[i]=object_start;
    space_count=space_count+2;
    i++;
    buff[i]='\n';
    i++;
    // after { and new line spaces need to be put

    putSpace(space_count,&buff[0],&i);

    // starting with writing Flight log
    int write_signature_help=i;
    char firstKey[]="FlightLog";

    writingKey(&buff[0],firstKey,&i);

    int str_start=i;

    buff[i]=object_start;
    i++;
    buff[i]='\n';
    i++;
    space_count=space_count+2;

    // after { and new line, spaces need to be put

    putSpace(space_count,&buff[0],&i);



    ////// LogEntries////
    char Log[]="LogEntries";

    writingKey(&buff[0],Log,&i);

    buff[i]=array_start;
    i++;
    buff[i]='\n';
    i++;
    space_count=space_count+2;

    // after { and new line, spaces need to be put

    putSpace(space_count,&buff[0],&i);
    // now a loop will run
    for (int counter=0;counter<length_dat;counter++){
        objectCreator(&buff[0],data_array[counter],0,space_count,&i);
        if(counter==length_dat-1)
        {
            buff[i]='\n';
            i++;
            space_count=space_count-2;
        }else{
            buff[i]=',';
            i++;
            buff[i]='\n';
            i++;
        }

        putSpace(space_count,&buff[0],&i);
    }
   /* space_count=space_count-2;
    putSpace(space_count,&buff[0],&i);*/
    buff[i]=']';
    i++;
    buff[i]=',';
    i++;
    buff[i]='\n';
    i++;
    putSpace(space_count,&buff[0],&i);


    printf("\nproblem at line 516\n");


    // writting permission artefact id from recentPA.txt


    char PA_id[70];//="ABDHUBDUBBEBDIBWBID";

    get_PA_ID(PA_id);
printf("\nproblem at 525\n");
    char PA[]="PermissionArtefact";
    writingKey_Value(&buff[0],PA,&i,PA_id,0);
    printf("\nproblem at 528\n");
    putSpace(space_count,&buff[0],&i);
    printf("\nproblem at 530\n");


    // fetching previous log hash from recentPA.txt
    // if None then its the first log for the particular PA
    // writing previous log hash
    char previous_log_hash[70];//="aaaanfbofburbfubvoinvknvofn";
    char previous_log_hash_2[70];//="aaaanfbofburbfubvoinvknvofn";
    printf("\nproblem at 536\n");
    int hash_status=fetch_previous_log_hash(previous_log_hash_2);

    if(hash_status==0){
        // no previous log hash is present.
        strcpy(previous_log_hash,"First_log");

    }else{
        strcpy(previous_log_hash,previous_log_hash_2);
    }

    char prev_hash[]="previous_log_hash";
    writingKey_Value(&buff[0],prev_hash,&i,previous_log_hash,1);



   // putSpace(space_count,&buff[0],&i);

  printf("\nproblem at line 554\n");
   space_count=space_count-2;
   putSpace(space_count,&buff[0],&i);
   int str_end=i;/////// the point where value of flightlog key ends
   buff[i]='}';
   i++;
   buff[i]='\n';
   i++;
   space_count=space_count-2;
   putSpace(space_count,&buff[0],&i);
   buff[i]='}';
   i++;
   buff[i]='\0';
    // At this point, json file is completed
    // next step is to canonicalize the value(object) under "FLightlog" key
   char *canonicalized_flight=(char*) malloc(sizeof(char)*14000);
  // char canonicalized_flight[14000];
   int canon_count=0;


   for(int h=str_start;h<=str_end;h++)
    {
    /* canonicalize hack : all the serialization under objects have been taken care of
    earlier only, now just keep on taking the elements one by one from one array to another
    except new lines and spaces
    */
        if (buff[h]!=' ' && buff[h]!='\n'){
            canonicalized_flight[canon_count]=buff[h];
            canon_count++;
        }

    }
    canonicalized_flight[canon_count]='\0';


    // at this point canonicalization has been done
    // our next step is to implement SHA256 algo

    signed char ch;
    unsigned char *st=new unsigned char[10000];
    unsigned int i_sha=0;
    //// assigning char to unsigned char, this has to be done to implement sha256
    while(1){
	   ch=canonicalized_flight[i_sha];
      // printf("%c",ch);
	   if(ch=='\0'){
		   break;
	   }
        *(st+i_sha)=ch;
		printf("%c",*(st+i_sha));
		i_sha++;

    }
  free(canonicalized_flight);
  printf("\nproblem at line 607\n");

   SHA256 sha;

	sha.update(st,i_sha);

	uint8_t * digest = sha.digest();
    int it=0;//just an iterator
	//printing digest
	printf("\n");

    while(it<32){
		printf("%02x ",*(digest+it));
		it++;
	}
    char HEX_format_Digest[65];

    // code for making hex string
    int hex_count=0;
    for(int yui=0;yui<32;yui++){
         //printf("%d",*(digest+yui)>>4||0x0f);
         HEX_format_Digest[hex_count]=HEX_ele[*(digest+yui)>>4 & 0x0f];
         hex_count++;
         HEX_format_Digest[hex_count]=HEX_ele[*(digest+yui) & 0x0f];
         hex_count++;
      }
   HEX_format_Digest[hex_count]='\0';

   printf("\nhere is the string:%s\n\n",HEX_format_Digest);

   // this is the hash of the current log file, this will be updated in recentPA.txt.
   //
   update_recentPA(2,HEX_format_Digest);

    //At this point we have got the string format for the digest


//	delete[] digest;

   delete []st;
   st=NULL;
  // printf("\n\n");
   // printf("%s",canonicalized_flight);
   //printf("\n\n");
    //printf("%s",buff);
  // printf("\n");

/*
   FILE *fptr;
   fptr=fopen("flight_log1.json","w");// this is the resultant json file
   fprintf(fptr, "%s", buff);
   fclose(fptr);
*/


    // Now RSA encryption is required
    /* following elements are required
    1)Private key
    2) above calculated digest HEX_format_digest
    */
    //////////////////////// RSA Encryption ///////////////
    char *cipher_string_hex=(char*) malloc(1013*sizeof(char));

    signing_log_content(HEX_format_Digest,cipher_string_hex);





    //char ty[]="734CB799A64A670D68F3AA84B2542BBEE1F7FD5AC022460CDD63E297DF55D8B04CBCB1103BDC720EB85D090CCA091D067A5F54BFADBECB09A2EBDDA9D92E00AF83D5E9FED29B38C295A315CBD52CB4ED1645BC3707D4ED165E5E4D906A1C149E5073AE6A8D5FC74CC11D63885CEFE236E4464C0D387A5861234814FF31A8A3A0B2D6022B8A8FAEC79DF61B638B0B85392F152D4743BFD2779C4472B73CF952027C83A8DB7EEA32D6A1B96ED8F92CF7576EBBB35F8CA341C69E60C352BD79628C17585030D8BF07C8E6D73FC0AB51D2226B1E7C42D5DA0BC4EAF5063C8F72601AC0090D92B74C17365463B3F07B0397692AF204F58A0ACB665CBA4C336121C4CA";
    /*
    mp_int Public_key;
    mp_read_radix(&Public_key,"65537",10);
    mp_int message;
    mp_exptmod(&cipher, &Public_key,&modulus,&message);                          this part for decrypting encrypted text

   char message_string_hex[513];
    mp_to_hex(&message,message_string_hex,sizeof(message_string_hex));
    printf("\n%s\n",message_string_hex);//this is the decrypted message
    */


   /////////////////////// base64 Encoded //////////////// hex to base64
   char *res=(char*) malloc(sizeof(char)*2500);
   //char res[2500];
   base64Encoder(cipher_string_hex,strlen(cipher_string_hex),res);
   free(cipher_string_hex);
   printf("\n\nbase 64 encoder %s\n",res);
   char signature[700]="\"signature\": \"";
   strcat(signature,res);
   free(res);
   char t[4];
   char *gio=&t[0];
   t[0]='"';
   t[1]=',';
   t[2]='\0';
   strcat(signature,gio);
   printf("\n%s\n",signature);
   printf("the length of string %d",(int)strlen(signature));


   /////////////////////// Writing Signature ///////////////// open file and rewrite
   /* string insertion
   1) length of string that has to be inserted  l1   res
   2) length of string  where we want to insert  l2  buff
   */
    int i_sign=0;
    char *buff1=(char*) malloc(sizeof(char)*15000);
    //char buff1[15000];
    space_count=0;
    buff1[i_sign]=object_start;
    space_count=space_count+4;
    i_sign++;
    buff1[i_sign]='\n';
    i_sign++;
    // after { and new line, spaces need to be put

    putSpace(space_count,&buff1[0],&i_sign);
    int you=0;
    while(signature[you]!='\0'){
      buff1[i_sign]=signature[you];
      i_sign++;
      you++;
   }
   buff1[i_sign]='\n';
   i_sign++;
   putSpace(4,&buff1[0],&i_sign);
   for(int iu=write_signature_help;iu<int(strlen(buff));iu++){
      buff1[i_sign]=buff[iu];
      i_sign++;

   }
   free(buff);
   buff1[i_sign]='\n';
   i_sign++;
   buff1[i_sign]='\0';
   printf("\n\n%s\n\n",buff1);

   //deciding file name for logs
   // current approach : take note of the done frequencies inside recentPA.txt
   //and last four letters of pa id first term.
   //filename ==== paidfirst_term-log-(done_frequency+1)

   char paID_firstTerm[20];
   char done_freq[10];

   log_naming_support(paID_firstTerm,done_freq);

   char file_directory[50]="/fs/microsd/log/";

   char filename[100];

   strcpy(filename,file_directory);
   strcat(filename,paID_firstTerm);
   strcat(filename,"_log_");
   char aux_start_filename[40];
   strcpy(aux_start_filename,filename);

   strcat(filename,done_freq);
   strcat(filename,".json");

   printf("%s\n\n",filename);

   update_log_of_logs(filename,HEX_format_Digest,done_freq,PA_id,aux_start_filename);

   FILE *fptr1;
   fptr1=fopen(filename,"w");// this is the resultant json file
   fprintf(fptr1, "%s", buff1);
   fclose(fptr1);
   free(buff1);

}



int check_Geobreach(pa_data_s data,vehicle_global_position_s vgp){
    // return 1 if geofence breached
    // return 0 if geofence not breached
    double lattitude1 =data.lattitude1;
    double longitude1 =data.longitude1;
    double lattitude2 =data.lattitude2;
    double longitude2 =data.longitude2;
    int home_altitude =data.home_altitude;

    double mini_lat=min(lattitude1,lattitude2);
    double maxi_lat=max(lattitude1,lattitude2);

    double mini_lon=min(longitude1,longitude2);
    double maxi_lon=max(longitude1,longitude2);
    double allowed_height=0.3048*(data.allowable_height);// converting feets to meters


    double main_lat=vgp.lat;
    double main_lon=vgp.lon;
    double main_alt=vgp.alt;

    if(main_lat>=mini_lat && main_lat<=maxi_lat){
        if(main_lon>=mini_lon && main_lon<=maxi_lon){
            if((main_alt-(double)home_altitude)>allowed_height){
                return 1;

            }else{
                return 0;
            }

        }else{
            return 1;
        }


    }else{
        return 1;
    }





}


int check_Timebreach(pa_data_s data,vehicle_gps_position_s vgp){
    // return 1 timebreach
    // return 0 no timebreach
    Date_time start;
    Date_time end;

    start.date=data.start_date;
    start.Hours=data.start_hours;
    start.Minutes=data.start_minutes;
    start.Month=data.start_month;
    start.Seconds=data.start_seconds;
    start.Year=data.start_year;

    end.date=data.end_date;
    end.Hours=data.end_hours;
    end.Minutes=data.end_minutes;
    end.Month=data.end_month;
    end.Seconds=data.end_seconds;
    end.Year=data.end_year;



    Date_time current;

    time_t timestamp;
    timestamp=(vgp.time_utc_usec)/1000000;//microsecond to seconds
	//printf("\ntimestamp is seconds: %u\n",timestamp);
    struct tm  ts;

    ts = *localtime(&timestamp);
    current.Year=ts.tm_year+1900;
    current.Month=ts.tm_mon+1;
    current.date=ts.tm_mday;
	current.Hours=ts.tm_hour;
	current.Minutes=ts.tm_min;
	current.Seconds=ts.tm_sec;

    int status=In_Time(current,start,end);

    if(status==0){
        return 1;

    }
    return 0;




}



void check_debug_mp(char *res){

         char Signature_Value_in_hex[556];//;=(char*) malloc(500*sizeof(char));
         char Modulus[554];
            //  mp_int ; dgca public key modulus
				strcpy(Modulus,"d6912a773335d3a193c742071762794e26bcac49d78b6b65a784e1c18d18f2f88b9f13d7fc41a5b9e11e3d75905b0f8373dbac5658962940d104711467814b7319774014644df82e9764b4e852cb7547d56de3224aada000db231b1356ffe86391b3eca2ddf67d44f40ea6e84092cc67387db3a2042487b8dacfe5b588738973a59f5fa813a33e4bb8fbafa407794db990a307201b0a0fbd92ddd181a868a6cfa64625353714e9b1de6f81f3addbf5b202030dbd9e51385db9314591f01af01f07cc92f18ebe60545cea53eb54438e251fd88e7380b5a0612c29ccb7dfd88f7a7d7f0dd07a9b596923602798ad9f1bbb89fffd3cdb8fb94c48640d27fa681e47");
				// char Modulus[513]="ab9d5c8d1fe67207749d63b7dcedd233ce32bb70d175a1bc38c612ab33e2c58e51f83f2788e4d52d9bceb5a1513929de3f526650071a067e6c161b05c60a495fc3ba79ed26f4fa8b2fe2ca8dec44b39759f39206f06a85f9424005a29f05e4cf3a0239340c28c993c1a61cf1b2b6b57c7d8e576ae86827f812b327625baec9ecbf55f1651d35600b9f955f6c2f3bea3aa5852ecdd36a0af818c19acc1030979bed3c89993faa92e0aa0502413b3ca86bbf63477f12ac069aff7137cb72c57f886da79033bbb3b4df0f6cc7fcc18e343aa76036681a566311e267c03b65c98abc91e58f090020c67f776199c0eb76d7e6363687475d3da36ff050f85275607fdd";
            // then processing is of PA
            char file_name[]="/fs/microsd/log/permission_artifact_breach.xml";
	         char Tag_Signed_Value0[17]="SignatureValue";
				//  char Signature_Value0[550];// in base 64
				char *Signature_Value0=(char*) malloc(550*sizeof(char));

				getTagvalue(Tag_Signed_Value0,Signature_Value0, file_name);

				base64decoder(Signature_Value0, Signature_Value_in_hex);

				free(Signature_Value0);
				//free( Signature_Value_in_hex);
				int cb;
				mp_int message,modulus,public_key,Decrypted;
				cb= mp_init_multi(&message,&modulus,&public_key,&Decrypted,NULL);
				cb= mp_read_radix(&message,Signature_Value_in_hex,16);
				// free(Signature_Value_in_hex);

				/// Public key of dgca: This has to be taken from a reserved file in directory./ firmware update

				cb= mp_read_radix(&modulus,Modulus,16);

				// mp_int ;
				cb= mp_read_radix(&public_key,"65537",10);

				//mp_int ;
				cb= mp_exptmod(&message, &public_key,&modulus,&Decrypted);   //                       this part for decrypting encrypted text
				printf("%d",cb);
            char message_string_hex[513];
            mp_to_hex(&Decrypted,message_string_hex,sizeof(message_string_hex));
            mp_clear_multi(&message,&modulus,&public_key,&Decrypted,NULL);

            strcpy(res,message_string_hex);
         //   strcat(res,"/0");

}


void valid_chunk_refe(FILE *fptr, char *re_sha, int type){

  // type 1== Reference section
  // type 2== Signed_info section

   char content_canonilized[2000];//=(char*) malloc(2000*sizeof(char));
   char xmlExeclusive[2000];//=(char*) malloc(2000*sizeof(char));
   char output[2000];//=(char*) malloc(2000*sizeof(char));
   char Sha_of_Reference[100];//=(char*) malloc(300*sizeof(char));
   char file_name[]="/fs/microsd/log/permission_artifact_breach.xml";

   if(type==1){

      Reference_canon(file_name,content_canonilized);
   }else{
       SignedInfo_canon(file_name,content_canonilized);
   }


   fprintf(fptr,"\n%s\n",content_canonilized);

   //free(Reference_canonilized);
   //fclose(fptr);
   cleanerXML(content_canonilized,output);


   fprintf(fptr,"\n%s\n",output);
   if(type==1){
      xmlExeclusiveCanon(output,xmlExeclusive);
      fprintf(fptr,"\nre %s\n",xmlExeclusive);
   }

   if(type==1){
      Sha256_implement(xmlExeclusive,Sha_of_Reference,NULL,0);
   }else{
      Sha256_implement(output,Sha_of_Reference,NULL,0);
   }


   fprintf(fptr,"\n%s\n",Sha_of_Reference);
   strcpy(re_sha,Sha_of_Reference);
}

/*
void valid_chunk_refe_2(FILE *fptr, char *re_sha){

   char s_canonilized[2000];//=(char*) malloc(2000*sizeof(char));
   char output[2000];//=(char*) malloc(2000*sizeof(char));
   char Sha_of_Reference[100];//=(char*) malloc(300*sizeof(char));
   char file_name[]="/fs/microsd/log/permission_artifact_breach.xml";
   SignedInfo_canon(file_name,s_canonilized);
   fprintf(fptr,"\n%s\n",s_canonilized);

   //free(Reference_canonilized);
   //fclose(fptr);
   cleanerXML(s_canonilized,output);
   fprintf(fptr,"\n%s\n",output);
   Sha256_implement(output,Sha_of_Reference);
   fprintf(fptr,"\n%s\n",Sha_of_Reference);
   strcpy(re_sha,Sha_of_Reference);

}*/


void pa_validation(){

      FILE *fptr;
      fptr=fopen("/fs/microsd/result.txt","w");

      char file_name[]="/fs/microsd/log/permission_artifact_breach.xml";
      char *Digest_Value0=(char*) malloc(100*sizeof(char));
      char Sha_of_Reference[70];
      char Sha_of_SignedInfo[70];
      valid_chunk_refe(fptr,Sha_of_Reference,1);
      fprintf(fptr,"\n\n%s\n",Sha_of_Reference);

      valid_chunk_refe(fptr,Sha_of_SignedInfo,2);
      fprintf(fptr,"\n\n%s\n",Sha_of_SignedInfo);

      char Tag_Digest_value0[12]="DigestValue";
      getTagvalue(Tag_Digest_value0,Digest_Value0,file_name);
      char *Digest_Value_in_hex=(char*) malloc(100*sizeof(char));
      base64decoder(Digest_Value0, Digest_Value_in_hex);
      fprintf(fptr,"\n\n Digest value %s\n",Digest_Value_in_hex);
      free(Digest_Value0);

      char *Extracted_sha_from_Signature_value=(char*) malloc(100*sizeof(char));
      char *message=(char*) malloc(750*sizeof(char));

      check_debug_mp(message);

      char padding_SHA256[447]="1ffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffffff003031300d060960864801650304020105000420";
      for(int o=0;o<444;o++){
            if(padding_SHA256[o]!=small_letter(message[o])){
               //return 2;
               fprintf(fptr,"%s","\n Bad SIGN \n");
               break;
            }
      }

      fprintf(fptr,"%s","\n SIGN is good \n");
      int k_cout=0;
      for(int i2=445;message[i2]!='\0';i2++){
            Extracted_sha_from_Signature_value[k_cout]=small_letter(message[i2]);
            k_cout++;
         // printf("%c",message_string_hex[i]);
      }
      Extracted_sha_from_Signature_value[k_cout]='\0';//
      fprintf(fptr,"\n%s\n",Extracted_sha_from_Signature_value);
      if(strcmp(Extracted_sha_from_Signature_value,Sha_of_SignedInfo)==0){

         if(strcmp(Sha_of_Reference,Digest_Value_in_hex)==0){
                  fprintf(fptr,"%s","\n 100000000000PA is valid and non tampered \n");

            }else{
                  fprintf(fptr,"%s","\n PA is not valid\n");

            }

         }else{
            fprintf(fptr,"%s","\n PA is not valid \n");
         //  PX4_INFO("PA is valid???????????");

         }
      free(Extracted_sha_from_Signature_value);
      free(message);
      free(Digest_Value_in_hex);
      fclose(fptr);

}
