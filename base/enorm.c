#include <math.h>

#include "minpack.h"

/*
 *     function enorm
 *
 *     given an n-vector x, this function calculates the
 *     euclidean norm of x.
 *
 *     the euclidean norm is computed by accumulating the sum of
 *     squares in three different sums. the sums of squares for the
 *     small and large components are scaled so that no overflows
 *     occur. non-destructive underflows are permitted. underflows
 *     and overflows do not occur in the computation of the unscaled
 *     sum of squares for the intermediate components.
 *     the definitions of small, intermediate and large components
 *     depend on two constants, rdwarf and rgiant. the main
 *     restrictions on these constants are that rdwarf**2 not
 *     underflow and rgiant**2 not overflow. the constants
 *     given here are suitable for every known computer.
 *
 *     the function statement is
 *
 *       double precision function enorm(n,x)
 *
 *     where
 *
 *       n is a positive integer input variable.
 *
 *       x is an input array of length n.
 *
 *     subprograms called
 * 
 *       fortran-supplied ... dabs,dsqrt
 *
 *     argonne national laboratory. minpack project. march 1980.
 *     burton s. garbow, kenneth e. hillstrom, jorge j. more
 */

double enorm_(const int *n, double const *x)
{
#ifdef USE_BLAS_NRM2
	int incx = 1;
	return dnrm2_(n, x, &incx);
#else
   double rdwarf = 3.834e-20;
   double rgiant = 1.304e19;

   /* Parameter adjustments */
   --x;

   /* Function Body */
   double s1 = 0;
   double s2 = 0;
   double s3 = 0;
   double x1max = 0;
   double x3max = 0;
   double floatn = (double) (*n);
   double agiant = rgiant / floatn;
   for (int i = 1; i <= *n; ++i) {
       double xabs = fabs(x[i]);
       if (xabs > rdwarf && xabs < agiant) {
           /* sum for intermediate components. */
           s2 += xabs * xabs;
           continue;
       }
       if (xabs <= rdwarf) {
           /* sum for small components. */
           if (xabs > x3max) {
               double d3 = x3max / xabs;
               s3 = 1 + s3 * (d3 * d3);
               x3max = xabs;
               continue;
           }
           if (xabs != 0) {
               double d4 = xabs / x3max;
               s3 += d4 * d4;
           }
           continue;
       }

       /* sum for large components. */
       if (xabs <= x1max) {
           double d2 = xabs / x1max;
           s1 += d2 * d2;
           continue;
       }
       double d1 = x1max / xabs;
       s1 = 1 + s1 * (d1 * d1);
       x1max = xabs;
   }

   /* calculation of norm. */
   if (s1 == 0) {
       if (s2 == 0)
           return x3max * sqrt(s3);
       if (s2 >= x3max)
           return sqrt(s2 * (1 + x3max / s2 * (x3max * s3)));
       if (s2 < x3max)
           return sqrt(x3max * (s2 / x3max + x3max * s3));
   }
   return x1max * sqrt(s1 + s2 / x1max / x1max);
#endif // USE_BLAS
}