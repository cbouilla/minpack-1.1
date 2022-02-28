#include <math.h>

#include "minpack.h"
/*
 *     subroutine r1mpyq 
 *
 *     given an m by n matrix a, this subroutine computes a*q where 
 *     q is the product of 2*(n - 1) transformations 
 *
 *           gv(n-1)*...*gv(1)*gw(1)*...*gw(n-1) 
 *
 *     and gv(i), gw(i) are givens rotations in the (i,n) plane which 
 *     eliminate elements in the i-th and n-th planes, respectively. 
 *     q itself is not given, rather the information to recover the 
 *     gv, gw rotations is supplied. 
 *
 *     the subroutine statement is 
 *
 *       subroutine r1mpyq(m,n,a,lda,v,w) 
 *
 *     where 
 *
 *       m is a positive integer input variable set to the number 
 *         of rows of a. 
 *
 *       n is a positive integer input variable set to the number 
 *         of columns of a. 
 *
 *       a is an m by n array. on input a must contain the matrix 
 *         to be postmultiplied by the orthogonal matrix q 
 *         described above. on output a*q has replaced a. 
 *
 *       lda is a positive integer input variable not less than m 
 *         which specifies the leading dimension of the array a. 
 *
 *       v is an input array of length n. v(i) must contain the 
 *         information necessary to recover the givens rotation gv(i) 
 *         described above. 
 *
 *       w is an input array of length n. w(i) must contain the 
 *         information necessary to recover the givens rotation gw(i) 
 *         described above. 
 *
 *     argonne national laboratory. minpack project. march 1980. 
 *     burton s. garbow, kenneth e. hillstrom, jorge j. more 
 */
void r1mpyq_(const int *m, const int *n, double *a, const int *lda, double *v, double *w)
{
	/* Parameter adjustments */
	--w;
	--v;
	int a_dim1 = *lda;
	int a_offset = 1 + a_dim1 * 1;
	a -= a_offset;

	/* apply the first set of givens rotations to a. */	
	for (int nmj = 1; nmj <= *n - 1; ++nmj) {
		int j = *n - nmj;
		double c, s;
		if (fabs(v[j]) > 1) {
			c = 1 / v[j];
			s = sqrt(1 - c * c);
		} else {
			s = v[j];
			c = sqrt(1 - s * s);
		}
		for (int i = 1; i <= *m; ++i) {
			double temp = c * a[i + j * a_dim1] - s * a[i + *n * a_dim1];
			a[i + *n * a_dim1] = s * a[i + j * a_dim1] + c * a[i + *n * a_dim1];
			a[i + j * a_dim1] = temp;
		}
	}

	/* apply the second set of givens rotations to a. */
	for (int j = 1; j <= *n - 1; ++j) {
		double c, s;
		if (fabs(w[j]) > 1) {
			c = 1 / w[j];
			s = sqrt(1 - c * c);
		} else {
			s = w[j];
			c = sqrt(1 - s * s);
		}
		for (int i = 1; i <= *m; ++i) {
			double temp = c * a[i + j * a_dim1] + s * a[i + *n * a_dim1];
			a[i + *n * a_dim1] = -s * a[i + j * a_dim1] + c * a[i + *n * a_dim1];
			a[i + j * a_dim1] = temp;
		}
	}
}