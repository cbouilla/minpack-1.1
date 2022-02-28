#include "minpack.h"
/*
 *     subroutine qform 
 *
 *     this subroutine proceeds from the computed QR factorization of 
 *     an m by n matrix A to accumulate the m by m orthogonal matrix 
 *     Q from its factored form. 
 *
 *     the subroutine statement is 
 *
 *       subroutine qform(m,n,q,ldq,wa) 
 *
 *     where 
 *
 *       m is a positive integer input variable set to the number 
 *         of rows of a and the order of q. 
 *
 *       n is a positive integer input variable set to the number 
 *         of columns of a. 
 *
 *       Q is an m by m array.  On input the full lower trapezoid in 
 *         the first min(m,n) columns of Q contains the factored form. 
 *         on output Q has been accumulated into a square matrix. 
 *
 *       ldq is a positive integer input variable not less than m 
 *         which specifies the leading dimension of the array Q. 
 *
 *       wa is a work array of length m. 
 * 
 *     argonne national laboratory. minpack project. march 1980. 
 *     burton s. garbow, kenneth e. hillstrom, jorge j. more 
 */ 
void qform_(const int * m, const int * n, double * q, const int *ldq, double * wa)
{
	/* Parameter adjustments */
	--wa;
	int q_dim1 = *ldq;
	int q_offset = 1 + q_dim1 * 1;
	q -= q_offset;

	/* zero out upper triangle of q in the first min(m,n) columns. */
	int minmn = (*n < *m) ? *n : *m;
	for (int j = 2; j <= minmn; ++j)
		for (int i = 1; i <= j - 1; ++i)
			q[i + j * q_dim1] = 0;

	/* initialize remaining columns to those of the identity matrix. */
	for (int j = *n + 1; j <= *m; ++j) {
		for (int i = 1; i <= *m; ++i)
			q[i + j * q_dim1] = 0;
		q[j + j * q_dim1] = 1;
	}

	/* Accumulate Q from its factored form. */
	for (int l = 1; l <= minmn; ++l) {
		int k = minmn - l + 1;
		for (int i = k; i <= *m; ++i) {
			wa[i] = q[i + k * q_dim1];
			q[i + k * q_dim1] = 0;
		}
		q[k + k * q_dim1] = 1;
		if (wa[k] == 0)
			return;
		for (int j = k; j <= *m; ++j) {
			double sum = 0;
			for (int i = k; i <= *m; ++i) 
				sum += q[i + j * q_dim1] * wa[i];
			double temp = sum / wa[k];
			for (int i = k; i <= *m; ++i)
				q[i + j * q_dim1] -= temp * wa[i];
		}
	}
}