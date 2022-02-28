#include <float.h>
#include <math.h>

#include "minpack.h"

/*     Subroutine qrfac
 *
 *     This subroutine uses householder transformations with column 
 *     pivoting (optional) to compute a QR factorization of the 
 *     m by n matrix A.  That is, qrfac determines an orthogonal 
 *     matrix Q, a permutation matrix P, and an upper trapezoidal 
 *     matrix R with diagonal elements of nonincreasing magnitude, 
 *     such that A*P = Q*R.  The householder transformation for 
 *     column k, k = 1,2,...,min(m,n), is of the form 
 *
 *                           t 
 *           i - (1/u(k))*u*u 
 *
 *     where u has zeros in the first k-1 positions.  The form of 
 *     this transformation and the method of pivoting first 
 *     appeared in the corresponding linpack subroutine. 
 *
 *     The subroutine statement is 
 *
 *       subroutine qrfac(m,n,a,lda,pivot,ipvt,lipvt,rdiag,acnorm,wa) 
 *
 *     where 
 *
 *       m is a positive integer input variable set to the number 
 *         of rows of A. 
 *
 *       n is a positive integer input variable set to the number 
 *         of columns of A. 
 *
 *       A is an m by n array.  On input a contains the matrix for 
 *         which the QR factorization is to be computed.  On output 
 *         the strict upper trapezoidal part of a contains the strict 
 *         upper trapezoidal part of R, and the lower trapezoidal 
 *         part of a contains a factored form of Q (the non-trivial 
 *         elements of the u vectors described above). 
 *
 *       lda is a positive integer input variable not less than m 
 *         which specifies the leading dimension of the array A. 
 *
 *       pivot is a logical input variable.  If pivot is set true, 
 *         then column pivoting is enforced. If pivot is set false, 
 *         then no column pivoting is done. 
 *
 *       ipvt is an integer output array of length lipvt.  ipvt 
 *         defines the permutation matrix P such that A*P = Q*R. 
 *         column j of P is column ipvt(j) of the identity matrix. 
 *         If pivot is false, ipvt is not referenced. 
 *
 *       lipvt is a positive integer input variable.  If pivot is false, 
 *         then lipvt may be as small as 1.  If pivot is true, then 
 *         lipvt must be at least n.
 *
 *       rdiag is an output array of length n which contains the 
 *         diagonal elements of R.
 *
 *       acnorm is an output array of length n which contains the 
 *         norms of the corresponding columns of the input matrix A. 
 *         If this information is not needed, then acnorm can coincide 
 *         with rdiag. 
 *
 *       wa is a work array of length n.  If pivot is false, then wa 
 *         can coincide with rdiag. 
 *
 *     Argonne national laboratory.  MINPACK project.  March 1980. 
 *     Burton s. Garbow, Kenneth e. Hillstrom, Jorge j. Moré 
 */
void qrfac_(const int *m, const int *n, double *a, const int *lda, int *pivot, int *ipvt, const int *lipvt, double *rdiag, double *acnorm, double *wa)
{
	(void)lipvt;
	/* Parameter adjustments */
	--wa;
	--acnorm;
	--rdiag;
	int a_dim1 = *lda;
	a -= 1 + a_dim1;
	--ipvt;

	/* epsmch is the machine precision. */
	double epsmch = DBL_EPSILON;

	/* compute the initial column norms and initialize several arrays. */
	for (int j = 1; j <= *n; ++j) {
		acnorm[j] = enorm_(m, &a[j * a_dim1 + 1]);
		rdiag[j] = acnorm[j];
		wa[j] = rdiag[j];
		if (*pivot)
			ipvt[j] = j;
	}

	/* reduce A to R with householder transformations. */
	int minmn = (*n < *m) ? *n : *m;
	for (int j = 1; j <= minmn; ++j) {
		if (*pivot) {
			/* bring the column of largest norm into the pivot position. */
			int kmax = j;
			for (int k = j; k <= *n; ++k)
				if (rdiag[k] > rdiag[kmax])
					kmax = k;
			if (kmax != j) {
				for (int i = 1; i <= *m; ++i) {
					double temp = a[i + j * a_dim1];
					a[i + j * a_dim1] = a[i + kmax * a_dim1];
					a[i + kmax * a_dim1] = temp;
				}
				rdiag[kmax] = rdiag[j];
				wa[kmax] = wa[j];
				int k = ipvt[j];
				ipvt[j] = ipvt[kmax];
				ipvt[kmax] = k;
			}
 		}

		/* Compute the householder transformation to reduce the */
		/* j-th column of a to a multiple of the j-th unit vector. */
		int i2 = *m - j + 1;
		double ajnorm = enorm_(&i2, &a[j + j * a_dim1]);
		if (ajnorm > 0) {
			if (a[j + j * a_dim1] < 0)
				ajnorm = -ajnorm;
			for (int i = j; i <= *m; ++i)
				a[i + j * a_dim1] /= ajnorm;
			a[j + j * a_dim1] += 1;
	
			/* Apply the transformation to the remaining columns */
			/* and update the norms. */
			for (int k = j + 1; k <= *n; ++k) {
				double sum = 0;
				for (int i = j; i <= *m; ++i)
					sum += a[i + j * a_dim1] * a[i + k * a_dim1];
				double temp = sum / a[j + j * a_dim1];
				for (int i = j; i <= *m; ++i)
					a[i + k * a_dim1] -= temp * a[i + j * a_dim1];
				if (!(*pivot) || rdiag[k] == 0)
					continue;
				double tmp2 = a[j + k * a_dim1] / rdiag[k];
				rdiag[k] *= sqrt(fmax(0, 1 - tmp2*tmp2));
				double d = rdiag[k] / wa[k];
				if (0.05 * (d * d) > epsmch)
					continue;
				int i3 = *m - j;
				rdiag[k] = enorm_(&i3, &a[j+1 + k * a_dim1]);
				wa[k] = rdiag[k];
			}
		}
		rdiag[j] = -ajnorm;
	}
}
