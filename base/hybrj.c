#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <assert.h>

#include "minpack.h"

/*
 *     Subroutine hybrj
 *
 *     The purpose of hybrj is to find a zero of a system of 
 *     n nonlinear functions in n variables by a modification 
 *     of the Powell hybrid method.  The user must provide a 
 *     subroutine which calculates the functions and the jacobian. 
 *
 *     The subroutine statement is 
 *
 *       subroutine hybrj(fcn,n,x,fvec,fjac,ldfjac,xtol,maxfev,diag, 
 *                        mode,factor,nprint,info,nfev,njev,r,lr,qtf, 
 *                        wa1,wa2,wa3,wa4) 
 *
 *     Where 
 *
 *       fcn is the name of the user-supplied subroutine which 
 *         calculates the functions and the jacobian.  fcn must 
 *         be declared in an external statement in the user 
 *         calling program, and should be written as follows. 
 *
 *         subroutine fcn(n,x,fvec,fjac,ldfjac,iflag) 
 *         integer n,ldfjac,iflag 
 *         double precision x(n),fvec(n),fjac(ldfjac,n) 
 *         ---------- 
 *         If iflag = 1 calculate the functions at x and 
 *         return this vector in fvec.  Do not alter fjac. 
 *         If iflag = 2 calculate the jacobian at x and 
 *         return this matrix in fjac.  Do not alter fvec. 
 *         --------- 
 *         return 
 *         end 
 *
 *         The value of iflag should not be changed by fcn unless 
 *         the user wants to terminate execution of hybrj. 
 *         In this case set iflag to a negative integer. 
 *
 *       n is a positive integer input variable set to the number 
 *         of functions and variables. 
 *
 *       x is an array of length n.  On input x must contain 
 *         an initial estimate of the solution vector.  On output x 
 *         contains the final estimate of the solution vector. 
 *
 *       fvec is an output array of length n which contains 
 *         the functions evaluated at the output x. 
 *
 *       fjac is an output n by n array which contains the 
 *         orthogonal matrix q produced by the QR factorization 
 *         of the final approximate jacobian. 
 *
 *       ldfjac is a positive integer input variable not less than n 
 *         which specifies the leading dimension of the array fjac. 
 *
 *       xtol is a nonnegative input variable.  Termination 
 *         occurs when the relative error between two consecutive 
 *         iterates is at most xtol. 
 *
 *       maxfev is a positive integer input variable.  Termination 
 *         occurs when the number of calls to fcn with iflag = 1 
 *         has reached maxfev. 
 *
 *       diag is an array of length n.  If mode = 1 (see 
 *         below), diag is internally set.  If mode = 2, diag 
 *         must contain positive entries that serve as 
 *         multiplicative scale factors for the variables. 
 *
 *       mode is an integer input variable.  If mode = 1, the 
 *         variables will be scaled internally.  If mode = 2, 
 *         the scaling is specified by the input diag.  Other 
 *         values of mode are equivalent to mode = 1. 
 *
 *       factor is a positive input variable used in determining the 
 *         initial step bound.  This bound is set to the product of 
 *         factor and the euclidean norm of diag*x if nonzero, or else 
 *         to factor itself. in most cases factor should lie in the 
 *         interval (0.1, 100). 100 is a generally recommended value. 
 *
 *       nprint is an integer input variable that enables controlled 
 *         printing of iterates if it is positive.  In this case, 
 *         fcn is called with iflag = 0 at the beginning of the first 
 *         iteration and every nprint iterations thereafter and 
 *         immediately prior to return, with x and fvec available 
 *         for printing.  fvec and fjac should not be altered. 
 *         if nprint is not positive, no special calls of fcn 
 *         with iflag = 0 are made. 
 *
 *       info is an integer output variable.  If the user has 
 *         terminated execution, info is set to the (negative) 
 *         value of iflag. see description of fcn.  Otherwise, 
 *         info is set as follows. 
 *
 *         info = 0   improper input parameters. 
 *
 *         info = 1   relative error between two consecutive iterates 
 *                    is at most xtol. 
 *
 *         info = 2   number of calls to fcn with iflag = 1 has 
 *                    reached maxfev. 
 *
 *         info = 3   xtol is too small.  No further improvement in 
 *                    the approximate solution x is possible. 
 *
 *         info = 4   iteration is not making good progress, as 
 *                    measured by the improvement from the last 
 *                    five jacobian evaluations. 
 *
 *         info = 5   iteration is not making good progress, as 
 *                    measured by the improvement from the last 
 *                    ten iterations. 
 *
 *       nfev is an integer output variable set to the number of 
 *         calls to fcn with iflag = 1. 
 *
 *       njev is an integer output variable set to the number of 
 *         calls to fcn with iflag = 2. 
 *
 *       r is an output array of length lr which contains the 
 *         upper triangular matrix produced by the QR factorization 
 *         of the final approximate jacobian, stored rowwise. 
 *
 *       lr is a positive integer input variable not less than 
 *         (n*(n+1))/2. 
 *
 *       qtf is an output array of length n which contains 
 *         the vector (q transpose)*fvec. 
 *
 *       wa1, wa2, wa3, and wa4 are work arrays of length n. 
 * 
 *     Argonne National Laboratory.  MINPACK project.  March 1980. 
 *     Burton S. Garbow, Kenneth E. Hillstrom, Jorge J. Moré 
 */
void hybrj_(minpack_func_nj fcn, const int *n, double *x, double *fvec,
	   double *fjac, const int *ldfjac, const double *xtol, const int *maxfev,
	   double *diag, const int *mode, const double *factor, const int *nprint,
	   int *info, int *nfev, int *njev, double *r, const int *lr, double *qtf, 
	   double *wa1, double *wa2, double *wa3, double *wa4)
{
	/* Parameter adjustments */
	--wa4;
	--wa3;
	--wa2;
	--wa1;
	--qtf;
	--diag;
	--fvec;
	--x;
	int fjac_dim1 = *ldfjac;
	int fjac_offset = 1 + fjac_dim1;
	fjac -= fjac_offset;
	--r;

	/* epsmch is the machine precision. */
	double epsmch = DBL_EPSILON;

	int iflag = 0;
	*info = 0;
	*nfev = 0;
	*njev = 0;

	/* check the input parameters for errors. */
	if (*n <= 0 || *ldfjac < *n || *xtol < 0 || *maxfev <= 0 || *factor <= 0 || *lr < *n * (*n + 1) / 2)
		goto fini;
	if (*mode == 2) {
		for (int j = 1; j <= *n; ++j)
			if (diag[j] <= 0)
				goto fini;
	}

	/* evaluate the function at the starting point */
	/* and calculate its norm. */
	iflag = 1;
	(*fcn) (n, &x[1], &fvec[1], &fjac[fjac_offset], ldfjac, &iflag);
	*nfev = 1;
	if (iflag < 0)
		goto fini;
	double fnorm = enorm_(n, &fvec[1]);

	/* initialize iteration counter and monitors. */
	int iter = 1;
	int ncsuc = 0;
	int ncfail = 0;
	int nslow1 = 0;
	int nslow2 = 0;

	/* beginning of the outer loop. */
	for (;;) {
		int jeval = 1;

		/* calculate the jacobian matrix. */
		iflag = 2;
		(*fcn) (n, &x[1], &fvec[1], &fjac[fjac_offset], ldfjac, &iflag);
		*njev += 1;
		if (iflag < 0)
			goto fini;

		/* on the first iteration and if mode is 1, scale according */
		/* to the norms of the columns of the initial jacobian. */
		double xnorm, delta;
		if (iter == 1) {
			if (*mode != 2)
				for (int j = 1; j <= *n; ++j) {
					diag[j] = wa2[j];
					if (wa2[j] == 0)
						diag[j] = 1;
				}

			/* on the first iteration, calculate the norm of the scaled x */
			/* and initialize the step bound delta. */
			for (int j = 1; j <= *n; ++j)
				wa3[j] = diag[j] * x[j];
			xnorm = enorm_(n, &wa3[1]);
			delta = *factor * xnorm;
			if (delta == 0)
				delta = *factor;
		}

#ifdef USE_LAPACK
		double *tau = wa1;

		/* query optimal size of work */
		int lapack_info = 0;
		int lwork = -1;
		dgeqrf_(n, n, &fjac[fjac_offset], ldfjac, &tau[1], &tau[1], &lwork, &lapack_info);	// LAPACK
		assert(lapack_info == 0);
		lwork =	tau[1];
		assert(lwork >= *n);
	
		/* alloc work area. TODO: move to start of function */
		double *work = malloc(lwork * sizeof(*work));
		assert(work != NULL);
	
		/* compute the QR factorization of the jacobian. */
		dgeqrf_(n, n, &fjac[fjac_offset], ldfjac, &tau[1], work, &lwork, &lapack_info);	// LAPACK
		assert(lapack_info == 0);

		/* qtf <-- (Q transpose)*fvec */
		for (int i = 1; i <= *n; ++i)
			wa4[i] = fvec[i];
		int c1 = 1;
		dormqr_("Left", "Transpose", n, &c1, n, &fjac[fjac_offset], ldfjac, &tau[1], &wa4[1], n, work, &lwork, &lapack_info);
		assert(lapack_info == 0);
		
		for (int j = 1; j <= *n; ++j)
			qtf[j] = wa4[j];

		/* copy the triangular factor of the QR factorization into R. */
		int sing = 0;
		for (int j = 1; j <= *n; ++j) {
			int l = j;
			for (int i = 1; i <= j - 1; ++i) {
				r[l] = fjac[i + j * fjac_dim1];
				l = l + *n - i;
			}
			r[l] = fjac[j + j * fjac_dim1];
			if (r[l] == 0)
				sing = 1;
		}

		/* accumulate the orthogonal factor in fjac. */
		dorgqr_(n, n, n, &fjac[fjac_offset], ldfjac, &tau[1], work, &lwork, &lapack_info);
		assert(lapack_info == 0);

		/* TODO: move to fini */
		free(work);
#else
		/* compute the qr factorization of the jacobian. */
		int pivot = 0;
		int c1 = 1;
		qrfac_(n, n, &fjac[fjac_offset], ldfjac, &pivot, NULL, &c1, &wa1[1], &wa2[1], &wa3[1]);

		/* form (Q transpose)*fvec and store in qtf. */
		for (int i = 1; i <= *n; ++i)
			qtf[i] = fvec[i];
		for (int j = 1; j <= *n; ++j) {
			if (fjac[j + j * fjac_dim1] != 0) {
				double sum = 0;
				for (int i = j; i <= *n; ++i)
					sum += fjac[i + j * fjac_dim1] * qtf[i];
				double temp = -sum / fjac[j + j * fjac_dim1];
				for (int i = j; i <= *n; ++i)
					qtf[i] += fjac[i + j * fjac_dim1] * temp;
			}
		}
		
		/* copy the triangular factor of the QR factorization into R. */
		int sing = 0;
		for (int j = 1; j <= *n; ++j) {
			int l = j;
			for (int i = 1; i <= j - 1; ++i) {
				r[l] = fjac[i + j * fjac_dim1];
				l = l + *n - i;
			}
			r[l] = wa1[j];
			if (r[l] == 0)
				sing = 1;
		}
		/* accumulate the orthogonal factor in fjac. */
		qform_(n, n, &fjac[fjac_offset], ldfjac, &wa1[1]);
#endif

		/* rescale if necessary. */
		if (*mode != 2)
			for (int j = 1; j <= *n; ++j)
				diag[j] = fmax(diag[j], wa2[j]);

		/* beginning of the inner loop. */
		for (;;) {

			/* if requested, call fcn to enable printing of iterates. */
			if (*nprint > 0) {
				iflag = 0;
				if ((iter - 1) % *nprint == 0)
					(*fcn) (n, &x[1], &fvec[1], &fjac[fjac_offset], ldfjac, &iflag);
				if (iflag < 0)
					goto fini;
			}

			/* determine the direction p. */
			dogleg_(n, &r[1], lr, &diag[1], &qtf[1], &delta, &wa1[1], &wa2[1], &wa3[1]);

			/* store the direction p and x + p. calculate the norm of p. */
			for (int j = 1; j <= *n; ++j) {
				wa1[j] = -wa1[j];
				wa2[j] = x[j] + wa1[j];
				wa3[j] = diag[j] * wa1[j];
			}
			double pnorm = enorm_(n, &wa3[1]);

			/* on the first iteration, adjust the initial step bound. */
			if (iter == 1)
				delta = fmin(delta, pnorm);

			/* evaluate the function at x + p and calculate its norm. */
			iflag = 1;
			(*fcn) (n, &wa2[1], &wa4[1], &fjac[fjac_offset], ldfjac, &iflag);
			++(*nfev);
			if (iflag < 0)
				goto fini;
			double fnorm1 = enorm_(n, &wa4[1]);

			/* compute the scaled actual reduction. */
			double actred = -1;
			if (fnorm1 < fnorm) {
				double d = fnorm1 / fnorm;
				actred = 1 - d * d;
			}

			/* compute the scaled predicted reduction. */
			int l = 1;
			for (int i = 1; i <= *n; ++i) {
				double sum = 0;
				for (int j = i; j <= *n; ++j) {
					sum += r[l] * wa1[j];
					++l;
				}
				wa3[i] = qtf[i] + sum;
			}
			double temp = enorm_(n, &wa3[1]);
			double prered = 0;
			if (temp < fnorm) {
				double d = temp / fnorm;
				prered = 1 - d * d;
			}

			/* compute the ratio of the actual to the predicted reduction. */
			double ratio = 0;
			if (prered > 0)
				ratio = actred / prered;

			/* update the step bound. */
			if (ratio < 0.1) {
				ncsuc = 0;
				++ncfail;
				delta = 0.5 * delta;
			} else {
				ncfail = 0;
				++ncsuc;
				if (ratio >= 0.5 || ncsuc > 1)
					delta = fmax(delta, pnorm / 0.5);
				if (fabs(ratio - 1) <= 0.1)
					delta = pnorm / 0.5;
			}

			/* test for successful iteration. */
			if (ratio >= 0.0001) {
				/* successful iteration. update x, fvec, and their norms. */
				for (int j = 1; j <= *n; ++j) {
					x[j] = wa2[j];
					wa2[j] = diag[j] * x[j];
					fvec[j] = wa4[j];
				}
				xnorm = enorm_(n, &wa2[1]);
				fnorm = fnorm1;
				++iter;
			}

			/* determine the progress of the iteration. */
			++nslow1;
			if (actred >= 0.001)
				nslow1 = 0;
			if (jeval)
				++nslow2;
			if (actred >= 0.1)
				nslow2 = 0;

			/* test for convergence. */
			if (delta <= *xtol * xnorm || fnorm == 0)
				*info = 1;
			if (*info != 0)
				goto fini;

			/* tests for termination and stringent tolerances. */
			if (*nfev >= *maxfev)
				*info = 2;
			if (0.1 * fmax(0.1 * delta, pnorm) <= epsmch * xnorm)
				*info = 3;
			if (nslow2 == 5)
				*info = 4;
			if (nslow1 == 10)
				*info = 5;
			if (*info != 0)
				goto fini;

			/* criterion for recalculating jacobian. */
			if (ncfail == 2)
				break;	/* exit the inner loop */

			/* calculate the rank one modification to the jacobian */
			/* and update qtf if necessary. */
			for (int j = 1; j <= *n; ++j) {
				double sum = 0;
				for (int i = 1; i <= *n; ++i)
					sum += fjac[i + j * fjac_dim1] * wa4[i];
				wa2[j] = (sum - wa3[j]) / pnorm;
				wa1[j] = diag[j] * (diag[j] * wa1[j] / pnorm);
				if (ratio >= 0.0001)
					qtf[j] = sum;
			}

			/* compute the QR factorization of the updated jacobian. */
			int c1 = 1;
			r1updt_(n, n, &r[1], lr, &wa1[1], &wa2[1], &wa3[1], &sing);
			r1mpyq_(n, n, &fjac[fjac_offset], ldfjac, &wa2[1], &wa3[1]);
			r1mpyq_(&c1, n, &qtf[1], &c1, &wa2[1], &wa3[1]);
			jeval = 0;
		}	/* inner loop. */
	}	/* outer loop. */
 
 fini:
	/* termination, either normal or user imposed. */
	if (iflag < 0)
		*info = iflag;
	iflag = 0;
	if (*nprint > 0)
		(*fcn) (n, &x[1], &fvec[1], &fjac[fjac_offset], ldfjac, &iflag);
}
