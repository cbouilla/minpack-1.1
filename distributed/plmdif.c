#include <math.h>
#include <stdlib.h>
#include <stddef.h>
#include <cblas.h>
#include <lapacke.h>

#include "pminpack.h"

static int query_LAPACK_worksize(int m, int n, int ldfjac)
{
	int lwork = -1;
	double work[1];
	int info;

	/* query LAPACK */
	info = LAPACKE_dgeqp3_work(LAPACK_COL_MAJOR, m, n, NULL, ldfjac, NULL, NULL, work, lwork);
	if (info != 0)
		return -1;
	int needed_dgeqp3 = work[0];
	if (needed_dgeqp3 < 3 * n + 1)
		return -1;

	info = LAPACKE_dormqr_work(LAPACK_COL_MAJOR, 'L', 'T', m, 1, n, NULL, ldfjac, NULL, NULL, m, work, lwork);
	if (info != 0)
		return -1;
	int needed_dormqr = work[0];
	if (needed_dormqr < n)
		return -1;
	/* take max */
	int max = needed_dgeqp3 > needed_dormqr ? needed_dgeqp3 : needed_dormqr;
	return max;
}


static int jac_qrfac(pminpack_func_mn fcn, void *farg, int n, int m, double *x, double *fvec, 
	double *fjac, ptrdiff_t ldfjac, int *ipvt, double *R, double *qtf, 
                 int *nfev, double *wa1, double *wa4, double *work, int lwork, int ictx)
{
	double eps = sqrt(MINPACK_EPSILON);

	for (int j = 0; j < n; ++j) {
		double temp = x[j];
		double h = eps * fabs(x[j]);
		if (h == 0)
			h = eps;
		x[j] += h;
		(*fcn)(farg, m, n, x, wa4);
		x[j] = temp;              // restore x[j]
		for (int i = 0; i < m; ++i)
			fjac[i + j * ldfjac] = (wa4[i] - fvec[i]) / h;
	}
	*nfev += n;

	/* set all columns free */
	for (int j = 0; j < n; j++)
		ipvt[j] = 0;

	/* compute the QR factorization of the jacobian. */
	double *tau = wa1;
	int lapack_info = LAPACKE_dgeqp3_work(LAPACK_COL_MAJOR, m, n, fjac, ldfjac, ipvt, tau, work, lwork);
	// dgeqp3_(&m, &n, fjac, &ldfjac, ipvt, tau, work, &lwork, &lapack_info);
	if (lapack_info != 0)
		return -1;

	/* qtf <-- (Q transpose)*fvec */
	for (int i = 0; i < m; ++i)
		wa4[i] = fvec[i];
	lapack_info = LAPACKE_dormqr_work(LAPACK_COL_MAJOR, 'L', 'T', m, 1, n, fjac, ldfjac, tau, wa4, m, work, lwork);
	if (lapack_info != 0)
		return -1;
	for (int j = 0; j < n; ++j)
		qtf[j] = wa4[j];

	/* copy R */
	for (int j = 0; j < n; j++)
		for (int i = 0; i <= j; i++)
			R[j * n + i] = fjac[j * ldfjac + i];
	return 0;
}

/*     subroutine plmdif
 *
 *     the purpose of plmdif is to minimize the sum of the squares of
 *     m nonlinear functions in n variables by a modification of
 *     the levenberg-marquardt algorithm. the user must provide a
 *     subroutine which calculates the functions. the jacobian is
 *     then calculated by a forward-difference approximation.
 *
 *       fcn is the name of the user-supplied subroutine which
 *         calculates the functions. fcn must be declared
 *         in an external statement in the user calling
 *         program, and should be written as follows.
 *
 *         subroutine fcn(m,n,x,fvec)
 *         int m,n
 *         double precision x(n),fvec(m)
 *         ----------
 *         calculate the functions at x and
 *         return this vector in fvec.
 *         ----------
 *         return
 *         end
 *
 *       m is a positive int input variable set to the number
 *         of functions.
 *
 *       n is a positive int input variable set to the number
 *         of variables. n must not exceed m.
 *
 *       x is an array of length n.  On input x must contain
 *         an initial estimate of the solution vector.  On output x
 *         contains the final estimate of the solution vector.
 *
 *       fvec is an output array of length m which contains
 *         the functions evaluated at the output x.
 *
 *       ftol is a nonnegative input variable.  Termination
 *         occurs when both the actual and predicted relative
 *         reductions in the sum of squares are at most ftol.
 *         Therefore, ftol measures the relative error desired
 *         in the sum of squares.
 *
 *       xtol is a nonnegative input variable.  Termination
 *         occurs when the relative error between two consecutive
 *         iterates is at most xtol.  Therefore, xtol measures the
 *         relative error desired in the approximate solution.
 *
 *       gtol is a nonnegative input variable.  Termination
 *         occurs when the cosine of the angle between fvec and
 *         any column of the jacobian is at most gtol in absolute
 *         value.  Therefore, gtol measures the orthogonality
 *         desired between the function vector and the columns
 *         of the jacobian.
 *
 *       maxfev is a positive int input variable. termination
 *         occurs when the number of calls to fcn is at least
 *         maxfev by the end of an iteration.
 *
 *       info is an int output variable. if the user has
 *         terminated execution, info is set to the (negative)
 *         value of iflag. see description of fcn. otherwise,
 *         info is set as follows.
 *
 *         info < 0  runtime error (scalapack or lapack).
 *
 *         info = 0  improper input parameters.
 *
 *         info = 1  both actual and predicted relative reductions
 *                   in the sum of squares are at most ftol.
 *
 *         info = 2  relative error between two consecutive iterates
 *                   is at most xtol.
 *
 *         info = 3  conditions for info = 1 and info = 2 both hold.
 *
 *         info = 4  the cosine of the angle between fvec and any
 *                   column of the jacobian is at most gtol in
 *                   absolute value.
 *
 *         info = 5  number of calls to fcn has reached or
 *                   exceeded maxfev.
 *
 *         info = 6  ftol is too small. no further reduction in
 *                   the sum of squares is possible.
 *
 *         info = 7  xtol is too small. no further improvement in
 *                   the approximate solution x is possible.
 *
 *         info = 8  gtol is too small. fvec is orthogonal to the
 *                   columns of the jacobian to machine precision.
 *
 *       nfev is an int output variable set to the number of
 *         calls to fcn.
 *
 *       ldfjac is a positive int input variable not less than m
 *         which specifies the leading dimension of the array fjac.
 *
 */

int plmdif(pminpack_func_mn fcn, void *farg, int m, int n, double *x, double *fvec, 
        double ftol, double xtol, double gtol, int maxfev, int *nfev, int ictx)
{
	double factor = 100.;
	int info = 0;
	*nfev = 0;

        /* this will be allocated dynamically by this function */ 
        double *work = NULL;
        double *fjac = NULL;
        double *R = NULL;
        double *wa4 = NULL;
        double diag[n], qtf[n], wa1[n], wa2[n], wa3[n];
        int ipvt[n];
	double fnorm, delta, xnorm, gnorm;

	/* check the input parameters for errors. */
	if (n <= 0 || m < n || ftol < 0 || xtol < 0 || gtol < 0 || maxfev <= 0)
		goto fini;

	/* setup pfjac */	
	int rank, nprocs;
	int nprow, npcol, myrow, mycol;
	int nb = 64;
	int mb = 64;
	Cblacs_pinfo(&rank, &nprocs);
	Cblacs_gridinfo(ictx, &nprow, &npcol, &myrow, &mycol);
	long pfjac_nrow = scalapack_numroc(m, mb, myrow, 0, nprow);
	long pfjac_ncol = scalapack_numroc(n, nb, mycol, 0, npcol);

	if (rank == 0) {
		printf("pLMDIF: %d functions in %d variables on a %d x %d process grid\n", m, n, nprow, npcol);
		char fjacB[16], pfjacB[16];
		human_format(fjacB, (long) 8 * n * m);
		human_format(pfjacB, (long) 8 * pfjac_ncol * pfjac_nrow);
		printf("pLMDIF: full jacobian is %sB. Each process holds %sB\n", fjacB, pfjacB);
	}

        /* allocate memory */
        ptrdiff_t ldfjac = m;
        int lwork = query_LAPACK_worksize(m, n, ldfjac);
        if (lwork < 0)
                goto fini;
        fjac = malloc(n * ldfjac * sizeof(*fjac));
	work = malloc((m + lwork) * sizeof(*work)); // extra space for wa4
        R = malloc(n * n * sizeof(*R));
        wa4 = malloc(m * sizeof(*wa4));
	if (fjac == NULL || work == NULL || R == NULL || wa4 == NULL)
		goto fini;

	/* evaluate the function at the starting point and calculate its norm */
	*nfev = 1;
	(*fcn) (farg, m, n, x, fvec);
	fnorm = enorm(m, fvec);

	/* initialize levenberg-marquardt parameter and iteration counter */
	double par = 0;
	int iter = 1;

	/* outer loop */
	for (;;) {
		if (rank == 0)
			printf("pLMDIF: begin outer iteration\n");

		int lres = jac_qrfac(fcn, farg, n, m, x, fvec, fjac, ldfjac, ipvt, R, qtf, nfev, wa1, wa4, work, lwork, ictx);
                if (lres < 0) {
                        info = -1;
                        goto fini;
                }

		double *acnorm = wa2;
		for (int j = 0; j < n; ++j) {
			int l = ipvt[j] - 1;
			acnorm[l] = enorm(j + 1, &R[j * n]);
		}

		/* on the first iteration and if mode is 1, scale according
		   to the norms of the columns of the initial jacobian. */
		if (iter == 1) {
			for (int j = 0; j < n; ++j) {
				diag[j] = acnorm[j];
				if (acnorm[j] == 0.0)
					diag[j] = 1.0;
			}

			/* On the first iteration, calculate the norm of the scaled x
			   and initialize the step bound delta. */
			for (int j = 0; j < n; ++j)
				wa3[j] = diag[j] * x[j];
			xnorm = enorm(n, wa3);
			delta = factor * xnorm;
			if (delta == 0)
				delta = factor;
		}

		/* Compute the norm of the scaled gradient. */
		gnorm = 0;
		if (fnorm != 0) {
			for (int j = 0; j < n; ++j) {
				int l = ipvt[j] - 1;
				if (acnorm[l] == 0)
					continue;
				double sum = 0;
				for (int i = 0; i <= j; ++i)
					sum += R[i + j * n] * (qtf[i] / fnorm);
				gnorm = fmax(gnorm, fabs(sum / acnorm[l]));
			}
		}

		/* test for convergence of the gradient norm. */
		if (gnorm <= gtol)
			info = 4;
		if (info != 0)
			break;

		/* rescale if necessary. */
		for (int j = 0; j < n; ++j)
			diag[j] = fmax(diag[j], acnorm[j]);

		/* inner loop. */
		for (;;) {
			if (rank == 0)
				printf("pLMDIF: begin inner iteration\n");

			/* determine the levenberg-marquardt parameter. */
			double *p = wa1;
			par = lmpar(n, R, n, ipvt, diag, qtf, delta, p, wa2, wa3, wa4);
			if (rank == 0)
				printf("pLMDIF: LM parameter %f\n", par);

			/* store the direction p and x + p. calculate the norm of p. */
			for (int j = 0; j < n; ++j) {
				p[j] = -p[j];
				wa2[j] = x[j] + p[j];
				wa3[j] = diag[j] * p[j];
			}
			double pnorm = enorm(n, wa3);

			/* On the first iteration, adjust the initial step bound. */
			if (iter == 1)
				delta = fmin(delta, pnorm);

			/* Evaluate the function at x + p and calculate its norm. */
			(*fcn) (farg, m, n, wa2, wa4);
			*nfev += 1;
			double fnorm1 = enorm(m, wa4);

			/* compute the scaled actual reduction. */
			double actred = -1;
			if (0.1 * fnorm1 < fnorm) {
				/* Computing 2nd power */
				double tmp = fnorm1 / fnorm;
				actred = 1.0 - tmp * tmp;
			}

			/* compute the scaled predicted reduction and
			   the scaled directional derivative. */
			for (int j = 0; j < n; ++j) {
				int l = ipvt[j] - 1;
				wa3[j] = p[l];
			}
			cblas_dtrmv(CblasColMajor, CblasUpper, CblasNoTrans, CblasNonUnit, n, R, n, wa3, 1);
			double temp1 = enorm(n, wa3) / fnorm;
			double temp2 = sqrt(par) * pnorm / fnorm;
			double prered = temp1 * temp1 + 2 * temp2 * temp2;
			double dirder = -(temp1 * temp1 + temp2 * temp2);

			/* compute the ratio of the actual to the predicted reduction. */
			double ratio = 0;
			if (prered != 0)
				ratio = actred / prered;

			/* update the step bound. */
			if (ratio > 0.25) {
				if (par == 0 || ratio >= 0.75) {
					delta = pnorm / 0.5;
					par = 0.5 * par;
				}
			} else {
				double temp;
				if (actred >= 0)
					temp = 0.5;
				else
					temp = 0.5 * dirder / (dirder + 0.5 * actred);
				if (0.1 * fnorm1 >= fnorm || temp < 0.1)
					temp = 0.1;
				delta = temp * fmin(delta, pnorm / 0.1);
				par /= temp;
			}

			/* test for successful iteration. */
			if (ratio >= 0.0001) {
				/* successful iteration. update x, fvec, and their norms. */
				for (int j = 0; j < n; ++j) {
					x[j] = wa2[j];
					wa2[j] = diag[j] * x[j];
				}
				for (int i = 0; i < m; ++i)
					fvec[i] = wa4[i];
				xnorm = enorm(n, wa2);
				fnorm = fnorm1;
				++iter;
			}

			/* tests for convergence */
			if (fabs(actred) <= ftol && prered <= ftol && 0.5 * ratio <= 1)
				info = 1;
			if (delta <= xtol * xnorm)
				info = 2;
			if (fabs(actred) <= ftol && prered <= ftol && 0.5 * ratio <= 1 && info == 2)
				info = 3;
			if (info != 0)
				goto fini;

			/* tests for termination and stringent tolerances. */
			if (*nfev >= maxfev)
				info = 5;
			if (fabs(actred) <= MINPACK_EPSILON && prered <= MINPACK_EPSILON && 0.5 * ratio <= 1)
				info = 6;
			if (delta <= MINPACK_EPSILON * xnorm)
				info = 7;
			if (gnorm <= MINPACK_EPSILON)
				info = 8;
			if (info != 0)
				goto fini;

			/* repeat if iteration unsuccessful. */
			if (ratio >= 0.0001)
				break;
		}		/* inner loop */
	}			/* outer loop. */

 fini:
	/* termination, either normal or user imposed. */
	free(work);
        free(fjac);
        free(R);
        free(wa4);
	return info;
}
