#include <math.h>
#include <stdlib.h>
#include <stddef.h>
#include <cblas.h>

#include <mpi.h>

#include "pminpack.h"

static int query_LAPACK_worksize(int * pfjac_desc, int * fvec_desc)
{
	int lwork = -1;
	double work[1];
	int info;

	/* query scaLAPACK */
	int m = pfjac_desc[M_];
	int n = pfjac_desc[N_];
	double *pfjac = NULL;
	int *ipiv = NULL;
	double *tau = NULL;
	info = scalapack_pdgeqpf(m, n, pfjac, 1, 1, pfjac_desc, ipiv, tau, work, lwork);
	if (info != 0)
		return -1;
	int needed_dgeqp3 = work[0];
	
	double *fvec = NULL;
	info = scalapack_pdormqr("Left", "Transpose", m, 1, n, pfjac, 1, 1, pfjac_desc, tau, 
		fvec, 1, 1, fvec_desc, work, lwork);
	if (info != 0)
		return -1;
	int needed_dormqr = work[0];

	/* take max */
	int max = needed_dgeqp3 > needed_dormqr ? needed_dgeqp3 : needed_dormqr;
	return max;
}


static int jac_qrfac(pminpack_func_mn fcn, void *farg, double *x, double *fvec, double *pfvec, int *pfvec_desc,
	double *pfjac, int * pfjac_desc, int *ipvt, double *R, double *qtf, 
                 int *nfev, double *tau, double *wa4, double *work, int lwork, int ctx)
{
	int nprow, npcol, myrow, mycol;
	Cblacs_gridinfo(ctx, &nprow, &npcol, &myrow, &mycol);
	
	int mb = pfjac_desc[MB_];
	int nb = pfjac_desc[NB_];
	int lld = pfjac_desc[LLD_];
	int n = pfjac_desc[N_];
	int m = pfjac_desc[M_];

	/* compute the jacobian by forward-difference approximation and store it
	   in 2D block-cyclic distribution */
	double eps = sqrt(MINPACK_EPSILON);
	for (int j = 0; j < n; ++j) {
		double temp = x[j];
		double h = eps * fabs(x[j]);
		if (h == 0)
			h = eps;
		x[j] += h;
		(*fcn)(farg, m, n, x, wa4);
		x[j] = temp;              // restore x[j]
		for (int i = 0; i < m; i++)
			wa4[i] = (wa4[i] - fvec[i]) / h;
		extrablacs_dgeld2d(wa4, m, pfjac, j, pfjac_desc);
	}
	*nfev += n;

	/* Compute the QR factorization of the jacobian. */
	int pfjac_ncol = scalapack_numroc(n, nb, mycol, 0, npcol);
	int local_ipvt[pfjac_ncol];

	int lapack_info = scalapack_pdgeqpf(m, n, pfjac, 1, 1, pfjac_desc, local_ipvt, tau, work, lwork);
	if (lapack_info != 0)
		return -1;

	/* in LAPACK    : On exit, if JPVT(J)=K, then the J-th column of A*P was the
                          the K-th column of A.
           in ScaLAPACK : On exit, if IPIV(I) = K, the local i-th column of sub( A )*P
*          was the global K-th column of sub( A )
        */

	/* distribute fvec to pfvec */
	extrablacs_dgeld2d(fvec, m, pfvec, 0, pfvec_desc);

	/* Compute qtf <-- (Q transpose)*fvec */	
	lapack_info = scalapack_pdormqr("Left", "Transpose", m, 1, n, pfjac, 1, 1, pfjac_desc, tau,
                                         pfvec, 1, 1, pfvec_desc, work, lwork);
	if (lapack_info != 0)
		return -1;

	/* recover qtf */
	// for (int j = 0; j < n; ++j)
	// 	qtf[j] = wa4[j];
	extrablacs_dgedl2d(n, 1, pfvec, 1, 1, pfvec_desc, qtf, n);

	for (int i = 0; i < n; i++)
		printf("\tpLMDIF: process (%d, %d) qtf[%d] = %f\n", myrow, mycol, i, qtf[i]); 
	fflush(stdout);

	/* copy R --- overkill, the wanted matrix is triangular */
	//for (int j = 0; j < n; j++)
	//	for (int i = 0; i <= j; i++)
	//		R[j * n + i] = pfjac[j * lld + i];
	extrablacs_dgedl2d(n, n, pfjac, 1, 1, pfjac_desc, R, n);

	/* recover ipvt */
	int ipvt_desc[9];
	lapack_info = scalapack_descinit(ipvt_desc, 1, n, mb, nb, 0, 0, ctx, 1);
	extrablacs_igedl2d(1, n, local_ipvt, 1, 1, ipvt_desc, ipvt, 1);

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
        double *pfjac = NULL;
        double *R = NULL;
        double *wa4 = NULL;
        double *pfvec = NULL;
        double diag[n], qtf[n], wa1[n], wa2[n], wa3[n];
        int ipvt[n];
	double fnorm, delta, xnorm, gnorm;

	/* check the input parameters for errors. */
	if (n <= 0 || m < n || ftol < 0 || xtol < 0 || gtol < 0 || maxfev <= 0)
		goto fini;

	/* setup pfjac and pfvec */	
	int rank, nprocs;
	int nprow, npcol, myrow, mycol;
	int nb = 5;
	int mb = 5;
	Cblacs_pinfo(&rank, &nprocs);
	Cblacs_gridinfo(ictx, &nprow, &npcol, &myrow, &mycol);
	int pfjac_nrow = scalapack_numroc(m, mb, myrow, 0, nprow);
	int pfjac_ncol = scalapack_numroc(n, nb, mycol, 0, npcol);
	if (pfjac_nrow == 0)
		pfjac_nrow = 1; 

	if (rank == 0) {
		printf("pLMDIF: %d functions in %d variables on a %d x %d process grid\n", m, n, nprow, npcol);
		char fjacB[16];
		human_format(fjacB, (long) 8 * n * m);
		printf("pLMDIF: full jacobian is %sB\n", fjacB);
	}
	char pfjacB[16];
	human_format(pfjacB, (long) 8 * pfjac_ncol * pfjac_nrow);
	printf("\tpLMDIF: process (%d, %d) owns local matrix of size %d x %d (%sB)\n", 
		myrow, mycol, pfjac_nrow, pfjac_ncol, pfjacB);
	fflush(stdout);

	int pfjac_desc[9];
	info = scalapack_descinit(pfjac_desc, m, n, mb, nb, 0, 0, ictx, pfjac_nrow);
	if (info < 0)
		goto fini;

	int pfvec_desc[9];
	info = scalapack_descinit(pfvec_desc, m, 1, mb, nb, 0, 0, ictx, m);
	if (info < 0)
		goto fini;

        /* allocate memory */
        int lwork = query_LAPACK_worksize(pfjac_desc, pfvec_desc);
	printf("\tpLMDIF: process (%d, %d) scalapack wants %d double for scratch\n", myrow, mycol, lwork); 
	fflush(stdout);
        if (lwork < 0)
                goto fini;
       
       	pfvec = malloc((long) pfjac_nrow * sizeof(*pfvec));
        pfjac = malloc((long) pfjac_nrow * pfjac_ncol * sizeof(*pfjac));
	work = malloc((m + lwork) * sizeof(*work));     // extra space for wa4
        R = malloc(n * n * sizeof(*R));
        wa4 = malloc(m * sizeof(*wa4));
	if (pfvec == NULL || pfjac == NULL || work == NULL || R == NULL || wa4 == NULL)
		goto fini;

	/* evaluate the function at the starting point and calculate its norm */
	*nfev = 1;
	(*fcn) (farg, m, n, x, fvec);
	fnorm = enorm(m, fvec);

	/* initialize levenberg-marquardt parameter and iteration counter */
	double par = 0;
	int iter = 1;

	fflush(stdout);
	/* outer loop */
	for (;;) {
		MPI_Barrier(MPI_COMM_WORLD);
		printf("LMDIF: process (%d, %d) begin outer iteration. nfev = %d / %d\n", myrow, mycol, *nfev, maxfev);

		int lres = jac_qrfac(fcn, farg, x, fvec, pfvec, pfvec_desc, pfjac, pfjac_desc, ipvt, R, qtf, nfev, wa1, wa4, work, lwork, ictx);
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

		printf("LMDIF: process (%d, %d) gnorm = %f, xnorm = %f, fnorm = %f\n", myrow, mycol, gnorm, xnorm, fnorm);

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
			printf("LMDIF: process (%d, %d) begin inner iteration\n", myrow, mycol);

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
		
				if (rank == 0)
					printf("pLMDIF: succesful iteration. fnorm = %f\n", fnorm);
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
 	MPI_Barrier(MPI_COMM_WORLD);
	printf("pLMDIF: rank %d, finished (info=%d)\n", rank, info);

	/* termination, either normal or user imposed. */
	free(work);
        free(pfjac);
        free(R);
        free(wa4);
	return info;
}
