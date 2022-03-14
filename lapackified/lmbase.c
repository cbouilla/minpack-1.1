#include <math.h>
#include <stdlib.h>
#include <stddef.h>
#include <cblas.h>

#include "cminpack.h"


static int query_worksize(int m, int n, int ldfjac)
{
        int c1 = 1;
        int lwork = -1;
        double work[1];
        int info;

        /* query LAPACK */
        dgeqp3_(&m, &n, NULL, &ldfjac, NULL, NULL, work, &lwork, &info);
        if (info != 0)
                return -1;
        int needed_dgeqp3 = work[0];
        if (needed_dgeqp3 < 3 * n + 1)
                return -1;
        dormqr_("Left", "Transpose", &m, &c1, &n, NULL, &ldfjac, NULL, NULL, &m, work, &lwork, &info);
        if (info != 0)
                return -1;
        int needed_dormqr = work[0];
        if (needed_dormqr < n)
                return -1;
        /* take max */
        int max = needed_dgeqp3 > needed_dormqr ? needed_dgeqp3 : needed_dormqr;
        return max;
}

/*
 * This function concentrates all the common code to lmdif and lmder
 */
int lmbase(cminpack_func_mnj fcn_der, cminpack_func_mn fcn_dif, void *farg, 
	int m, int n, double *x, double *fvec, double *fjac, int ldfjac, double ftol,
	    double xtol, double gtol, int maxfev, double epsfcn, double *diag, 
	    int mode, double factor, int nprint, int *nfev, int *njev, int *ipvt, 
	    double *qtf, double *wa1, double *wa2, double *wa3, double *wa4)
{
	/* Parameter adjustments */
	ptrdiff_t fjac_dim1 = ldfjac;
	int c1 = 1;

	/* epsmch is the machine precision. */
	double epsmch = MINPACK_EPSILON;
	int iflag = 0;
	int info = 0;
	*nfev = 0;
        double * work = NULL;

	/* check the input parameters for errors. */
	if (n <= 0 || m < n || ldfjac < m || ftol < 0 || xtol < 0 || gtol < 0
	    || maxfev <= 0 || factor <= 0)
		goto fini;
	if ((fcn_dif == NULL) && (fcn_der == NULL))
		goto fini;
	if ((fcn_dif != NULL) && (fcn_der != NULL))
		goto fini;
	if (mode == 2) {
		for (int j = 0; j < n; ++j)
			if (diag[j] <= 0)
				goto fini;
	}

        int lwork = query_worksize(m, n, ldfjac);
        if (lwork < 0)
                goto fini;
        work = malloc(lwork * sizeof(*work));
        if (work == NULL)
                goto fini;

	/* evaluate the function at the starting point and calculate its norm */
	*nfev = 1;
	if (fcn_der)
		iflag = (*fcn_der) (farg, m, n, x, fvec, fjac, ldfjac, 1);
	if (fcn_dif)
		iflag = (*fcn_dif) (farg, m, n, x, fvec, 1);
	if (iflag < 0)
		goto fini;
	double fnorm = enorm(m, fvec);
	double delta, xnorm, gnorm;

	/* initialize levenberg-marquardt parameter and iteration counter */
	double par = 0;
	int iter = 1;

	/* outer loop */
	for (;;) {
		/* calculate the jacobian matrix. */
		if (fcn_der) {
			iflag = (*fcn_der) (farg, m, n, x, fvec, fjac, ldfjac, 2);
			*njev += 1;
		}
		if (fcn_dif) {
			iflag = fdjac2(fcn_dif, farg, m, n, x, fvec, fjac, ldfjac, 2, epsfcn, wa4);
			*nfev += n;
		}
		if (iflag < 0)
			break;

		/* if requested, call fcn to enable printing of iterates. */
		if (nprint > 0) {
			if ((iter - 1) % nprint == 0) {
				if (fcn_der)
					iflag = (*fcn_der) (farg, m, n, x, fvec, fjac, ldfjac, 0);
				if (fcn_dif)
					iflag = (*fcn_dif) (farg, m, n, x, fvec, 0);
			}
			if (iflag < 0)
				break;
		}

		double *acnorm = wa2;
		for (int j = 0; j < n; ++j)
			acnorm[j] = enorm(m, &fjac[j * fjac_dim1]);

		/* on the first iteration and if mode is 1, scale according
		   to the norms of the columns of the initial jacobian. */
		if (iter == 1) {
			if (mode != 2) {
				for (int j = 0; j < n; ++j) {
					diag[j] = acnorm[j];
					if (acnorm[j] == 0.0)
						diag[j] = 1.0;
				}
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

		/* set all columns free */
		for (int j = 0; j < n; j++)
			ipvt[j] = 0;

		/* compute the QR factorization of the jacobian. */
                int lapack_info;
                double *tau = wa1;
		dgeqp3_(&m, &n, fjac, &ldfjac, ipvt, tau, work, &lwork, &lapack_info);
		if (lapack_info != 0)
                        goto fini;

		/* qtf <-- (Q transpose)*fvec */
		for (int i = 0; i < m; ++i)
			wa4[i] = fvec[i];
		dormqr_("Left", "Transpose", &m, &c1, &n, fjac, &ldfjac, tau, wa4, &m, work, &lwork, &lapack_info);
                if (lapack_info != 0)
                        goto fini;
		for (int j = 0; j < n; ++j)
			qtf[j] = wa4[j];

		/* Compute the norm of the scaled gradient. */
		gnorm = 0;
		if (fnorm != 0) {
			for (int j = 0; j < n; ++j) {
				int l = ipvt[j] - 1;
				if (acnorm[l] == 0)
					continue;
				double sum = 0;
				for (int i = 0; i <= j; ++i)
					sum += fjac[i + j * fjac_dim1] * (qtf[i] / fnorm);
				gnorm = fmax(gnorm, fabs(sum / acnorm[l]));
			}
		}

		/* test for convergence of the gradient norm. */
		if (gnorm <= gtol)
			info = 4;
		if (info != 0)
			break;

		/* rescale if necessary. */
		if (mode != 2)
			for (int j = 0; j < n; ++j)
				diag[j] = fmax(diag[j], wa2[j]);

		/* inner loop. */
		for (;;) {
			/* determine the levenberg-marquardt parameter. */
			par = lmpar(n, fjac, ldfjac, ipvt, diag, qtf, delta, wa1, wa2, wa3, wa4);

			/* store the direction p and x + p. calculate the norm of p. */
			for (int j = 0; j < n; ++j) {
				wa1[j] = -wa1[j];
				wa2[j] = x[j] + wa1[j];
				wa3[j] = diag[j] * wa1[j];
			}
			double pnorm = enorm(n, wa3);

			/* On the first iteration, adjust the initial step bound. */
			if (iter == 1)
				delta = fmin(delta, pnorm);

			/* Evaluate the function at x + p and calculate its norm. */
			if (fcn_der)
				iflag = (*fcn_der) (farg, m, n, wa2, wa4, fjac, ldfjac, 1);
			if (fcn_dif)
				iflag = (*fcn_dif) (farg, m, n, wa2, wa4, 1);
			*nfev += 1;
			if (iflag < 0)
				goto fini;
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
				wa3[j] = wa1[l];
			}
			cblas_dtrmv(CblasColMajor, CblasUpper, CblasNoTrans, CblasNonUnit, n, fjac, ldfjac, wa3, 1);
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
			if (fabs(actred) <= epsmch && prered <= epsmch && 0.5 * ratio <= 1)
				info = 6;
			if (delta <= epsmch * xnorm)
				info = 7;
			if (gnorm <= epsmch)
				info = 8;
			if (info != 0)
				goto fini;

			/* repeat if iteration unsuccessful. */
			if (ratio >= 0.0001)
				break;
		}		/* inner loop */
	}			/* outer loop. */

 fini:
        free(work);

	/* termination, either normal or user imposed. */
	if (iflag < 0)
		info = iflag;
	iflag = 0;
	if (nprint > 0) {
		if (fcn_der)
			(*fcn_der) (farg, m, n, x, fvec, fjac, ldfjac, 0);
		if (fcn_dif)
			(*fcn_dif) (farg, m, n, x, fvec, 0);
	}
	return info;
}
