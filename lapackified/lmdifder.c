#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <assert.h>

#include "minpack.h"

/*
 * This function concentrates all the common code to lmdif and lmder
 */
void lmddifer_(minpack_func_mnj fcn_der, minpack_func_mn fcn_dif, const int *m, const int *n, double *x, 
	double *fvec, double *fjac, const int *ldfjac, const double *ftol,
	const double *xtol, const double *gtol, const int *maxfev, const double * epsfcn, 
    double *diag, const int *mode, const double *factor, const int *nprint, 
    int *info, int *nfev, int *njev, int *ipvt, double *qtf, 
	double *wa1, double *wa2, double *wa3, double *wa4)
{
    /* Parameter adjustments */
    int fjac_dim1 = *ldfjac;
        int c1 = 1;

    /* epsmch is the machine precision. */
    double epsmch = DBL_EPSILON;
    int iflag = 0;
    *info = 0;
    *nfev = 0;

    /* check the input parameters for errors. */
    if (*n <= 0 || *m < *n || *ldfjac < *m || *ftol < 0 || *xtol < 0 || *gtol < 0 || *maxfev <= 0 || *factor <= 0)
        goto fini;
    if ((fcn_dif == NULL) && (fcn_der == NULL))
        goto fini;
    if ((fcn_dif != NULL) && (fcn_der != NULL))
        goto fini;
    if (*mode == 2) {
        for (int j = 0; j < *n; ++j)
            if (diag[j] <= 0)
                goto fini;
    }

    /* evaluate the function at the starting point and calculate its norm */
    *nfev = 1;
    iflag = 1;
    if (fcn_der)
        (*fcn_der)(m, n, x, fvec, fjac, ldfjac, &iflag);
    if (fcn_dif)
        (*fcn_dif)(m, n, x, fvec, &iflag);
    if (iflag < 0)
        goto fini;
    double fnorm = enorm_(m, fvec);
    double delta, xnorm, gnorm;

    /* initialize levenberg-marquardt parameter and iteration counter */
    double par = 0;
    int iter = 1;

    /* outer loop */
    for(;;) {
        /* calculate the jacobian matrix. */
        iflag = 2;
        if (fcn_der) {
            (*fcn_der)(m, n, x, fvec, fjac, ldfjac, &iflag);
            *njev += 1;
        }
        if (fcn_dif) {
            fdjac2_(fcn_dif, m, n, x, fvec, fjac, ldfjac, &iflag, epsfcn, wa4);
            *nfev += *n;
        }
        if (iflag < 0)
            break;

        /* if requested, call fcn to enable printing of iterates. */
        if (*nprint > 0) {
            iflag = 0;
            if ((iter - 1) % *nprint == 0) {
                if (fcn_der)
                    (*fcn_der)(m, n, x, fvec, fjac, ldfjac, &iflag);
                if (fcn_dif)
                    (*fcn_dif)(m, n, x, fvec, &iflag);
            }
            if (iflag < 0)
                break;
        }

        double *acnorm = wa2;
        for (int j = 0; j < *n; ++j)
            acnorm[j] = enorm_(m, &fjac[j * fjac_dim1]);

        /* on the first iteration and if mode is 1, scale according
           to the norms of the columns of the initial jacobian. */
        if (iter == 1) {
            if (*mode != 2) {
                for (int j = 0; j < *n; ++j) {
                    diag[j] = acnorm[j];
                    if (acnorm[j] == 0.0)
                        diag[j] = 1.0;
                }
            }

            /* On the first iteration, calculate the norm of the scaled x
               and initialize the step bound delta. */
            for (int j = 0; j < *n; ++j)
                wa3[j] = diag[j] * x[j];
            xnorm = enorm_(n, wa3);
            delta = *factor * xnorm;
            if (delta == 0)
                delta = *factor;
        }

        /* prepare QR factorization */
        double *tau = wa1;

        /* set all columns free */
        for (int j = 0; j < *n; j++)
            ipvt[j] = 0;
    
        /* query optimal size of work */
        int lapack_info = 0;
        int lwork = -1;
        dgeqp3_(m, n, fjac, ldfjac, ipvt, tau, tau, &lwork, &lapack_info);    // LAPACK
        assert(lapack_info == 0);
        lwork = tau[0];
        assert(lwork >= 3 * (*n) + 1);
    
        /* alloc work area. TODO: move to start of function */
        double *work = malloc(lwork * sizeof(*work));
        assert(work != NULL);
    
        /* compute the QR factorization of the jacobian. */
        dgeqp3_(m, n, fjac, ldfjac, ipvt, tau, work, &lwork, &lapack_info);   // LAPACK
        assert(lapack_info == 0);

        /* qtf <-- (Q transpose)*fvec */
        for (int i = 0; i < *m; ++i)
            wa4[i] = fvec[i];
        dormqr_("Left", "Transpose", m, &c1, n, fjac, ldfjac, tau, wa4, m, work, &lwork, &lapack_info);
        assert(lapack_info == 0);
        for (int j = 0; j < *n; ++j)
            qtf[j] = wa4[j];

        /* TODO: move to fini */
        free(work);

        /* Compute the norm of the scaled gradient. */
        gnorm = 0;
        if (fnorm != 0) {
            for (int j = 0; j < *n; ++j) {
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
        if (gnorm <= *gtol)
            *info = 4;
        if (*info != 0)
            break;

        /* rescale if necessary. */
        if (*mode != 2)
            for (int j = 0; j < *n; ++j)
                diag[j] = fmax(diag[j], wa2[j]);

        /* inner loop. */
        for (;;) {
            /* determine the levenberg-marquardt parameter. */
            lmpar_(n, fjac, ldfjac, ipvt, diag, qtf, &delta, &par, wa1, wa2, wa3, wa4);

            /* store the direction p and x + p. calculate the norm of p. */
            for (int j = 0; j < *n; ++j) {
                wa1[j] = -wa1[j];
                wa2[j] = x[j] + wa1[j];
                wa3[j] = diag[j] * wa1[j];
            }
            double pnorm = enorm_(n, wa3);

            /* On the first iteration, adjust the initial step bound. */
            if (iter == 1)
                delta = fmin(delta, pnorm);

            /* Evaluate the function at x + p and calculate its norm. */
            iflag = 1;
            if (fcn_der)
                (*fcn_der)(m, n, wa2, wa4, fjac, ldfjac, &iflag);
            if (fcn_dif)
                (*fcn_dif)(m, n, wa2, wa4, &iflag);
            *nfev += 1;
            if (iflag < 0)
                goto fini;
            double fnorm1 = enorm_(m, wa4);

            /* compute the scaled actual reduction. */
            double actred = -1;
            if (0.1 * fnorm1 < fnorm) {
                /* Computing 2nd power */
                double tmp = fnorm1 / fnorm;
                actred = 1.0 - tmp * tmp;
            }

            /* compute the scaled predicted reduction and
               the scaled directional derivative. */
            for (int j = 0; j < *n; ++j) {
                int l = ipvt[j] - 1;
                wa3[j] =  wa1[l];
            }
            dtrmv_("Upper", "Notrans", "non-unit", n, fjac, ldfjac, wa3, &c1);
            double temp1 = enorm_(n, wa3) / fnorm;
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
                for (int j = 0; j < *n; ++j) {
                    x[j] = wa2[j];
                    wa2[j] = diag[j] * x[j];
                }
                for (int i = 0; i < *m; ++i)
                    fvec[i] = wa4[i];
                xnorm = enorm_(n, wa2);
                fnorm = fnorm1;
                ++iter;
            }

            /* tests for convergence */
            if (fabs(actred) <= *ftol && prered <= *ftol && 0.5 * ratio <= 1)
                *info = 1;
            if (delta <= *xtol * xnorm)
                *info = 2;
            if (fabs(actred) <= *ftol && prered <= *ftol && 0.5 * ratio <= 1 && *info == 2)
                *info = 3;
            if (*info != 0)
                goto fini;

            /* tests for termination and stringent tolerances. */
            if (*nfev >= *maxfev)
                *info = 5;
            if (fabs(actred) <= epsmch && prered <= epsmch && 0.5 * ratio <= 1)
                *info = 6;
            if (delta <= epsmch * xnorm)
                *info = 7;
            if (gnorm <= epsmch)
                *info = 8;
            if (*info != 0)
                goto fini;

            /* repeat if iteration unsuccessful. */
            if (ratio >= 0.0001)
                break;
        } /* inner loop */
    } /* outer loop. */
 
 fini:
    /* termination, either normal or user imposed. */
    if (iflag < 0)
        *info = iflag;
    iflag = 0;
    if (*nprint > 0) {
        if (fcn_der)
            (*fcn_der)(m, n, x, fvec, fjac, ldfjac, &iflag);
        if (fcn_dif)
            (*fcn_dif)(m, n, x, fvec, &iflag);
    }
}