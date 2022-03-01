#include <float.h>
#include <math.h>
#include <stdlib.h>

#include "minpack.h"

/*
 *     Subroutine lmder 
 *
 *     The purpose of lmder is to minimize the sum of the squares of 
 *     m nonlinear functions in n variables by a modification of 
 *     the Levenberg-Marquardt algorithm.  The user must provide a 
 *     subroutine which calculates the functions and the jacobian. 
 *
 *     The subroutine statement is 
 *
 *       subroutine lmder(fcn,m,n,x,fvec,fjac,ldfjac,ftol,xtol,gtol, 
 *                        maxfev,diag,mode,factor,nprint,info,nfev, 
 *                        njev,ipvt,qtf,wa1,wa2,wa3,wa4) 
 *
 *     where 
 *
 *       fcn is the name of the user-supplied subroutine which 
 *         calculates the functions and the jacobian.  fcn must 
 *         be declared in an external statement in the user 
 *         calling program, and should be written as follows. 
 *
 *         subroutine fcn(m,n,x,fvec,fjac,ldfjac,iflag) 
 *         integer m,n,ldfjac,iflag 
 *         double precision x(n),fvec(m),fjac(ldfjac,n) 
 *         ---------- 
 *         If iflag = 1 calculate the functions at x and 
 *         return this vector in fvec.  Do not alter fjac. 
 *         If iflag = 2 calculate the jacobian at x and 
 *         return this matrix in fjac.  Do not alter fvec. 
 *         ---------- 
 *         return 
 *         end 
 *
 *         The value of iflag should not be changed by fcn unless 
 *         the user wants to terminate execution of lmder. 
 *         In this case set iflag to a negative integer. 
 *
 *       m is a positive integer input variable set to the number 
 *         of functions. 
 *
 *       n is a positive integer input variable set to the number 
 *         of variables. n must not exceed m. 
 *
 *       x is an array of length n.  On input x must contain 
 *         an initial estimate of the solution vector.  On output x 
 *         contains the final estimate of the solution vector. 
 *
 *       fvec is an output array of length m which contains 
 *         the functions evaluated at the output x. 
 *
 *       fjac is an output m by n array.  The upper n by n submatrix 
 *         of fjac contains an upper triangular matrix R with 
 *         diagonal elements of nonincreasing magnitude such that 
 *
 *                t     t           t 
 *               P *(jac *jac)*P = R *R, 
 *
 *         where P is a permutation matrix and jac is the final 
 *         calculated jacobian.  Column j of P is column ipvt(j) 
 *         (see below) of the identity matrix.  The lower trapezoidal 
 *         part of fjac contains information generated during 
 *         the computation of R. 
 *
 *       ldfjac is a positive integer input variable not less than m 
 *         which specifies the leading dimension of the array fjac. 
 *
 *       ftol is a nonnegative input variable.  Termination 
 *         occurs when both the actual and predicted relative 
 *         reductions in the sum of squares are at most ftol. 
 *         Therefore, ftol measures the relative error desired 
 *         in the sum of squares. 
 *
 *       xtol is a nonnegative input variable. termination 
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
 *         initial step bound. this bound is set to the product of 
 *         factor and the euclidean norm of diag*x if nonzero, or else 
 *         to factor itself. in most cases factor should lie in the 
 *         interval (.1,100.).100. is a generally recommended value. 
 *
 *       nprint is an integer input variable that enables controlled 
 *         printing of iterates if it is positive.  In this case, 
 *         fcn is called with iflag = 0 at the beginning of the first 
 *         iteration and every nprint iterations thereafter and 
 *         immediately prior to return, with x, fvec, and fjac 
 *         available for printing.  fvec and fjac should not be 
 *         altered.  If nprint is not positive, no special calls 
 *         of fcn with iflag = 0 are made. 
 *
 *       info is an integer output variable.  If the user has 
 *         terminated execution, info is set to the (negative) 
 *         value of iflag. see description of fcn. otherwise, 
 *         info is set as follows. 
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
 *         info = 5  number of calls to fcn with iflag = 1 has 
 *                   reached maxfev. 
 *
 *         info = 6  ftol is too small.  No further reduction in 
 *                   the sum of squares is possible. 
 *
 *         info = 7  xtol is too small.  No further improvement in 
 *                   the approximate solution x is possible. 
 *
 *         info = 8  gtol is too small.  fvec is orthogonal to the 
 *                   columns of the jacobian to machine precision. 
 *
 *       nfev is an integer output variable set to the number of 
 *         calls to fcn with iflag = 1. 
 *
 *       njev is an integer output variable set to the number of 
 *         calls to fcn with iflag = 2. 
 *
 *       ipvt is an integer output array of length n.  ipvt 
 *         defines a permutation matrix p such that jac*P = Q*R, 
 *         where jac is the final calculated jacobian, Q is 
 *         orthogonal (not stored), and R is upper triangular 
 *         with diagonal elements of nonincreasing magnitude. 
 *         column j of P is column ipvt(j) of the identity matrix. 
 *
 *       qtf is an output array of length n which contains 
 *         the first n elements of the vector (Q transpose)*fvec. 
 *
 *       wa1, wa2, and wa3 are work arrays of length n. 
 *
 *       wa4 is a work array of length m. 
 *
 *     Argonne National Laboratory.  MINPACK project.  March 1980. 
 *     Burton S. Garbow, Kenneth E. Hillstrom, Jorge J. Moré 
*/
void lmder_(minpack_func_mnj fcn, const int *m, const int *n, double *x, 
	double *fvec, double *fjac, const int *ldfjac, const double *ftol,
	const double *xtol, const double *gtol, const int *maxfev, 
    double *diag, const int *mode, const double *factor, const int *nprint, 
    int *info, int *nfev, int *njev, int *ipvt, double *qtf, 
	double *wa1, double *wa2, double *wa3, double *wa4)
{
    /* Parameter adjustments */
    --wa4;
    --fvec;
    --wa3;
    --wa2;
    --wa1;
    --qtf;
    --ipvt;
    --diag;
    --x;
    int fjac_dim1 = *ldfjac;
    int fjac_offset = 1 + fjac_dim1;
    fjac -= fjac_offset;

    /* epsmch is the machine precision. */
    double epsmch = MINPACK_EPSILON;
    int iflag = 0;
    *info = 0;
    *nfev = 0;

    /* check the input parameters for errors. */
    if (*n <= 0 || *m < *n || *ldfjac < *m || *ftol < 0 || *xtol < 0 || *gtol < 0 || *maxfev <= 0 || *factor <= 0)
        goto fini;
    if (*mode == 2) {
        for (int j = 1; j <= *n; ++j)
            if (diag[j] <= 0)
                goto fini;
    }

    /* evaluate the function at the starting point and calculate its norm */
    iflag = 1;
    (*fcn)(m, n, &x[1], &fvec[1], &fjac[fjac_offset], ldfjac, &iflag);
    *nfev = 1;
    if (iflag < 0)
        goto fini;
    double fnorm = enorm_(m, &fvec[1]);
    double delta, xnorm, gnorm;

    /* initialize levenberg-marquardt parameter and iteration counter */
    double par = 0;
    int iter = 1;

    /* outer loop */
    for(;;) {
        /* calculate the jacobian matrix. */
        iflag = 2;
        (*fcn)(m, n, &x[1], &fvec[1], &fjac[fjac_offset], ldfjac, &iflag);
        *njev += 1;
        if (iflag < 0)
            break;

        /* if requested, call fcn to enable printing of iterates. */
        if (*nprint > 0) {
            iflag = 0;
            if ((iter - 1) % *nprint == 0)
                (*fcn)(m, n, &x[1], &fvec[1], &fjac[fjac_offset], ldfjac, &iflag);
            if (iflag < 0)
                break;
        }

        double *acnorm = wa2;
        for (int j = 1; j <= *n; ++j)
            acnorm[j] = enorm_(m, &fjac[j * fjac_dim1 + 1]);

        /* on the first iteration and if mode is 1, scale according
           to the norms of the columns of the initial jacobian. */
        if (iter == 1) {
            if (*mode != 2) {
                for (int j = 1; j <= *n; ++j) {
                    diag[j] = acnorm[j];
                    if (acnorm[j] == 0.0)
                        diag[j] = 1.0;
                }
            }

            /* On the first iteration, calculate the norm of the scaled x
               and initialize the step bound delta. */
            for (int j = 1; j <= *n; ++j)
                wa3[j] = diag[j] * x[j];
            xnorm = enorm_(n, &wa3[1]);
            delta = *factor * xnorm;
            if (delta == 0)
                delta = *factor;
        }

        double *rdiag = wa1;
        int pivot = 1;
        qrfac_(m, n, &fjac[fjac_offset], ldfjac, &pivot, &ipvt[1], n, &rdiag[1], &acnorm[1], &wa3[1]);

        /* form (q transpose)*fvec and store the first n components in qtf. */
        for (int i = 1; i <= *m; ++i)
            wa4[i] = fvec[i];
        for (int j = 1; j <= *n; ++j) {
            if (fjac[j + j * fjac_dim1] != 0) {
                double sum = 0;
                for (int i = j; i <= *m; ++i)
                    sum += fjac[i + j * fjac_dim1] * wa4[i];
                double temp = -sum / fjac[j + j * fjac_dim1];
                for (int i = j; i <= *m; ++i)
                    wa4[i] += fjac[i + j * fjac_dim1] * temp;
            }
            fjac[j + j * fjac_dim1] = wa1[j];
            qtf[j] = wa4[j];
        }

        /* Compute the norm of the scaled gradient. */
        gnorm = 0;
        if (fnorm != 0) {
            for (int j = 1; j <= *n; ++j) {
                int l = ipvt[j];
                if (acnorm[l] == 0)
                    continue;
                double sum = 0;
                for (int i = 1; i <= j; ++i)
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
            for (int j = 1; j <= *n; ++j)
                diag[j] = fmax(diag[j], wa2[j]);

        /* inner loop. */
        for (;;) {
            /* determine the levenberg-marquardt parameter. */
            lmpar_(n, &fjac[fjac_offset], ldfjac, &ipvt[1], &diag[1], &qtf[1], &delta, &par, &wa1[1], &wa2[1], &wa3[1], &wa4[1]);

            /* store the direction p and x + p. calculate the norm of p. */
            for (int j = 1; j <= *n; ++j) {
                wa1[j] = -wa1[j];
                wa2[j] = x[j] + wa1[j];
                wa3[j] = diag[j] * wa1[j];
            }
            double pnorm = enorm_(n, &wa3[1]);

            /* On the first iteration, adjust the initial step bound. */
            if (iter == 1)
                delta = fmin(delta, pnorm);

            /* Evaluate the function at x + p and calculate its norm. */
            iflag = 1;
            (*fcn)(m, n, &wa2[1], &wa4[1], &fjac[fjac_offset], ldfjac, &iflag);
            ++(*nfev);
            if (iflag < 0)
                goto fini;
            double fnorm1 = enorm_(m, &wa4[1]);

            /* compute the scaled actual reduction. */
            double actred = -1;
            if (0.1 * fnorm1 < fnorm) {
                /* Computing 2nd power */
                double tmp = fnorm1 / fnorm;
                actred = 1.0 - tmp * tmp;
            }

            /* compute the scaled predicted reduction and
               the scaled directional derivative. */
            for (int j = 1; j <= *n; ++j) {
                wa3[j] = 0.0;
                int l = ipvt[j];
                double temp = wa1[l];
                for (int i = 1; i <= j; ++i)
                    wa3[i] += fjac[i + j * fjac_dim1] * temp;
            }

            double temp1 = enorm_(n, &wa3[1]) / fnorm;
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
                for (int j = 1; j <= *n; ++j) {
                    x[j] = wa2[j];
                    wa2[j] = diag[j] * x[j];
                }
                for (int i = 1; i <= *m; ++i)
                    fvec[i] = wa4[i];
                xnorm = enorm_(n, &wa2[1]);
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
    if (*nprint > 0)
        (*fcn)(m, n, &x[1], &fvec[1], &fjac[fjac_offset], ldfjac, &iflag);
}
