#include "minpack.h"

/*
 *     Subroutine hybrd1
 *
 *     Rhe purpose of hybrd1 is to find a zero of a system of 
 *     n nonlinear functions in n variables by a modification 
 *     of the Powell hybrid method.  This is done by using the 
 *     more general nonlinear equation solver hybrd.  The user 
 *     must provide a subroutine which calculates the functions. 
 *     The jacobian is then calculated by a forward-difference 
 *     approximation. 
 *
 *     The subroutine statement is 
 *
 *       subroutine hybrd1(fcn,n,x,fvec,tol,info,wa,lwa) 
 *
 *     where 
 *
 *       fcn is the name of the user-supplied subroutine which 
 *         calculates the functions.  fcn must be declared 
 *         in an external statement in the user calling 
 *         program, and should be written as follows. 
 *
 *         subroutine fcn(n,x,fvec,iflag) 
 *         integer n,iflag 
 *         double precision x(n),fvec(n) 
 *         ---------- 
 *         calculate the functions at x and 
 *         return this vector in fvec. 
 *         --------- 
 *         return 
 *         end 
 *
 *         The value of iflag should not be changed by fcn unless 
 *         the user wants to terminate execution of hybrd1. 
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
 *       tol is a nonnegative input variable.  Termination occurs 
 *         when the algorithm estimates that the relative error 
 *         between x and the solution is at most tol. 
 *
 *       info is an integer output variable.  If the user has 
 *         terminated execution, info is set to the (negative) 
 *         value of iflag.  See description of fcn.  Otherwise, 
 *         info is set as follows. 
 *
 *         info = 0   improper input parameters. 
 *
 *         info = 1   algorithm estimates that the relative error 
 *                    between x and the solution is at most tol. 
 *
 *         info = 2   number of calls to fcn has reached or exceeded 
 *                    200*(n+1). 
 *
 *         info = 3   tol is too small.  No further improvement in 
 *                    the approximate solution x is possible. 
 *
 *         info = 4   iteration is not making good progress. 
 *
 *       wa is a work array of length lwa. 
 *
 *       lwa is a positive integer input variable not less than 
 *         (n*(3*n+13))/2. 
 *
 *     argonne national laboratory. minpack project. march 1980. 
 *     burton s. garbow, kenneth e. hillstrom, jorge j. more 
 */
void hybrd1_(minpack_func_n fcn, const int *n, double *x, double *fvec, 
    const double *tol, int *info, double *wa, const int *lwa)
{
    *info = 0;

    /* check the input parameters for errors. */
    if (*n <= 0 || *tol < 0 || *lwa < *n * (*n * 3 + 13) / 2)
	   return;

    /* call hybrd. */
    double factor = 100;
    int nfev = 0;
    int maxfev = (*n + 1) * 200;
    double xtol = *tol;
    int ml = *n - 1;
    int mu = *n - 1;
    double epsfcn = 0;
    int mode = 2;
    for (int j = 0; j < *n; ++j)
	   wa[j] = 1;
    int nprint = 0;
    int lr = *n * (*n + 1) / 2;
    int index = *n * 6 + lr;

    hybrd_(fcn, n, x, fvec, &xtol, &maxfev, &ml, &mu, &epsfcn, wa, 
        &mode, &factor, &nprint, info, &nfev, &wa[index], n, 
        &wa[*n * 6], &lr, &wa[*n], &wa[*n * 2], &wa[*n * 3], &wa[*n * 4], &wa[*n * 5]);
    if (*info == 5)
	   *info = 4;
}