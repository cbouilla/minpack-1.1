#include "pminpack.h"

/*     Subroutine plmdif1
 *
 *     The purpose of plmdif1 is to minimize the sum of the squares of
 *     m nonlinear functions in n variables by a modification of the
 *     levenberg-marquardt algorithm.  This is done by using the more
 *     general least-squares solver plmdif.  The user must provide a
 *     subroutine which calculates the functions.  The jacobian is
 *     then calculated by a forward-difference approximation.
 *
 *       fcn is the name of the user-supplied subroutine which
 *         calculates the functions.  fcn must be declared
 *         in an external statement in the user calling
 *         program, and should be written as follows.
 *
 *         subroutine fcn(m,n,x,fvec)
 *         integer m,n
 *         double precision x(n),fvec(m)
 *         ----------
 *         calculate the functions at x and
 *         return this vector in fvec.
 *         ----------
 *         return
 *         end
 *
 *       m is a positive integer input variable set to the number
 *         of functions.
 *
 *       n is a positive integer input variable set to the number
 *         of variables.  n must not exceed m.
 *
 *       x is an array of length n.  On input x must contain
 *         an initial estimate of the solution vector.  On output x
 *         contains the final estimate of the solution vector.
 *
 *       fvec is an output array of length m which contains
 *         the functions evaluated at the output x.
 *
 *       tol is a nonnegative input variable.  Termination occurs
 *         when the algorithm estimates either that the relative
 *         error in the sum of squares is at most tol or that
 *         the relative error between x and the solution is at
 *         most tol.
 *
 *       info is an integer output variable set as follows.
 *
 *         info = 0  improper input parameters.
 *
 *         info = 1  algorithm estimates that the relative error
 *                   in the sum of squares is at most tol.
 *
 *         info = 2  algorithm estimates that the relative error
 *                   between x and the solution is at most tol.
 *
 *         info = 3  conditions for info = 1 and info = 2 both hold.
 *
 *         info = 4  fvec is orthogonal to the columns of the
 *                   jacobian to machine precision.
 *
 *         info = 5  number of calls to fcn has reached or
 *                   exceeded 200*(n+1).
 *
 *         info = 6  tol is too small. no further reduction in
 *                   the sum of squares is possible.
 *
 *         info = 7  tol is too small. no further improvement in
 *                   the approximate solution x is possible.
 *
 *       ctx is a BLACS grid context.
 */

int plmdif1(pminpack_func_mn fcn, void *farg, int m, int n, double *x, double *fvec, double tol, int ctx)
{
	int info = 0;

	/* check the input parameters for errors. */
	if (n <= 0 || m < n || tol < 0)
		return info;

	/* call plmdif */
	int maxfev = (n + 1) * 200;
	double ftol = tol;
	double xtol = tol;
	double gtol = 0;
    int nfev = 0;
	info = plmdif(fcn, farg, m, n, x, fvec, ftol, xtol, gtol, maxfev, &nfev, ctx);
	if (info == 8)
		info = 4;
    return info;
}