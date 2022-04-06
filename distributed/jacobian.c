#include <stddef.h>
#include <math.h>

#include "pminpack.h"

 /*
  *     this subroutine computes a forward-difference approximation 
  *     to the m by n jacobian matrix associated with a specified 
  *     problem of m functions in n variables. 
  *
  *       fcn is the name of the user-supplied subroutine which 
  *         calculates the functions. fcn must be declared 
  *         in an external statement in the user calling 
  *         program, and should be written as follows. 
  *
  *         subroutine fcn(m,n,x,fvec,iflag) 
  *         int m,n,iflag 
  *         double precision x(n),fvec(m) 
  *         ---------- 
  *         calculate the functions at x and 
  *         return this vector in fvec. 
  *         ---------- 
  *         return 
  *         end 
  *
  *         the value of iflag should not be changed by fcn unless 
  *         the user wants to terminate execution of fdjac2. 
  *         in this case set iflag to a negative int. 
  *
  *       m is a positive int input variable set to the number 
  *         of functions. 
  *
  *       n is a positive int input variable set to the number 
  *         of variables. n must not exceed m. 
  *
  *       x is an input array of length n. 
  *
  *       fvec is an input array of length m which must contain the 
  *         functions evaluated at x. 
  *
  *       fjac is an output m by n array which contains the 
  *         approximation to the jacobian matrix evaluated at x. 
  *
  *       ldfjac is a positive int input variable not less than m 
  *         which specifies the leading dimension of the array fjac. 
  *
  *       iflag is an int variable which can be used to terminate 
  *         the execution of fdjac2. see description of fcn. 
  *
  *       wa is a work array of length m. 
  */

void forward_difference_jacobian(pminpack_func_mn fcn, void *farg, int m, int n, double *x,
	    double const *fvec, double *fjac, int ldfjac, double *wa, int ictx)
{
	ptrdiff_t fjac_dim1 = ldfjac;
	double eps = sqrt(MINPACK_EPSILON);

	for (int j = 0; j < n; ++j) {
		double temp = x[j];
		double h = eps * fabs(x[j]);
		if (h == 0)
			h = eps;
		x[j] += h;
		(*fcn)(farg, m, n, x, wa);
		x[j] = temp;              // restore x[j]
		for (int i = 0; i < m; ++i)
			fjac[i + j * fjac_dim1] = (wa[i] - fvec[i]) / h;
	}
}
