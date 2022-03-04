#include <float.h>
#include <math.h>

#include "minpack.h"
/*
 *     subroutine r1updt 
 *
 *     given an m by n lower trapezoidal matrix S, an m-vector u, 
 *     and an n-vector v, the problem is to determine an 
 *     orthogonal matrix Q such that 
 *
 *                   t 
 *           (S + u*v )*Q
 *
 *     is again lower trapezoidal. 
 *
 *     this subroutine determines Q as the product of 2*(n - 1) 
 *     transformations 
 *
 *           gv(n-1)*...*gv(1)*gw(1)*...*gw(n-1) 
 *
 *     where gv(i), gw(i) are givens rotations in the (i,n) plane 
 *     which eliminate elements in the i-th and n-th planes, 
 *     respectively.  Q itself is not accumulated, rather the 
 *     information to recover the gv, gw rotations is returned. 
 *
 *     the subroutine statement is 
 *
 *       subroutine r1updt(m,n,s,ls,u,v,w,sing) 
 *
 *     where 
 *
 *       m is a positive integer input variable set to the number 
 *         of rows of S. 
 *
 *       n is a positive integer input variable set to the number 
 *         of columns of S.  n must not exceed m. 
 *
 *       S is an array of length ls.  On input S must contain the lower 
 *         trapezoidal matrix S stored by columns.  On output S contains 
 *         the lower trapezoidal matrix produced as described above. 
 *
 *       ls is a positive integer input variable not less than 
 *         (n*(2*m-n+1))/2. 
 *
 *       u is an input array of length m which must contain the 
 *         vector u. 
 *
 *       v is an array of length n.  On input v must contain the vector 
 *         v.  On output v(i) contains the information necessary to 
 *         recover the givens rotation gv(i) described above. 
 *
 *       w is an output array of length m.  w(i) contains information 
 *         necessary to recover the givens rotation gw(i) described 
 *         above. 
 *
 *       sing is a logical output variable.  sing is set true if any 
 *         of the diagonal elements of the output s are zero.  Otherwise 
 *         sing is set false. 
 *
 *     Argonne National Laboratory.  MINPACK project.  March 1980. 
 *     Burton S. Garbow, Kenneth E. Hillstrom, Jorge J. Moré, John L. Nazareth 
 */
void r1updt_(const int *m, const int *n, double *s, const int *ls, const double *u, double *v,
	     double *w, int *sing)
{
	(void)ls;

	/* Parameter adjustments */
	--w;
	--u;
	--v;
	--s;

	/* giant is the largest magnitude. */
	double giant = MINPACK_GIANT;

	/* initialize the diagonal element pointer. */
	int jj = *n * ((*m * 2) - *n + 1) / 2 - (*m - *n);

	/* move the nontrivial part of the last column of s into w. */
	int l = jj;
	for (int i = *n; i <= *m; ++i) {
		w[i] = s[l];
		++l;
	}

	/* rotate the vector v into a multiple of the n-th unit vector */
	/* in such a way that a spike is introduced into w. */
	for (int nmj = 1; nmj <= *n - 1; ++nmj) {
		int j = *n - nmj;
		jj -= *m - j + 1;
		w[j] = 0;
		if (v[j] == 0)
			continue;

		/* determine a givens rotation which eliminates the */
		/* j-th element of v. */
		double sin_, cos_, tau;
		if (fabs(v[*n]) < fabs(v[j])) {
			double cotan_ = v[*n] / v[j];
			sin_ = 0.5 / sqrt(0.25 + 0.25 * (cotan_ * cotan_));
			cos_ = sin_ * cotan_;
			tau = 1;
			if (fabs(cos_) * giant > 1)
				tau = 1 / cos_;
		} else {
			double tan_ = v[j] / v[*n];
			cos_ = 0.5 / sqrt(0.25 + 0.25 * (tan_ * tan_));
			sin_ = cos_ * tan_;
			tau = sin_;
		}

		/* apply the transformation to v and store the information */
		/* necessary to recover the givens rotation. */
		v[*n] = sin_ * v[j] + cos_ * v[*n];
		v[j] = tau;

		/* apply the transformation to s and extend the spike in w. */
		l = jj;
		for (int i = j; i <= *m; ++i) {
			double temp = cos_ * s[l] - sin_ * w[i];
			w[i] = sin_ * s[l] + cos_ * w[i];
			s[l] = temp;
			++l;
		}
	}

	/* add the spike from the rank 1 update to w. */
	for (int i = 1; i <= *m; ++i)
		w[i] += v[*n] * u[i];

	/* eliminate the spike. */
	*sing = 0;
	for (int j = 1; j <= *n - 1; ++j) {
		if (w[j] != 0) {
			/* determine a givens rotation which eliminates the */
			/* j-th element of the spike. */
			double cos_, sin_, tau;
			if (fabs(s[jj]) < fabs(w[j])) {
				double cotan_ = s[jj] / w[j];
				sin_ = 0.5 / sqrt(0.25 + 0.25 * (cotan_ * cotan_));
				cos_ = sin_ * cotan_;
				tau = 1;
				if (fabs(cos_) * giant > 1)
					tau = 1 / cos_;
			} else {
				double tan_ = w[j] / s[jj];
				cos_ = 0.5 / sqrt(0.25 + 0.25 * (tan_ * tan_));
				sin_ = cos_ * tan_;
				tau = sin_;
			}

			/* apply the transformation to s and reduce the spike in w. */
			l = jj;
			for (int i = j; i <= *m; ++i) {
				double temp = cos_ * s[l] + sin_ * w[i];
				w[i] = -sin_ * s[l] + cos_ * w[i];
				s[l] = temp;
				++l;
			}

			/* store the information necessary to recover the */
			/* givens rotation. */
			w[j] = tau;
		}

		/* test for zero diagonal elements in the output s. */
		if (s[jj] == 0)
			*sing = 1;
		jj += *m - j + 1;
	}

	/* move w back into the last column of the output s. */
	l = jj;
	for (int i = *n; i <= *m; ++i) {
		s[l] = w[i];
		++l;
	}
	if (s[jj] == 0)
		*sing = 1;
}
