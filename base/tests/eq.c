#include "f2c.h"
#include <math.h>

#include "eq.h"

double d_sign(const double *a, const double *b)
{
	double x = fabs(*a);
	return *b >= 0 ? x : -x;
}

char *problem_name[] = {
	"Rosenbrock function",
	"Powell singular function",
	"Powell badly scaled function",
	"Wood function",
	"Helical valley function",
	"Watson function",
	"Chebyquad function",
	"Brown almost-linear function",
	"Discrete boundary value function",
	"Discrete integral equation function",
	"Trigonometric function",
	"Variably dimensioned function",
	"Broyden tridiagonal function",
	"Broyden banded function",
};

struct test_case tests[55] = {
	{1, 2, 1,     22, 1, 0           },  // 01  
	{1, 2, 10,     9, 1, 0           },  // 02
	{1, 2, 100,    9, 1, .5373479E-12},  // 03
	{2, 4, 1,    123, 4, .1736920E-34},  // 04
	{2, 4, 10,   110, 4, .1828554E-32},  // 05
	{2, 4, 100,  156, 4, .2539466E-34},  // 06
	{3, 2, 1,    181, 1, .1712494E-08},  // 07
	{3, 2, 10,    11, 1, .3744485E-07},  // 08
	{4, 4, 1,     94, 1, .4005067E-10},  // 09
	{4, 4, 10,   234, 1, .5154548E-09},  // 10
	{4, 4, 100,  514, 1, .2491133E-10},  // 11
	{5, 3, 1,     27, 1, .2753457E-12},  // 12
	{5, 3, 10,    31, 1, .4427907E-13},  // 13
	{5, 3, 100,   40, 1, .1610836E-12},  // 14
	{6, 6, 1,     96, 1, .3191646E-12},  // 15
	{6, 6, 10,   310, 1, .2336128E-10},  // 16
	{6, 9, 1,    167, 1, .1817744E-13},  // 17
	{6, 9, 10,   167, 1, .1138869E-11},  // 18
	{7, 5, 1,     17, 1, .3540292E-11},  // 19
	{7, 5, 10,   256, 1, .1902384E-09},  // 20
	{7, 5, 100,  522, 1, .2269130E-10},  // 21
	{7, 6, 1,     25, 1, .6841874E-09},  // 22
	{7, 6, 10,   172, 1, .2281996E-10},  // 23
	{7, 6, 100,  274, 1, .1260220E-10},  // 24
	{7, 7, 1,     20, 1, .2398700E-08},  // 25
	{7, 7, 10,   525, 1, .1624356E-09},  // 26
	{7, 7, 100,  657, 4, .2416003E+00},  // 27
	{7, 8, 1,    120, 4, .6440508E-01},  // 28
	{7, 9, 1,     41, 1, .1889808E-08},  // 29
	{8, 10, 1,    31, 1, .1312131E-13},  // 30
	{8, 10, 10,   31, 1, .1601186E-13},  // 31
	{8, 10, 100,  38, 1, .4213000E-14},  // 32
	{8, 30, 1,   113, 1, .2073414E-12},  // 33
	{8, 40, 1,   196, 1, .8885114E-13},  // 34
	{9, 10, 1,    16, 1, .2497793E-14},  // 35
	{9, 10, 10,   19, 1, .1725954E-12},  // 36
	{9, 10, 100,  52, 1, .4174414E-09},  // 37
	{10, 1, 1,     7, 1, .5551115E-16},  // 38
	{10, 1, 10,    9, 1, .5551115E-16},  // 39
	{10, 1, 100,  16, 1, .2775558E-16},  // 40
	{10, 10, 1,   16, 1, .5009132E-14},  // 41
	{10, 10, 10,  19, 1, .2188310E-12},  // 42
	{10, 10, 100, 39, 1, .3034804E-14},  // 43
	{11, 10, 1,  130, 4, .5296377E-02},  // 44
	{11, 10, 10,  84, 1, .5917302E-10},  // 45
	{11, 10, 100, 85, 1, .1856423E-08},  // 46
	{12, 10, 1,   31, 1, .5119467E-11},  // 47
	{12, 10, 10,  35, 1, .7399640E-10},  // 48
	{12, 10, 100, 66, 1, 0           },  // 49
	{13, 10, 1,   21, 1, .1493879E-07},  // 50
	{13, 10, 10,  59, 1, .5337442E-08},  // 51
	{13, 10, 100, 42, 1, .9878880E-10},  // 52
	{14, 10, 1,   30, 1, .2057898E-08},  // 53
	{14, 10, 10,  45, 1, .7953613E-08},  // 54
	{14, 10, 100, 58, 1, .4526424E-09},  // 55
};

/*
 *     Subroutine vecfcn
 *
 *     This subroutine defines fourteen test functions.  The first 
 *     five test functions are of dimensions 2,4,2,4,3, respectively, 
 *     while the remaining test functions are of variable dimension 
 *     n for any n greater than or equal to 1 (problem 6 is an 
 *     exception to this, since it does not allow n = 1). 
 *
 *     The subroutine statement is 
 *
 *       subroutine vecfcn(n,x,fvec,nprob) 
 *
 *     where 
 *
 *       n is a positive int input variable. 
 *
 *       x is an input array of length n. 
 *
 *       fvec is an output array of length n which contains the nprob 
 *         function vector evaluated at x. 
 *
 *       nprob is a positive int input variable which defines the 
 *         number of the problem. nprob must not exceed 14. 
 */
void vecfcn(int n, double *x, double *fvec, int nprob)
{
	/* Local variables */

	/* Parameter adjustments */
	--fvec;
	--x;

	/* Function Body */

/*     PROBLEM SELECTOR. */

	switch (nprob) {
	case 1:
		/* Rosenbrock function. */
		fvec[1] = 10 * (x[2] - x[1] * x[1]);
		fvec[2] = 1 - x[1];
		break;

	case 2:
		/* Powell singular function. */
		fvec[1] = x[1] + 10 * x[2];
		fvec[2] = sqrt(5) * (x[3] - x[4]);
		double a = x[2] - 2 * x[3];
		fvec[3] = a * a;
		double b = x[1] - x[4];
		fvec[4] = sqrt(10) * (b * b);
		break;

	case 3:
		/* Powell badly scaled function. */
		fvec[1] = 1e4 * x[1] * x[2] - 1;
		fvec[2] = exp(-x[1]) + exp(-x[2]) - 1.0001;
		break;

	case 4:
		/* WOOD FUNCTION. */
		double c = x[2] - x[1] * x[1];
		double d = x[4] - x[3] * x[3];
		fvec[1] = -200 * x[1] * c - (1 - x[1]);
		fvec[2] = 200 * c + 20.2 * (x[2] - 1) + 19.8 * (x[4] - 1);
		fvec[3] = -180 * x[3] * d - (1 - x[3]);
		fvec[4] = 180 * d + 20.2 * (x[4] - 1) + 19.8 * (x[2] - 1);
		break;

	case 5:
		/* HELICAL VALLEY FUNCTION. */
		double twopi = 2 * M_PI;
		double tmp1 = x[2] < 0 ? -0.25 : 0.25;
		if (x[1] > 0)
			tmp1 = atan(x[2] / x[1]) / twopi;
		if (x[1] < 0)
			tmp1 = atan(x[2] / x[1]) / twopi + .5;
		double tmp2 = sqrt(x[1] * x[1] + x[2] * x[2]);
		fvec[1] = 10 * (x[3] - 10 * tmp1);
		fvec[2] = 10 * (tmp2 - 1);
		fvec[3] = x[3];
		break;

	case 6:
		/* WATSON FUNCTION. */
		for (int k = 1; k <= n; ++k)
			fvec[k] = 0;
		for (int i = 1; i <= 29; ++i) {
			double ti = i / 29.;
			double s1 = 0;
			double dx = 1;
			for (int j = 2; j <= n; ++j) {
				s1 += (j - 1) * dx * x[j];
				dx = ti * dx;
			}
			double s2 = 0;
			dx = 1;
			for (int j = 1; j <= n; ++j) {
				s2 += dx * x[j];
				dx = ti * dx;
			}
			double temp1 = s1 - s2 * s2 - 1;
			double temp2 = 2 * ti * s2;
			double temp = 1 / ti;
			for (int k = 1; k <= n; ++k) {
				fvec[k] += temp * (k - 1 - temp2) * temp1;
				temp = ti * temp;
			}
		}
		double temp = x[2] - x[1] * x[1] - 1;
		fvec[1] += x[1] * (1 - 2 * temp);
		fvec[2] += temp;
		break;

	case 7:
		/* chebyquad function. */
		for (int i = 1; i <= n; ++i)
			fvec[i] = 0;
		for (int j = 1; j <= n; ++j) {
			tmp1 = 1;
			tmp2 = 2 * x[j] - 1;
			temp = 2 * tmp2;
			for (int i = 1; i <= n; ++i) {
				fvec[i] += tmp2;
				double ti = temp * tmp2 - tmp1;
				tmp1 = tmp2;
				tmp2 = ti;
			}
		}
		double dx = 1. / n;
		double iev = -1;
		for (int i = 1; i <= n; ++i) {
			fvec[i] = dx * fvec[i];
			if (iev > 0)
				fvec[i] += 1. / (i * i - 1);
			iev = -iev;
		}
		break;

	case 8:
		/* brown almost-linear function. */
		double sum = -(n + 1);
		double prod = 1;
		for (int j = 1; j <= n; ++j) {
			sum += x[j];
			prod *= x[j];
		}
		for (int i = 1; i <= n; ++i)
			fvec[i] = x[i] + sum;
		fvec[n] = prod - 1;
		break;

	case 9:
		/* DISCRETE BOUNDARY VALUE FUNCTION. */
		double h = 1. / (n + 1);
		for (int k = 1; k <= n; ++k) {
			double a = x[k] + k * h + 1;
			double temp = a * (a * a);
			double temp1 = 0;
			if (k != 1)
				temp1 = x[k - 1];
			double temp2 = 0;
			if (k != n)
				temp2 = x[k + 1];
			fvec[k] = 2 * x[k] - temp1 - temp2 + temp * (h * h) / 2;
		}
		break;

	case 10:
		/*     DISCRETE INTEGRAL EQUATION FUNCTION. */
		h = 1. / (n + 1);
		for (int k = 1; k <= n; ++k) {
			double tk = k * h;
			double sum1 = 0;
			for (int j = 1; j <= k; ++j) {
				double tj = j * h;
				double a = x[j] + tj + 1;
				temp = a * (a * a);
				sum1 += tj * temp;
			}
			double sum2 = 0;
			for (int j = k + 1; j <= n; ++j) {
				double tj = j * h;
				double a = x[j] + tj + 1;
				temp = a * (a * a);
				sum2 += (1 - tj) * temp;
			}
			fvec[k] = x[k] + h * ((1 - tk) * sum1 + tk * sum2) / 2;
		}
		break;

	case 11:
		/*     TRIGONOMETRIC FUNCTION. */
		sum = 0;
		for (int j = 1; j <= n; ++j) {
			fvec[j] = cos(x[j]);
			sum += fvec[j];
		}
		for (int k = 1; k <= n; ++k) {
			fvec[k] = (n + k) - sin(x[k]) - sum - k * fvec[k];
		}
		break;

	case 12:
		/*     VARIABLY DIMENSIONED FUNCTION. */
		sum = 0;
		for (int j = 1; j <= n; ++j)
			sum += j * (x[j] - 1.);
		temp = sum * (1. + 2. * (sum * sum));
		for (int k = 1; k <= n; ++k)
			fvec[k] = x[k] - 1 + k * temp;
		break;

	case 13:
		/*     BROYDEN TRIDIAGONAL FUNCTION. */
		for (int k = 1; k <= n; ++k) {
			double temp = (3 - 2 * x[k]) * x[k];
			double temp1 = 0;
			if (k != 1)
				temp1 = x[k - 1];
			double temp2 = 0;
			if (k != n)
				temp2 = x[k + 1];
			fvec[k] = temp - temp1 - 2 * temp2 + 1;
		}
		break;

	case 14:
		/*     BROYDEN BANDED FUNCTION. */
		int ml = 5;
		int mu = 1;
		for (int k = 1; k <= n; ++k) {
			double k1 = fmax(1, k - ml);
			double k2 = fmin(k + mu, n);
			double temp = 0;
			for (int j = k1; j <= k2; ++j)
				if (j != k)
					temp += x[j] * (1 + x[j]);
			fvec[k] = x[k] * (2 + 5 * (x[k] * x[k])) + 1 - temp;
		}
	}
}

/*
 *     Subroutine initpt
 *
 *     This subroutine specifies the standard starting points for
 *     the functions defined by subroutine vecfcn.  The subroutine
 *     returns in x a multiple (factor) of the standard starting
 *     point.  For the sixth function the standard starting point is
 *     zero, so in this case, if factor is not unity, then the
 *     subroutine returns the vector  x(j) = factor, j=1,...,n.
 *
 *     The subroutine statement is
 *
 *       subroutine initpt(n,x,nprob,factor)
 *
 *     where
 *
 *       n is a positive int input variable.
 *
 *       x is an output array of length n which contains the standard
 *         starting point for problem nprob multiplied by factor.
 *
 *       nprob is a positive int input variable which defines the
 *         number of the problem.  nprob must not exceed 14.
 *
 *       factor is an input variable which specifies the multiple of
 *         the standard starting point.  If factor is unity, no
 *         multiplication is performed.
 */
void initpt(int n, double *x, int nprob, double factor)
{
	--x;

	/*     SELECTION OF INITIAL POINT. */
	switch (nprob) {
	case 1:
		/* ROSENBROCK FUNCTION. */
		x[1] = -1.2;
		x[2] = 1;
		break;
	case 2:
		/* POWELL SINGULAR FUNCTION. */
		x[1] = 3;
		x[2] = -1;
		x[3] = 0;
		x[4] = 1;
		break;
	case 3:
		/* POWELL BADLY SCALED FUNCTION. */
		x[1] = 0;
		x[2] = 1;
		break;
	case 4:
		/* WOOD FUNCTION. */
		x[1] = -3;
		x[2] = -1;
		x[3] = -3;
		x[4] = -1;
		break;
	case 5:
		/* HELICAL VALLEY FUNCTION. */
		x[1] = -1;
		x[2] = 0;
		x[3] = 0;
		break;
	case 6:
		/* WATSON FUNCTION. */
		for (int j = 1; j <= n; ++j)
			x[j] = 0;
		break;
	case 7:
		/* CHEBYQUAD FUNCTION. */
		double h = 1. / (n + 1);
		for (int j = 1; j <= n; ++j)
			x[j] = j * h;
		break;
	case 8:
		/* BROWN ALMOST-LINEAR FUNCTION. */
		for (int j = 1; j <= n; ++j)
			x[j] = 0.5;
		break;
	case 9:
	case 10:
		/* DISCRETE BOUNDARY VALUE AND INTEGRAL EQUATION FUNCTIONS. */
		h = 1. / (n + 1);
		for (int j = 1; j <= n; ++j) {
			double tj = j * h;
			x[j] = tj * (tj - 1);
		}
		break;
	case 11:
		/* TRIGONOMETRIC FUNCTION. */
		h = 1. / n;
		for (int j = 1; j <= n; ++j)
			x[j] = h;
		break;
	case 12:
		/* VARIABLY DIMENSIONED FUNCTION. */
		h = 1. / n;
		for (int j = 1; j <= n; ++j)
			x[j] = 1. - j * h;
		break;
	case 13:
	case 14:
		/* BROYDEN TRIDIAGONAL AND BANDED FUNCTIONS. */
		for (int j = 1; j <= n; ++j)
			x[j] = -1;
	}

	/* COMPUTE MULTIPLE OF INITIAL POINT. */
	if (factor == 1)
		return;
	if (nprob == 6) {
		for (int j = 1; j <= n; ++j)
			x[j] = factor;
	} else {
		for (int j = 1; j <= n; ++j)
			x[j] = factor * x[j];
	}
}

/*     Subroutine vecjac
 *
 *     This subroutine defines the jacobian matrices of fourteen
 *     test functions.  The problem dimensions are as described
 *     in the prologue comments of vecfcn.
 *
 *     The subroutine statement is
 *
 *       subroutine vecjac(n,x,fjac,ldfjac,nprob)
 *
 *     where
 *
 *       n is a positive int variable.
 *
 *       x is an array of length n.
 *
 *       fjac is an n by n array.  On output fjac contains the
 *         jacobian matrix of the nprob function evaluated at x.
 *
 *       ldfjac is a positive int variable not less than n 
 *         which specifies the leading dimension of the array fjac.
 *
 *       nprob is a positive int variable which defines the
 *         number of the problem.  nprob must not exceed 14.
 */
void vecjac(int n, double *x, double *fjac, int ldfjac, int nprob)
{
	/* Initialized data */

	static double zero = 0.;
	static double fiftn = 15.;
	static double twenty = 20.;
	static double hundrd = 100.;
	static double c1 = 1e4;
	static double c3 = 200.;
	static double c4 = 20.2;
	static double c5 = 19.8;
	static double c6 = 180.;
	static double c9 = 29.;
	static double one = 1.;
	static double two = 2.;
	static double three = 3.;
	static double four = 4.;
	static double five = 5.;
	static double six = 6.;
	static double eight = 8.;
	static double ten = 10.;

	/* System generated locals */
	int fjac_offset, i__1, i__2, i__3, i__4;
	double d__1, d__2;

	/* Local variables */
	static double h__;
	static int i__, j, k, k1, k2, ml;
	static double ti, tj, tk;
	static int mu;
	static double tpi, sum, sum1, sum2, prod, temp, temp1, temp2, temp3, temp4;

	/* Parameter adjustments */
	--x;
	ldfjac = ldfjac;
	fjac_offset = 1 + ldfjac;
	fjac -= fjac_offset;

	/* Function Body */

/*     JACOBIAN ROUTINE SELECTOR. */

	switch (nprob) {
	case 1:
		/*     ROSENBROCK FUNCTION. */
		fjac[ldfjac + 1] = -one;
		fjac[(ldfjac << 1) + 1] = zero;
		fjac[ldfjac + 2] = -twenty * x[1];
		fjac[(ldfjac << 1) + 2] = ten;
		break;
	case 2:
		goto L20;
	case 3:
		goto L50;
	case 4:
		goto L60;
	case 5:
		goto L90;
	case 6:
		goto L100;
	case 7:
		goto L200;
	case 8:
		goto L230;
	case 9:
		goto L290;
	case 10:
		goto L320;
	case 11:
		goto L350;
	case 12:
		goto L380;
	case 13:
		goto L420;
	case 14:
		goto L450;
	}


 L10:
	goto L490;

/*     POWELL SINGULAR FUNCTION. */

 L20:
	for (k = 1; k <= 4; ++k) {
		for (j = 1; j <= 4; ++j) {
			fjac[k + j * ldfjac] = zero;
/* L30: */
		}
/* L40: */
	}
	fjac[ldfjac + 1] = one;
	fjac[(ldfjac << 1) + 1] = ten;
	fjac[ldfjac * 3 + 2] = sqrt(five);
	fjac[(ldfjac << 2) + 2] = -fjac[ldfjac * 3 + 2];
	fjac[(ldfjac << 1) + 3] = two * (x[2] - two * x[3]);
	fjac[ldfjac * 3 + 3] = -two * fjac[(ldfjac << 1) + 3];
	fjac[ldfjac + 4] = two * sqrt(ten) * (x[1] - x[4]);
	fjac[(ldfjac << 2) + 4] = -fjac[ldfjac + 4];
	goto L490;

/*     POWELL BADLY SCALED FUNCTION. */

 L50:
	fjac[ldfjac + 1] = c1 * x[2];
	fjac[(ldfjac << 1) + 1] = c1 * x[1];
	fjac[ldfjac + 2] = -exp(-x[1]);
	fjac[(ldfjac << 1) + 2] = -exp(-x[2]);
	goto L490;

/*     WOOD FUNCTION. */

 L60:
	for (k = 1; k <= 4; ++k) {
		for (j = 1; j <= 4; ++j) {
			fjac[k + j * ldfjac] = zero;
/* L70: */
		}
/* L80: */
	}
/* Computing 2nd power */
	d__1 = x[1];
	temp1 = x[2] - three * (d__1 * d__1);
/* Computing 2nd power */
	d__1 = x[3];
	temp2 = x[4] - three * (d__1 * d__1);
	fjac[ldfjac + 1] = -c3 * temp1 + one;
	fjac[(ldfjac << 1) + 1] = -c3 * x[1];
	fjac[ldfjac + 2] = -two * c3 * x[1];
	fjac[(ldfjac << 1) + 2] = c3 + c4;
	fjac[(ldfjac << 2) + 2] = c5;
	fjac[ldfjac * 3 + 3] = -c6 * temp2 + one;
	fjac[(ldfjac << 2) + 3] = -c6 * x[3];
	fjac[(ldfjac << 1) + 4] = c5;
	fjac[ldfjac * 3 + 4] = -two * c6 * x[3];
	fjac[(ldfjac << 2) + 4] = c6 + c4;
	goto L490;

/*     HELICAL VALLEY FUNCTION. */

 L90:
	tpi = eight * atan(one);
/* Computing 2nd power */
	d__1 = x[1];
/* Computing 2nd power */
	d__2 = x[2];
	temp = d__1 * d__1 + d__2 * d__2;
	temp1 = tpi * temp;
	temp2 = sqrt(temp);
	fjac[ldfjac + 1] = hundrd * x[2] / temp1;
	fjac[(ldfjac << 1) + 1] = -hundrd * x[1] / temp1;
	fjac[ldfjac * 3 + 1] = ten;
	fjac[ldfjac + 2] = ten * x[1] / temp2;
	fjac[(ldfjac << 1) + 2] = ten * x[2] / temp2;
	fjac[ldfjac * 3 + 2] = zero;
	fjac[ldfjac + 3] = zero;
	fjac[(ldfjac << 1) + 3] = zero;
	fjac[ldfjac * 3 + 3] = one;
	goto L490;

/*     WATSON FUNCTION. */

 L100:
	i__1 = n;
	for (k = 1; k <= i__1; ++k) {
		i__2 = n;
		for (j = k; j <= i__2; ++j) {
			fjac[k + j * ldfjac] = zero;
/* L110: */
		}
/* L120: */
	}
	for (i__ = 1; i__ <= 29; ++i__) {
		ti = (double)i__ / c9;
		sum1 = zero;
		temp = one;
		i__1 = n;
		for (j = 2; j <= i__1; ++j) {
			i__2 = j - 1;
			sum1 += (double)i__2 *temp * x[j];
			temp = ti * temp;
/* L130: */
		}
		sum2 = zero;
		temp = one;
		i__1 = n;
		for (j = 1; j <= i__1; ++j) {
			sum2 += temp * x[j];
			temp = ti * temp;
/* L140: */
		}
/* Computing 2nd power */
		d__1 = sum2;
		temp1 = two * (sum1 - d__1 * d__1 - one);
		temp2 = two * sum2;
/* Computing 2nd power */
		d__1 = ti;
		temp = d__1 * d__1;
		tk = one;
		i__1 = n;
		for (k = 1; k <= i__1; ++k) {
			tj = tk;
			i__2 = n;
			for (j = k; j <= i__2; ++j) {
				i__3 = k - 1;
				i__4 = j - 1;
				fjac[k + j * ldfjac] += tj * (((double)i__3 / ti - temp2) * ((double)i__4 / ti - temp2) - temp1);
				tj = ti * tj;
/* L150: */
			}
			tk = temp * tk;
/* L160: */
		}
/* L170: */
	}
/* Computing 2nd power */
	d__1 = x[1];
	fjac[ldfjac + 1] = fjac[ldfjac + 1] + six * (d__1 * d__1) - two * x[2] + three;
	fjac[(ldfjac << 1) + 1] -= two * x[1];
	fjac[(ldfjac << 1) + 2] += one;
	i__1 = n;
	for (k = 1; k <= i__1; ++k) {
		i__2 = n;
		for (j = k; j <= i__2; ++j) {
			fjac[j + k * ldfjac] = fjac[k + j * ldfjac];
/* L180: */
		}
/* L190: */
	}
	goto L490;

/*     CHEBYQUAD FUNCTION. */

 L200:
	tk = one / (double)(n);
	i__1 = n;
	for (j = 1; j <= i__1; ++j) {
		temp1 = one;
		temp2 = two * x[j] - one;
		temp = two * temp2;
		temp3 = zero;
		temp4 = two;
		i__2 = n;
		for (k = 1; k <= i__2; ++k) {
			fjac[k + j * ldfjac] = tk * temp4;
			ti = four * temp2 + temp * temp4 - temp3;
			temp3 = temp4;
			temp4 = ti;
			ti = temp * temp2 - temp1;
			temp1 = temp2;
			temp2 = ti;
/* L210: */
		}
/* L220: */
	}
	goto L490;

/*     BROWN ALMOST-LINEAR FUNCTION. */

 L230:
	prod = one;
	i__1 = n;
	for (j = 1; j <= i__1; ++j) {
		prod = x[j] * prod;
		i__2 = n;
		for (k = 1; k <= i__2; ++k) {
			fjac[k + j * ldfjac] = one;
/* L240: */
		}
		fjac[j + j * ldfjac] = two;
/* L250: */
	}
	i__1 = n;
	for (j = 1; j <= i__1; ++j) {
		temp = x[j];
		if (temp != zero) {
			goto L270;
		}
		temp = one;
		prod = one;
		i__2 = n;
		for (k = 1; k <= i__2; ++k) {
			if (k != j) {
				prod = x[k] * prod;
			}
/* L260: */
		}
 L270:
		fjac[n + j * ldfjac] = prod / temp;
/* L280: */
	}
	goto L490;

/*     DISCRETE BOUNDARY VALUE FUNCTION. */

 L290:
	i__1 = n + 1;
	h__ = one / (double)i__1;
	i__1 = n;
	for (k = 1; k <= i__1; ++k) {
/* Computing 2nd power */
		d__1 = x[k] + (double)k *h__ + one;
		temp = three * (d__1 * d__1);
		i__2 = n;
		for (j = 1; j <= i__2; ++j) {
			fjac[k + j * ldfjac] = zero;
/* L300: */
		}
/* Computing 2nd power */
		d__1 = h__;
		fjac[k + k * ldfjac] = two + temp * (d__1 * d__1) / two;
		if (k != 1) {
			fjac[k + (k - 1) * ldfjac] = -one;
		}
		if (k != n) {
			fjac[k + (k + 1) * ldfjac] = -one;
		}
/* L310: */
	}
	goto L490;

/*     DISCRETE INTEGRAL EQUATION FUNCTION. */

 L320:
	i__1 = n + 1;
	h__ = one / (double)i__1;
	i__1 = n;
	for (k = 1; k <= i__1; ++k) {
		tk = (double)k *h__;
		i__2 = n;
		for (j = 1; j <= i__2; ++j) {
			tj = (double)j *h__;
/* Computing 2nd power */
			d__1 = x[j] + tj + one;
			temp = three * (d__1 * d__1);
/* Computing MIN */
			d__1 = tj * (one - tk), d__2 = tk * (one - tj);
			fjac[k + j * ldfjac] = h__ * min(d__1, d__2) * temp / two;
/* L330: */
		}
		fjac[k + k * ldfjac] += one;
/* L340: */
	}
	goto L490;

/*     TRIGONOMETRIC FUNCTION. */

 L350:
	i__1 = n;
	for (j = 1; j <= i__1; ++j) {
		temp = sin(x[j]);
		i__2 = n;
		for (k = 1; k <= i__2; ++k) {
			fjac[k + j * ldfjac] = temp;
/* L360: */
		}
		i__2 = j + 1;
		fjac[j + j * ldfjac] = (double)i__2 *temp - cos(x[j]);
/* L370: */
	}
	goto L490;

/*     VARIABLY DIMENSIONED FUNCTION. */

 L380:
	sum = zero;
	i__1 = n;
	for (j = 1; j <= i__1; ++j) {
		sum += (double)j *(x[j] - one);
/* L390: */
	}
/* Computing 2nd power */
	d__1 = sum;
	temp = one + six * (d__1 * d__1);
	i__1 = n;
	for (k = 1; k <= i__1; ++k) {
		i__2 = n;
		for (j = k; j <= i__2; ++j) {
			i__3 = k * j;
			fjac[k + j * ldfjac] = (double)i__3 *temp;
			fjac[j + k * ldfjac] = fjac[k + j * ldfjac];
/* L400: */
		}
		fjac[k + k * ldfjac] += one;
/* L410: */
	}
	goto L490;

/*     BROYDEN TRIDIAGONAL FUNCTION. */

 L420:
	i__1 = n;
	for (k = 1; k <= i__1; ++k) {
		i__2 = n;
		for (j = 1; j <= i__2; ++j) {
			fjac[k + j * ldfjac] = zero;
/* L430: */
		}
		fjac[k + k * ldfjac] = three - four * x[k];
		if (k != 1) {
			fjac[k + (k - 1) * ldfjac] = -one;
		}
		if (k != n) {
			fjac[k + (k + 1) * ldfjac] = -two;
		}
/* L440: */
	}
	goto L490;

/*     BROYDEN BANDED FUNCTION. */

 L450:
	ml = 5;
	mu = 1;
	i__1 = n;
	for (k = 1; k <= i__1; ++k) {
		i__2 = n;
		for (j = 1; j <= i__2; ++j) {
			fjac[k + j * ldfjac] = zero;
/* L460: */
		}
/* Computing MAX */
		i__2 = 1, i__3 = k - ml;
		k1 = max(i__2, i__3);
/* Computing MIN */
		i__2 = k + mu;
		k2 = min(i__2, n);
		i__2 = k2;
		for (j = k1; j <= i__2; ++j) {
			if (j != k) {
				fjac[k + j * ldfjac] = -(one + two * x[j]);
			}
/* L470: */
		}
/* Computing 2nd power */
		d__1 = x[k];
		fjac[k + k * ldfjac] = two + fiftn * (d__1 * d__1);
/* L480: */
	}
 L490:
}

 /*
 *     subroutine errjac 
 *
 *     this subroutine is derived from vecjac which defines the 
 *     jacobian matrices of fourteen test functions. the problem 
 *     dimensions are as described in the prologue comments of vecfcn. 
 *     various errors are deliberately introduced to provide a test 
 *     for chkder. 
 *
 *     the subroutine statement is 
 *
 *       subroutine errjac(n,x,fjac,ldfjac,nprob) 
 *
 *     where 
 *
 *       n is a positive integer variable. 
 *
 *       x is an array of length n. 
 *
 *       fjac is an n by n array. on output fjac contains the 
 *         jacobian matrix, with various errors deliberately 
 *         introduced, of the nprob function evaluated at x. 
 *
 *       ldfjac is a positive integer variable not less than n 
 *         which specifies the leading dimension of the array fjac. 
 *
 *       nprob is a positive integer variable which defines the 
 *         number of the problem. nprob must not exceed 14. 
 *
 *     subprograms called 
 *
 *       fortran-supplied ... datan,dcos,dexp,dmin1,dsin,dsqrt, 
 *                            max0,min0 
 *
 *     argonne national laboratory. minpack project. march 1980. 
 *     burton s. garbow, kenneth e. hillstrom, jorge j. more 
 */
void errjac(int n, double *x, double *fjac, int ldfjac, int nprob)
{
    /* Initialized data */
    static double zero = 0.;
    static double one = 1.;
    static double two = 2.;
    static double three = 3.;
    static double four = 4.;
    static double five = 5.;
    static double six = 6.;
    static double eight = 8.;
    static double ten = 10.;
    static double fiftn = 15.;
    static double twenty = 20.;
    static double hundrd = 100.;
    static double c1 = 1e4;
    static double c3 = 200.;
    static double c4 = 20.2;
    static double c5 = 19.8;
    static double c6 = 180.;
    static double c9 = 29.;

    /* System generated locals */
    int fjac_offset, i__1, i__2, i__3, i__4;
    double d__1, d__2;

    /* Local variables */
    static double h__;
    static int i__, j, k, k1, k2, ml;
    static double ti, tj, tk;
    static int mu;
    static double tpi, sum, sum1, sum2, prod, temp, temp1, temp2, temp3, 
	    temp4;

    /* Parameter adjustments */
    --x;
    fjac_offset = 1 + ldfjac;
    fjac -= fjac_offset;

    /* Function Body */

/*     JACOBIAN ROUTINE SELECTOR. */

    switch (nprob) {
	case 1:  goto L10;
	case 2:  goto L20;
	case 3:  goto L50;
	case 4:  goto L60;
	case 5:  goto L90;
	case 6:  goto L100;
	case 7:  goto L200;
	case 8:  goto L230;
	case 9:  goto L290;
	case 10:  goto L320;
	case 11:  goto L350;
	case 12:  goto L380;
	case 13:  goto L420;
	case 14:  goto L450;
    }

/*     ROSENBROCK FUNCTION WITH SIGN REVERSAL AFFECTING ELEMENT (1,1). */

L10:
    fjac[ldfjac + 1] = one;
    fjac[(ldfjac << 1) + 1] = zero;
    fjac[ldfjac + 2] = -twenty * x[1];
    fjac[(ldfjac << 1) + 2] = ten;
    goto L490;

/*     POWELL SINGULAR FUNCTION WITH SIGN REVERSAL AFFECTING ELEMENT */
/*     (3,3). */

L20:
    for (k = 1; k <= 4; ++k) {
	for (j = 1; j <= 4; ++j) {
	    fjac[k + j * ldfjac] = zero;
/* L30: */
	}
/* L40: */
    }
    fjac[ldfjac + 1] = one;
    fjac[(ldfjac << 1) + 1] = ten;
    fjac[ldfjac * 3 + 2] = sqrt(five);
    fjac[(ldfjac << 2) + 2] = -fjac[ldfjac * 3 + 2];
    fjac[(ldfjac << 1) + 3] = two * (x[2] - two * x[3]);
    fjac[ldfjac * 3 + 3] = two * fjac[(ldfjac << 1) + 3];
    fjac[ldfjac + 4] = two * sqrt(ten) * (x[1] - x[4]);
    fjac[(ldfjac << 2) + 4] = -fjac[ldfjac + 4];
    goto L490;

/*     POWELL BADLY SCALED FUNCTION WITH THE SIGN OF THE JACOBIAN */
/*     REVERSED. */

L50:
    fjac[ldfjac + 1] = -c1 * x[2];
    fjac[(ldfjac << 1) + 1] = -c1 * x[1];
    fjac[ldfjac + 2] = exp(-x[1]);
    fjac[(ldfjac << 1) + 2] = exp(-x[2]);
    goto L490;

/*     WOOD FUNCTION WITHOUT ERROR. */

L60:
    for (k = 1; k <= 4; ++k) {
	for (j = 1; j <= 4; ++j) {
	    fjac[k + j * ldfjac] = zero;
/* L70: */
	}
/* L80: */
    }
/* Computing 2nd power */
    d__1 = x[1];
    temp1 = x[2] - three * (d__1 * d__1);
/* Computing 2nd power */
    d__1 = x[3];
    temp2 = x[4] - three * (d__1 * d__1);
    fjac[ldfjac + 1] = -c3 * temp1 + one;
    fjac[(ldfjac << 1) + 1] = -c3 * x[1];
    fjac[ldfjac + 2] = -two * c3 * x[1];
    fjac[(ldfjac << 1) + 2] = c3 + c4;
    fjac[(ldfjac << 2) + 2] = c5;
    fjac[ldfjac * 3 + 3] = -c6 * temp2 + one;
    fjac[(ldfjac << 2) + 3] = -c6 * x[3];
    fjac[(ldfjac << 1) + 4] = c5;
    fjac[ldfjac * 3 + 4] = -two * c6 * x[3];
    fjac[(ldfjac << 2) + 4] = c6 + c4;
    goto L490;

/*     HELICAL VALLEY FUNCTION WITH MULTIPLICATIVE ERROR AFFECTING */
/*     ELEMENTS (2,1) AND (2,2). */

L90:
    tpi = eight * atan(one);
/* Computing 2nd power */
    d__1 = x[1];
/* Computing 2nd power */
    d__2 = x[2];
    temp = d__1 * d__1 + d__2 * d__2;
    temp1 = tpi * temp;
    temp2 = sqrt(temp);
    fjac[ldfjac + 1] = hundrd * x[2] / temp1;
    fjac[(ldfjac << 1) + 1] = -hundrd * x[1] / temp1;
    fjac[ldfjac * 3 + 1] = ten;
    fjac[ldfjac + 2] = five * x[1] / temp2;
    fjac[(ldfjac << 1) + 2] = five * x[2] / temp2;
    fjac[ldfjac * 3 + 2] = zero;
    fjac[ldfjac + 3] = zero;
    fjac[(ldfjac << 1) + 3] = zero;
    fjac[ldfjac * 3 + 3] = one;
    goto L490;

/*     WATSON FUNCTION WITH SIGN REVERSALS AFFECTING THE COMPUTATION OF */
/*     TEMP1. */

L100:
    i__1 = n;
    for (k = 1; k <= i__1; ++k) {
	i__2 = n;
	for (j = k; j <= i__2; ++j) {
	    fjac[k + j * ldfjac] = zero;
/* L110: */
	}
/* L120: */
    }
    for (i__ = 1; i__ <= 29; ++i__) {
	ti = (double) i__ / c9;
	sum1 = zero;
	temp = one;
	i__1 = n;
	for (j = 2; j <= i__1; ++j) {
	    i__2 = j - 1;
	    sum1 += (double) i__2 * temp * x[j];
	    temp = ti * temp;
/* L130: */
	}
	sum2 = zero;
	temp = one;
	i__1 = n;
	for (j = 1; j <= i__1; ++j) {
	    sum2 += temp * x[j];
	    temp = ti * temp;
/* L140: */
	}
/* Computing 2nd power */
	d__1 = sum2;
	temp1 = two * (sum1 + d__1 * d__1 + one);
	temp2 = two * sum2;
/* Computing 2nd power */
	d__1 = ti;
	temp = d__1 * d__1;
	tk = one;
	i__1 = n;
	for (k = 1; k <= i__1; ++k) {
	    tj = tk;
	    i__2 = n;
	    for (j = k; j <= i__2; ++j) {
		i__3 = k - 1;
		i__4 = j - 1;
		fjac[k + j * ldfjac] += tj * (((double) i__3 / ti - 
			temp2) * ((double) i__4 / ti - temp2) - temp1);
		tj = ti * tj;
/* L150: */
	    }
	    tk = temp * tk;
/* L160: */
	}
/* L170: */
    }
/* Computing 2nd power */
    d__1 = x[1];
    fjac[ldfjac + 1] = fjac[ldfjac + 1] + six * (d__1 * d__1) - two * x[
	    2] + three;
    fjac[(ldfjac << 1) + 1] -= two * x[1];
    fjac[(ldfjac << 1) + 2] += one;
    i__1 = n;
    for (k = 1; k <= i__1; ++k) {
	i__2 = n;
	for (j = k; j <= i__2; ++j) {
	    fjac[j + k * ldfjac] = fjac[k + j * ldfjac];
/* L180: */
	}
/* L190: */
    }
    goto L490;

/*     CHEBYQUAD FUNCTION WITH JACOBIAN TWICE CORRECT SIZE. */

L200:
    tk = one / (double) (n);
    i__1 = n;
    for (j = 1; j <= i__1; ++j) {
	temp1 = one;
	temp2 = two * x[j] - one;
	temp = two * temp2;
	temp3 = zero;
	temp4 = two;
	i__2 = n;
	for (k = 1; k <= i__2; ++k) {
	    fjac[k + j * ldfjac] = two * tk * temp4;
	    ti = four * temp2 + temp * temp4 - temp3;
	    temp3 = temp4;
	    temp4 = ti;
	    ti = temp * temp2 - temp1;
	    temp1 = temp2;
	    temp2 = ti;
/* L210: */
	}
/* L220: */
    }
    goto L490;

/*     BROWN ALMOST-LINEAR FUNCTION WITHOUT ERROR. */

L230:
    prod = one;
    i__1 = n;
    for (j = 1; j <= i__1; ++j) {
	prod = x[j] * prod;
	i__2 = n;
	for (k = 1; k <= i__2; ++k) {
	    fjac[k + j * ldfjac] = one;
/* L240: */
	}
	fjac[j + j * ldfjac] = two;
/* L250: */
    }
    i__1 = n;
    for (j = 1; j <= i__1; ++j) {
	temp = x[j];
	if (temp != zero) {
	    goto L270;
	}
	temp = one;
	prod = one;
	i__2 = n;
	for (k = 1; k <= i__2; ++k) {
	    if (k != j) {
		prod = x[k] * prod;
	    }
/* L260: */
	}
L270:
	fjac[n + j * ldfjac] = prod / temp;
/* L280: */
    }
    goto L490;

/*     DISCRETE BOUNDARY VALUE FUNCTION WITH MULTIPLICATIVE ERROR */
/*     AFFECTING THE JACOBIAN DIAGONAL. */

L290:
    i__1 = n + 1;
    h__ = one / (double) i__1;
    i__1 = n;
    for (k = 1; k <= i__1; ++k) {
/* Computing 2nd power */
	d__1 = x[k] + (double) k * h__ + one;
	temp = three * (d__1 * d__1);
	i__2 = n;
	for (j = 1; j <= i__2; ++j) {
	    fjac[k + j * ldfjac] = zero;
/* L300: */
	}
/* Computing 2nd power */
	d__1 = h__;
	fjac[k + k * ldfjac] = four + temp * (d__1 * d__1);
	if (k != 1) {
	    fjac[k + (k - 1) * ldfjac] = -one;
	}
	if (k != n) {
	    fjac[k + (k + 1) * ldfjac] = -one;
	}
/* L310: */
    }
    goto L490;

/*     DISCRETE INTEGRAL EQUATION FUNCTION WITH SIGN ERROR AFFECTING */
/*     THE JACOBIAN DIAGONAL. */

L320:
    i__1 = n + 1;
    h__ = one / (double) i__1;
    i__1 = n;
    for (k = 1; k <= i__1; ++k) {
	tk = (double) k * h__;
	i__2 = n;
	for (j = 1; j <= i__2; ++j) {
	    tj = (double) j * h__;
/* Computing 2nd power */
	    d__1 = x[j] + tj + one;
	    temp = three * (d__1 * d__1);
/* Computing MIN */
	    d__1 = tj * (one - tk), d__2 = tk * (one - tj);
	    fjac[k + j * ldfjac] = h__ * min(d__1,d__2) * temp / two;
/* L330: */
	}
	fjac[k + k * ldfjac] -= one;
/* L340: */
    }
    goto L490;

/*     TRIGONOMETRIC FUNCTION WITH SIGN ERRORS AFFECTING THE */
/*     OFFDIAGONAL ELEMENTS OF THE JACOBIAN. */

L350:
    i__1 = n;
    for (j = 1; j <= i__1; ++j) {
	temp = sin(x[j]);
	i__2 = n;
	for (k = 1; k <= i__2; ++k) {
	    fjac[k + j * ldfjac] = -temp;
/* L360: */
	}
	i__2 = j + 1;
	fjac[j + j * ldfjac] = (double) i__2 * temp - cos(x[j]);
/* L370: */
    }
    goto L490;

/*     VARIABLY DIMENSIONED FUNCTION WITH OPERATION ERROR AFFECTING */
/*     THE UPPER TRIANGULAR ELEMENTS OF THE JACOBIAN. */

L380:
    sum = zero;
    i__1 = n;
    for (j = 1; j <= i__1; ++j) {
	sum += (double) j * (x[j] - one);
/* L390: */
    }
/* Computing 2nd power */
    d__1 = sum;
    temp = one + six * (d__1 * d__1);
    i__1 = n;
    for (k = 1; k <= i__1; ++k) {
	i__2 = n;
	for (j = k; j <= i__2; ++j) {
	    i__3 = k * j;
	    fjac[k + j * ldfjac] = (double) i__3 / temp;
	    fjac[j + k * ldfjac] = fjac[k + j * ldfjac];
/* L400: */
	}
	fjac[k + k * ldfjac] += one;
/* L410: */
    }
    goto L490;

/*     BROYDEN TRIDIAGONAL FUNCTION WITHOUT ERROR. */

L420:
    i__1 = n;
    for (k = 1; k <= i__1; ++k) {
	i__2 = n;
	for (j = 1; j <= i__2; ++j) {
	    fjac[k + j * ldfjac] = zero;
/* L430: */
	}
	fjac[k + k * ldfjac] = three - four * x[k];
	if (k != 1) {
	    fjac[k + (k - 1) * ldfjac] = -one;
	}
	if (k != n) {
	    fjac[k + (k + 1) * ldfjac] = -two;
	}
/* L440: */
    }
    goto L490;

/*     BROYDEN BANDED FUNCTION WITH SIGN ERROR AFFECTING THE JACOBIAN */
/*     DIAGONAL. */

L450:
    ml = 5;
    mu = 1;
    i__1 = n;
    for (k = 1; k <= i__1; ++k) {
	i__2 = n;
	for (j = 1; j <= i__2; ++j) {
	    fjac[k + j * ldfjac] = zero;
/* L460: */
	}
/* Computing MAX */
	i__2 = 1, i__3 = k - ml;
	k1 = max(i__2,i__3);
/* Computing MIN */
	i__2 = k + mu;
	k2 = min(i__2,n);
	i__2 = k2;
	for (j = k1; j <= i__2; ++j) {
	    if (j != k) {
		fjac[k + j * ldfjac] = -(one + two * x[j]);
	    }
/* L470: */
	}
/* Computing 2nd power */
	d__1 = x[k];
	fjac[k + k * ldfjac] = two - fiftn * (d__1 * d__1);
/* L480: */
    }
L490:
}