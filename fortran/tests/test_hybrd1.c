/*
 *     this program tests codes for the solution of n nonlinear
 *     equations in n variables.
 *
 *     Argonne National Laboratory. MINPACK project. march 1980.
 *     Burton S. Garbow, Kenneth E. Hillstrom, Jorge J. Moré
 */
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <stdbool.h>

#include "minpack.h"
#include "eq.h"

int hybrd1_known_failures[] = { 27, 28, 44, -1 };

/* global variables */
int nprob;
int nfev;
int njev;
int na[60];
int nf[60];
int np[60];
int nx[60];
double fnm[60];

/* This function is called by the solver and obeys the Fortran calling convention */
void fcn(const int *n, const double *x, double *fvec, int *iflag)
{
	(void)iflag;
	vecfcn(*n, x, fvec, nprob);
	nfev += 1;
}

void do_test(int ic)
{
	double tol = sqrt(MINPACK_EPSILON);
	double ftol = 1e-7;

	double x[50];		// solution
	double fvec[50];	// residuals
	double wa[4000];	// work array
	int lwa = 4000;

	// set problem instance
	nprob = tests[ic].nprob;
	int n = tests[ic].n;
	double factor = tests[ic].factor;

	initpt(n, x, nprob, factor);	// set initial point

	int info = 0;
	nfev = 0;
	hybrd1_(fcn, &n, x, fvec, &tol, &info, wa, &lwa);
	double fnorm2 = enorm_(&n, fvec);	// evaluate residuals

	np[ic] = nprob;
	na[ic] = n;
	nf[ic] = nfev;
	nx[ic] = info;
	fnm[ic] = fnorm2;

	// determine status and print it
	if (fnorm2 < ftol) {
		printf("ok %d - %s (n=%d, factor=%.0f)\n", ic + 1, problem_name[nprob - 1], n, factor);
	} else {
		bool TODO = false;
		for (int i = 0; hybrd1_known_failures[i] >= 0; i++)
			if (ic + 1 == hybrd1_known_failures[i])
				TODO = true;
		if (TODO)
			printf("not ok %d - %s (n=%d, factor=%.0f) # TODO hybrd1 known failure\n", ic + 1, problem_name[nprob - 1], n, factor);
		else
			printf("not ok %d - %s (n=%d, factor=%.0f)\n", ic + 1, problem_name[nprob - 1], n, factor);
	}
	printf("\t# norm of residual %15.7e\n", fnorm2);
}

/* runs the usual collection of tests */
int main(void)
{
	/* TAP protocol */
	printf("1..55\n");
	for (int ic = 0; ic < 55; ic++)
		do_test(ic);

	printf("\n\n# Summary of 55 calls to hybrd1: \n\n");
	printf("#  test  nprob   n    nfev   info  final l2 norm \n");
	for (int i = 0; i < 55; i++)
		printf("# %5d%5d%5d%5d%5d%16.7e\n", i + 1, np[i], na[i], nf[i], nx[i], fnm[i]);
	return EXIT_SUCCESS;
}
