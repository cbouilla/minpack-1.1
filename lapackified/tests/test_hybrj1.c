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

#include <minpack.h>
#include "eq.h"

int hybrj1_known_failures[] = { 2, 3, 18, 27, 28, 44, -1 };

/* global variables */
int nprob;
int nfev;
int njev;
int na[60];
int nf[60];
int nj[60];
int np[60];
int nx[60];
double fnm[60];

/* This function is called by the solver and obeys the Fortran calling convention */
void fcn(int *n, double *x, double *fvec, double *fjac, int *ldfjac, int *iflag)
{
	if (*iflag == 1) {
		vecfcn(*n, x, fvec, nprob);
		nfev += 1;
	}
	if (*iflag == 2) {
		vecjac(*n, x, fjac, *ldfjac, nprob);
		njev += 1;
	}
}

void do_test(int ic)
{
	double tol = sqrt(DBL_EPSILON);
	double ftol = 1e-7;

	double x[40];		// solution
	double fvec[40];	// residuals
	double wa[1060];	// work array
	double fjac[1600];	// jacobian
	int lwa = 1060;

	// set problem instance
	nprob = tests[ic].nprob;
	int n = tests[ic].n;
	double factor = tests[ic].factor;

	initpt(n, x, nprob, factor);	// set initial point

	int info = 0;
	nfev = 0;
	hybrj1_(fcn, &n, x, fvec, fjac, &n, &tol, &info, wa, &lwa);
	double fnorm2 = enorm_(&n, fvec);	// evaluate residuals

	np[ic] = nprob;
	na[ic] = n;
	nf[ic] = nfev;
	nj[ic] = njev;
	nx[ic] = info;
	fnm[ic] = fnorm2;

	// determine status and print it
	if (fnorm2 < ftol) {
		printf("ok %d - %s (n=%d, factor=%.0f)\n", ic + 1, problem_name[nprob - 1], n, factor);
	} else {
		bool TODO = false;
		for (int i = 0; hybrj1_known_failures[i] >= 0; i++)
			if (ic + 1 == hybrj1_known_failures[i])
				TODO = true;
		if (TODO)
			printf("not ok %d - %s (n=%d, factor=%.0f) # TODO hybrj1 known failure\n", ic + 1, problem_name[nprob - 1], n, factor);
		else
			printf("not ok %d - %s (n=%d, factor=%.0f)\n", ic + 1, problem_name[nprob - 1], n, factor);
	}
	if (info != tests[ic].info)
		printf("\t# NOTICE: mismatching status. Got %d, expected %d\n", info, tests[ic].info);
	if (info != 1)
		printf("\t# WARNING: solver claimed failure (code %d)\n", info);

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
	printf("#  test  nprob   n    nfev   njev   info  final l2 norm \n");
	for (int i = 0; i < 55; i++)
		printf("# %5d%5d%5d%5d%5d%5d%16.7e\n", i + 1, np[i], na[i], nf[i], nj[i], nx[i], fnm[i]);
	return EXIT_SUCCESS;
}
