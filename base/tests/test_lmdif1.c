/*
 * This program tests codes for the least-squares solution of
 * m nonlinear equations in n variables.
 *
 * Argonne National Laboratory. MINPACK project. march 1980.
 * Burton S. Garbow, Kenneth E. Hillstrom, Jorge j. Moré 
 */
#include <stdio.h>
#include <assert.h>
#include <float.h>
#include <math.h>
#include <stdbool.h>
#include <stdlib.h>

#include <minpack.h>
#include "ls.h"

int lmdif1_known_failures[] = { 26, 27, 38, 40, -1 };

/* global variables */
int nprob;
int nfev;
int njev;
int ma[60];
int na[60];
int nf[60];
int nj[60];
int np[60];
int nx[60];
double fnm[60];

/* This function is called by the solver and obeys the Fortran calling convention */
void fcn(int *m, int *n, double *x, double *fvec, int *iflag)
{
	ssqfcn(*m, *n, x, fvec, nprob);
	if (*iflag == 1)
		nfev += 1;
	if (*iflag == 2)
		njev += 1;
}

void do_test(int ic)
{
	double tol = sqrt(DBL_EPSILON);
	double ftol = 2e-5;

	// set global variables
	nprob = tests[ic].nprob;
	int n = tests[ic].n;
	int m = tests[ic].m;
	double factor = tests[ic].factor;

	double x[40];		// solution
	double fvec[100];	// residuals
	int iwa[40];		// integer work array for lmdif1
	double wa[2865];	// work array for lmdif1
	int lwa = 2865;

	initpt(n, x, nprob, factor);	// set initial point
	ssqfcn(m, n, x, fvec, nprob);	// evaluate residuals

	int info = 0;
	nfev = 0;
	njev = 0;
	lmdif1_(fcn, &m, &n, x, fvec, &tol, &info, iwa, wa, &lwa);	// find solution
	ssqfcn(m, n, x, fvec, nprob);	// evaluate residuals
	double fnorm2 = enorm_(&m, fvec);
	njev /= n;

	np[ic] = nprob;
	na[ic] = n;
	ma[ic] = m;
	nf[ic] = nfev;
	nj[ic] = njev;
	nx[ic] = info;
	fnm[ic] = fnorm2;

	// determine status and print it
	bool ok = (fnorm2 - tests[ic].fnorm2 < ftol);
	bool ok2 = (tests[ic].fnorm2_lastchance >= 0) && (fnorm2 - tests[ic].fnorm2_lastchance < ftol);

	if (ok || ok2) {
		printf("ok %d - %s (n = %d, m = %d, factor = %f)\n", ic + 1, problem_name[nprob - 1], n, m, factor);
	} else {
		bool TODO = false;
		for (int i = 0; lmdif1_known_failures[i] >= 0; i++)
			if (ic + 1 == lmdif1_known_failures[i])
				TODO = true;
		if (TODO)
			printf("not ok %d - %s (n = %d, m = %d, factor = %f) # TODO lmdif1 known failure\n", ic + 1, problem_name[nprob - 1], n, m, factor);
		else
			printf("not ok %d - %s (n = %d, m = %d, factor = %f)\n", ic + 1, problem_name[nprob - 1], n, m, factor);
	}
	commentator(ic, x, fvec, ftol, ftol, nfev, njev, info);
}

/* runs the usual collection of tests */
int main()
{
	/* TAP protocol */
	printf("1..53\n");

	for (int i = 0; i < 53; i++)
		do_test(i);

	printf("\n\n# Summary of 53 calls to lmdif1: \n\n");
	printf("#  test  nprob   n    m   nfev  njev  info  final l2 norm \n");
	for (int i = 0; i < 53; i++)
		printf("# %5d%5d%5d%5d%6d%6d%6d%16.7e\n", i + 1, np[i], na[i], ma[i], nf[i], nj[i], nx[i], fnm[i]);
	return EXIT_SUCCESS;
}
