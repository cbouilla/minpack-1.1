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
#include <sys/time.h>

#include "minpack.h"
#include "ls.h"

double wtime()
{
	struct timeval ts;
	gettimeofday(&ts, NULL);
	return (double)ts.tv_sec + ts.tv_usec / 1e6;
}

int nfev;
int njev;
int nprob;

/* This function is called by the solver and obeys the Fortran calling convention */
void fcn(int *m, int *n, double *x, double *fvec, double *fjac, int *ldfjac, int *iflag)
{
	if (*iflag == 1) {
		ssqfcn(*m, *n, x, fvec, nprob);
		nfev += 1;
	}
	if (*iflag == 2) {
		ssqjac(*m, *n, x, fjac, *ldfjac, nprob);
		njev += 1;
	}
}

void do_test(int nprob_, int n, int m, double factor)
{
	double tol = sqrt(DBL_EPSILON);

	nprob = nprob_;
	nfev = 0;
	njev = 0;

	int lwa = 5*n + m;
	printf("# %s with n=%d, m=%d\n", problem_name[nprob - 1], n, m);

	double x[n];		// solution
	int iwa[n];		    // integer work array for lmdif1

	double * wa = malloc(lwa * sizeof(*wa));	    // work array for lmdif1
	double * fvec = malloc(m * sizeof(*fvec));	    // residuals
	double * fjac = malloc(m * n * sizeof(*fjac));  // jacobian

	initpt(n, x, nprob, factor);	// set initial point

	int info = 0;
	double start = wtime();
	lmder1_(fcn, &m, &n, x, fvec, fjac, &m, &tol, &info, iwa, wa, &lwa);	// find solution
	double stop = wtime();

	ssqfcn(m, n, x, fvec, nprob);	// evaluate residuals
	double fnorm2 = enorm_(&m, fvec);
	njev /= n;
	free(wa);
	free(fvec);
	free(fjac);

	printf("# LMDER1: %.1fs\n", stop - start);
	printf("# function evaluations: %d\n", nfev);
	printf("# jacobian evaluations: %d\n", nfev);
	printf("# Final norm of residual: %15.7e\n", fnorm2);
	printf("# info: %d\n", info);
	printf("\n");
}


int main()
{
	double start = wtime();
	do_test(1, 1000, 1000, 1);
	do_test(1, 1000, 10000, 1);
	do_test(1, 1000, 100000, 1);
	do_test(3, 1000, 1000, 1);
	do_test(3, 1000, 10000, 1);
	do_test(16, 1000, 1000, 1);
	do_test(16, 3000, 3000, 1);
	do_test(12, 3, 1000000, 1);
	printf("# total time: %.1fs\n", wtime() - start);
}