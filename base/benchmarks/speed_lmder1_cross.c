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
#include "eq.h"

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
	(void) m;
	if (*iflag == 1) {
		vecfcn(*n, x, fvec, nprob);
		nfev += 1;
	}
	if (*iflag == 2) {
		vecjac(*n, x, fjac, *ldfjac, nprob);
		njev += 1;
	}
}

void do_test(int nprob_, int n, double factor)
{
	double tol = sqrt(DBL_EPSILON);

	nprob = nprob_;
	nfev = 0;
	njev = 0;
	
	int lwa = 5*n + n;
	printf("# %s with n=%d\n", problem_name[nprob - 1], n);

	double x[n];		// solution
	int iwa[n];		    // integer work array for lmdif1

	double * wa = malloc(lwa * sizeof(*wa));	    // work array for lmdif1
	double * fvec = malloc(n * sizeof(*fvec));	    // residuals
	double * fjac = malloc(n * n * sizeof(*fjac));  // jacobian

	assert(wa != NULL);
	initpt(n, x, nprob, factor);	// set initial point

	int info = 0;
	double start = wtime();
	lmder1_(fcn, &n, &n, x, fvec, fjac, &n, &tol, &info, iwa, wa, &lwa);	// find solution
	double stop = wtime();

	vecfcn(n, x, fvec, nprob);	// evaluate residuals
	double fnorm2 = enorm_(&n, fvec);
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
	do_test(13, 3000, 1);
	do_test(11, 3200, 1);
	do_test(12, 1000, 1);
	printf("# total time: %.1fs\n", wtime() - start);
}