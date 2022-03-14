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
void fcn(const int *m, const int *n, const double *x, double *fvec, int *iflag)
{
	ssqfcn(*m, *n, x, fvec, nprob);
	if (*iflag == 1)
		nfev += 1;
	if (*iflag == 2)
		njev += 1;
}

void do_test(int nprob_, int n, int m, double factor)
{
	double tol = sqrt(MINPACK_EPSILON);

	nprob = nprob_;
	nfev = 0;
	njev = 0;
	
	int lwa = m*n + 5*n + m;
	printf("# %s with n=%d, m=%d\n", problem_name[nprob - 1], n, m);

	double x[n];		// solution
	int iwa[n];		    // integer work array for lmdif1

	double * wa = malloc(lwa * sizeof(*wa));	    // work array for lmdif1
	double * fvec = malloc(m * sizeof(*fvec));	    // residuals
	assert(wa != NULL);
	initpt(n, x, nprob, factor);	// set initial point

	int info = 0;
	double start = wtime();
	lmdif1_(fcn, &m, &n, x, fvec, &tol, &info, iwa, wa, &lwa);	// find solution
	double stop = wtime();

	ssqfcn(m, n, x, fvec, nprob);	// evaluate residuals
	double fnorm2 = enorm_(&m, fvec);
	njev /= n;
	free(wa);
	free(fvec);

	printf("# LMDIF1: %.1fs\n", stop - start);
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
	// do_test(3, 1000, 1000, 1);
	// do_test(3, 1000, 10000, 1);
	do_test(16, 1000, 1000, 1);
	do_test(16, 3000, 3000, 1);
	do_test(12, 3, 1000000, 1);
	printf("# total time: %.1fs\n", wtime() - start);
}