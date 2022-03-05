#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "minpack.h"
#include "eq.h"

/*
 *     this program tests the ability of chkder to detect
 *     inconsistencies between functions and their first derivatives.
 *     fourteen test function vectors and jacobians are used. eleven of
 *     the tests are false(f), i.e. there are inconsistencies between
 *     the function vectors and the corresponding jacobians. three of
 *     the tests are true(t), i.e. there are no inconsistencies. the
 *     driver reads in data, calls chkder and prints out information
 *     required by and received from chkder.
 *
 *     argonne national laboratory. minpack project. march 1980.
 *     burton s. garbow, kenneth e. hillstrom, jorge j. more
 */


/* test cases */
struct test_case_chkder {
	int nprob;
	int n;
	int correct;
};

struct test_case_chkder tests_chkder[] = {
	{    1,    2, 	0},
	{    2,    4, 	0},
	{    3,    2, 	0},
	{    4,    4, 	1},
	{    5,    3, 	0},
	{    6,    9, 	0},
	{    7,    7, 	0},
	{    8,   10, 	1},
	{    9,   10, 	0},
	{   10,   10, 	0},
	{   11,   10, 	0},
	{   12,   10, 	0},
	{   13,   10, 	1},
	{   14,   10, 	0},
};


int main()
{
	printf("1..14\n"); // TAP protocol
	int ldfjac = 10;
	double cp = .123;
	double err[10];
	double fjac[100];
	double fvec1[10], fvec2[10];
	double x1[10], x2[10];
	 
	double errmax[14];
	double errmin[14];
                      
	for (int ic = 0; ic < 14; ic++) {
		int n = tests_chkder[ic].n;
		int nprob = tests_chkder[ic].nprob;

		initpt(n, x1, nprob, 1.);
		for (int i = 0; i < n; ++i) {
			x1[i] += cp;
			cp = -cp;
		}
		int mode = 1;
		chkder_(&n, &n, x1, fvec1, fjac, &ldfjac, x2, fvec2, &mode, err);
		mode = 2;
		vecfcn(n, x1, fvec1, nprob);
		errjac(n, x1, fjac, ldfjac, nprob);
		vecfcn(n, x2, fvec2, nprob);
		chkder_(&n, &n, x1, fvec1, fjac, &ldfjac, x2, fvec2, &mode, err);
		
		// compute minimum and maximum error
		errmin[ic] = err[0];
		errmax[ic] = err[0];
		for (int i = 0; i < n; ++i) {
			errmin[ic] = fmin(errmin[ic], err[i]);
			errmax[ic] = fmax(errmax[ic], err[i]);
		}

		int correct = tests_chkder[ic].correct;
		if (correct && errmin[ic] > 0.8)
			printf("ok %d - %s, n=%d, correct\n", ic+1, problem_name[nprob - 1], n);
		else if (!correct && errmin[ic] < 0.1)
			printf("ok %d - %s, n=%d, incorrect\n", ic+1, problem_name[nprob - 1], n);
		else
			printf("not ok %d - %s, n=%d, unexpected\n", ic+1, problem_name[nprob - 1], n);
	
		for (int i=0; i < n; i++)
			printf("\t# err[%2d] = %.20g\n", i, err[i]);
	}

	printf("\n");
	printf("################################\n");
	printf("# Summary of 14 tests of chkder\n");
	printf("# nprob   n    status     errmin         errmax\n");
	for (int i = 0; i < 14; ++i)
		printf("# %5d%5d%5d%15.3f%15.3f\n", tests_chkder[i].nprob, tests_chkder[i].n, 
			tests_chkder[i].correct, errmin[i], errmax[i]);
	return EXIT_SUCCESS;
}
