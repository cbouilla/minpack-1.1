#include <assert.h>
// #include <stdio.h>

#include "pminpack.h"

/*
 * Given a replicated vector of size n (available on all processes), copy it
 * to A[:, jA], where A is a distributed matrix.  The column index jA is
 * zero-based.
 */
void extrablacs_rvec2dmat(double *v, int vlen, double *A, int jA, int *descA)
{
	int ctx = descA[CTXT_];
	int mb = descA[MB_];
	int nb = descA[NB_];
	int lld = descA[LLD_];
	// int n = descA[N_];
	int m = descA[M_];

	assert(vlen == m);
	
	int nprow, npcol, myrow, mycol;
	Cblacs_gridinfo(ctx, &nprow, &npcol, &myrow, &mycol);

	int target_pcol = (jA / nb) % npcol; // process column that stores A[:,j]
	if (mycol != target_pcol)
		return;

	int localj = nb * (jA / (nb * npcol)) + (jA % nb); 
	int locali = 0;
	for (int i = myrow * mb; i < m; i += nprow * mb) {
		int hi = mb;
		if (i + hi > m)
			hi = m - i;
		for (int ii = 0; ii < hi; ii++)
			A[locali + ii + localj * lld] = v[i + ii];
		locali += mb;
	}
}


/*
 * Copy a  distributed (sub-)matrix A[iA:iA+m, jA:jA+n] to a local array B on all processes
 */
void extrablacs_dmat2rmat(int m, int n, double *A, int iA, int jA, int *descA, 
						  double *B, int ldB)
{
	int ctxA = descA[CTXT_];
	int mA = descA[M_];
	int nA = descA[N_];
	int mb = descA[MB_];
	int nb = descA[NB_];
	
	assert(ldB >= m);
	assert(iA + m - 1 <= mA);
	assert(jA + n - 1 <= nA);

	/* obtain the system context associated with context ctxA */
	int sysctx;
	Cblacs_get(ctxA, 10, &sysctx);

	/* create a "grid" for the single (0, 0) process */
	int ctxB = sysctx;
	Cblacs_gridinit(&ctxB, "Row-Major", 1, 1);

	int nprowB, npcolB, myrowB, mycolB;
	Cblacs_gridinfo(ctxB, &nprowB, &npcolB, &myrowB, &mycolB);
	int root = (myrowB == 0 && mycolB == 0);

	/* descriptor for the local array on the (0, 0) process */
	int descB[9];
	if (root)
		scalapack_descinit(descB, m, n, mb, nb, 0, 0, ctxB, ldB);
	else
		descB[CTXT_] = -1;

	/* copy the matrix to the process on the top-left corner of the grid */
	Cpdgemr2d(m, n, A, iA, jA, descA, B, 1, 1, descB, ctxA);

	/* broadcast the data to all other processes */
	if (root)
		Cdgebs2d(ctxA, "All", " ", m, n, B, ldB);
	else
		Cdgebr2d(ctxA, "All", " ", m, n, B, ldB, 0, 0);

	/* release the 1x1 "grid" */
	if (root)
		Cblacs_gridexit(ctxB);
}

void extrablacs_idmat2rmat(int m, int n, int *A, int iA, int jA, int *descA, 
						  int *B, int ldB)
{
	int ctxA = descA[CTXT_];
	int mA = descA[M_];
	int nA = descA[N_];
	int mb = descA[MB_];
	int nb = descA[NB_];
	
	assert(ldB >= m);
	assert(iA + m - 1 <= mA);
	assert(jA + n - 1 <= nA);

	/* obtain the system context associated with context ctxA */
	int sysctx;
	Cblacs_get(ctxA, 10, &sysctx);

	/* create a "grid" for the single (0, 0) process */
	int ctxB = sysctx;
	Cblacs_gridinit(&ctxB, "Row-Major", 1, 1);

	int nprowB, npcolB, myrowB, mycolB;
	Cblacs_gridinfo(ctxB, &nprowB, &npcolB, &myrowB, &mycolB);
	int root = (myrowB == 0 && mycolB == 0);

	/* descriptor for the local array on the (0, 0) process */
	int descB[9];
	if (root)
		scalapack_descinit(descB, m, n, mb, nb, 0, 0, ctxB, ldB);
	else
		descB[CTXT_] = -1;

	/* copy the matrix to the process on the top-left corner of the grid */
	Cpigemr2d(m, n, A, iA, jA, descA, B, 1, 1, descB, ctxA);

	/* broadcast the data to all other processes */
	if (root)
		Cigebs2d(ctxA, "All", " ", m, n, B, ldB);
	else
		Cigebr2d(ctxA, "All", " ", m, n, B, ldB, 0, 0);

	/* release the 1x1 "grid" */
	if (root)
		Cblacs_gridexit(ctxB);
}
