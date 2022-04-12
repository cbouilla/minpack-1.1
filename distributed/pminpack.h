#include <float.h>

/*
 * for lmdif1 and lmdif
 *         calculate the functions at x and return this vector in fvec.
 *         One jacobian computation requires exactly n function calls
 */
typedef void (*pminpack_func_mn)(void *arg, int m, int n, double const *x, double *fvec);

int plmdif1(pminpack_func_mn fcn, void *farg, int m, int n, double *x, double *fvec, double tol, int ictx);


int plmdif(pminpack_func_mn fcn, void *farg, int m, int n, double *x, double *fvec, 
        double ftol, double xtol, double gtol, int maxfev, int *nfev, int ictx);

/* This replaces dpmpar */
#define MINPACK_EPSILON DBL_EPSILON
#define MINPACK_DWARF   DBL_MIN
#define MINPACK_GIANT   DBL_MAX

double enorm(int n, const double * x);

void forward_difference_jacobian(pminpack_func_mn fcn, void *farg, int m, int n, double *x,
	    double const *fvec, double *fjac, int ldfjac, double *wa, int ictx);

double lmpar(int n, double *r, int ldr, int *ipvt, double *diag, 
	double *qtb, double delta, double *x, double *sdiag, double *wa1, double *wa2);

void qrsolv(int n, double *r, int ldr, int *ipvt, double *diag, double *qtb, double *x, double *sdiag, double *wa);

double wtime();
void human_format(char * target, long n);
int scalapack_numroc(int n, int nb, int rank, int srcrank, int nprocs);
int scalapack_descinit(int *desc, int m, int n, int mb, int nb, int irsrc, int icsrc, int ictx, int lld);
int scalapack_pdgeqpf(int M, int N, double * A, int IA, int JA, int * DESCA, int * IPIV, double * TAU, double * WORK, int LWORK);
int scalapack_pdormqr(char * SIDE, char * TRANS, int M, int N, int K, double * A, int IA, int JA, int * DESCA, double *TAU,
                      double * C, int IC, int JC, int * DESCC, double * WORK, int LWORK);

void extrablacs_rvec2dmat(double *v, int vlen, double *A, int jA, int *descA);
void extrablacs_dmat2rmat(int m, int n, double *A, int iA, int jA, int *descA, double *B, int ldB);
void extrablacs_idmat2rmat(int m, int n, int *A, int iA, int jA, int *descA, int *B, int ldB);

/* scalapack functions */

#define DTYPE_ 0
#define CTXT_  1 
#define M_     2
#define N_     3 
#define MB_    4 
#define NB_    5
#define RSRC_  6 
#define CSRC_  7 
#define LLD_   8

/* ScaLAPACK functions */
extern int numroc_(const int *n, const int *nb, const int *rank, const int *srcrank, const int *nprocs);
extern void descinit_(int * desc, const int * m, const int * n, const int * mb, const int * nb, 
                      const int * irsrc, const int * icsrc, const int * ictx, const int * lld, int * info);
extern void pdgeqpf_(int * M, int * N, double * A, int * IA, int * JA, int * DESCA, int * IPIV, double * TAU, 
                     double * WORK, int * LWORK, int *INFO);
extern void pdormqr_(const char * SIDE, const char * TRANS, const int * M, const int *N, const int *K, const double * A, 
                     int * IA, int * JA, int * DESCA, double * TAU, double * C, int * IC, int * JC, int * DESCC, 
                     double * WORK, int * LWORK, int * INFO);
extern void Cpdgemr2d(int m, int n, double *A, int IA, int JA, int *descA, 
                      double *B, int IB, int JB, int *descB, int gcontext);

/* BLACS functions */
extern void Cblacs_get(int icontxt, int what, int *val);
extern void Cblacs_pinfo(int * rank, int * nprocs);
extern void Cblacs_gridinfo(int ictx, int * nprow, int * npcol, int * myrow, int * mycol);
extern void Cblacs_gridinit(int * ictx, const char * order, int nprow, int npcol);
extern void Cblacs_gridexit(int ictx);
extern void Cdgebs2d(int ctx, const char *scope, const char *top, int m, int n, double *A, int ldA);
extern void Cdgebr2d(int ctx, const char *scope, const char *top, int m, int n, double *A, int ldA, int rsrc, int csrc);