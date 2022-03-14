#include <float.h>

/* 
 * for hybrd1 and hybrd:
 *    calculate the functions at x and return this vector in fvec.
 *   return a negative value to terminate hybrd1/hybrd.
 */
typedef int (*cminpack_func_n)(void *arg, int n, const double * x, double *fvec, int iflag);

/*
 * for hybrj1 and hybrj
 *         if iflag = 1 calculate the functions at x and return this vector in fvec. do not alter fjac.
 *         if iflag = 2 calculate the jacobian  at x and return this matrix in fjac. do not alter fvec.
 * return a negative value to terminate hybrj1/hybrj
 */
typedef int (*cminpack_func_nj)(void *arg, int n, const double *x, double *fvec, double *fjac, int ldfjac, int iflag);

/*
 * for lmdif1 and lmdif
 *         calculate the functions at x and return this vector in fvec.
 *         if iflag = 1 the result is used to compute the residuals.
 *         if iflag = 2 the result is used to compute the Jacobian by finite differences.
 *         One jacobian computation requires exactly n function calls with iflag = 2.
 * return a negative value to terminate lmdif1/lmdif
 */
typedef int (*cminpack_func_mn)(void *arg, int m, int n, double const *x, double *fvec, int iflag);

/*
 * for lmder1 and lmder
 *         if iflag = 1 calculate the functions at x and return this vector in fvec. do not alter fjac.
 *         if iflag = 2 calculate the jacobian  at x and return this matrix in fjac. do not alter fvec.
 * return a negative value to terminate lmder1/lmder
 */
typedef int (*cminpack_func_mnj)(void *arg, int m, int n, const double *x, double *fvec, double *fjac, int ldfjac, int iflag);

int lmder1(cminpack_func_mnj fcn, void *farg, int m, int n, double *x, 
	double *fvec, double *fjac, int ldfjac, double tol, 
	int *ipvt, double *wa, int lwa);

int lmder(cminpack_func_mnj fcn, void *fargs, int m, int n, double *x, 
	double *fvec, double *fjac, int ldfjac, double ftol,
	double xtol, double gtol, int maxfev, 
    double *diag, int mode, double factor, int nprint, 
     int *nfev, int *njev, int *ipvt, double *qtf, 
	double *wa1, double *wa2, double *wa3, double *wa4);

int lmdif1(cminpack_func_mn fcn, void *farg, int m, int n, double *x, 
	double *fvec, double tol, int *iwa, double *wa, int lwa);

int lmdif(cminpack_func_mn fcn, void *farg, int m, int n, double * x,
	   double * fvec, double ftol, double xtol, double gtol, int maxfev, double epsfcn, 
       double * diag,  int mode, double factor, int nprint,   int * nfev, double * fjac, int ldfjac, 
	   int * ipvt, double * qtf, double * wa1, double * wa2, double * wa3, double * wa4);

int hybrj1(cminpack_func_nj fcn, void *farg, int n, double *x, double *fvec, double *fjac,
	     int ldfjac, double tol, double *wa, int lwa);

int hybrj(cminpack_func_nj fcn, void *farg, int n, double *x, double *fvec,
	   double *fjac, int ldfjac, double xtol, int maxfev,
	   double *diag, int mode, double factor, int nprint,
	   int *nfev, int *njev, double *r, int lr, double *qtf, 
	   double *wa1, double *wa2, double *wa3, double *wa4);

int hybrd1(cminpack_func_n fcn, void *farg, int n, double *x, double *fvec, 
    double tol, double *wa, int lwa);

int hybrd(cminpack_func_n fcn, void *farg,
		  int n, double *x, double *fvec, double xtol, int maxfev,
		  double epsfcn, double *diag, int mode,
		  double factor, int nprint, int *nfev,
		  double *fjac, int ldfjac, double *r, int lr, double *qtf,
		  double *wa1, double *wa2, double *wa3, double *wa4);


/* This replaces dpmpar */
#define MINPACK_EPSILON DBL_EPSILON
#define MINPACK_DWARF   DBL_MIN
#define MINPACK_GIANT   DBL_MAX

void chkder(int m, int n, double * x, double * fvec, double * fjac, 
	int ldfjac, double * xp, double * fvecp, int mode, double * err);

void dogleg(int n, double * r, int lr, const double *diag, const double *qtb, 
		double delta, double *x, double *wa1, double *wa2);

double enorm(int n, const double * x);

int fdjac1(cminpack_func_n fcn, void *farg, int n, double *x, double *fvec,
	    double *fjac, int ldfjac, int iflag, double epsfcn, double *wa1, double *wa2);

int fdjac2(cminpack_func_mn fcn, void *farg, int m, int n, double *x,
	    double const *fvec, double *fjac, int ldfjac, int iflag,
	    double epsfcn, double *wa);

int hybrbase(cminpack_func_n fcn_dif, cminpack_func_nj fcn_der, void *farg,
	      int n, double *x, double *fvec, double *fjac, int ldfjac, 
	      double xtol, int maxfev, double epsfcn, double *diag, int mode,
	      double factor, int nprint, int *nfev, int *njev,
	      double *r, int lr, double *qtf,
	      double *wa1, double *wa2, double *wa3, double *wa4);

int r1updt(int m, int n, double *s, int ls, const double *u, double *v, double *w);
void r1mpyq(int m, int n, double *a, int lda, double *v, double *w);

double lmpar(int n, double *r, int ldr, int *ipvt, double *diag, 
	double *qtb, double delta, double *x, double *sdiag, double *wa1, double *wa2);

void qrsolv(int n, double *r, int ldr, int *ipvt, double *diag, double *qtb, double *x, double *sdiag, double *wa);

int lmbase(cminpack_func_mnj fcn_der, cminpack_func_mn fcn_dif, void *farg, 
	int m, int n, double *x, double *fvec, double *fjac, int ldfjac, double ftol,
	    double xtol, double gtol, int maxfev, double epsfcn, double *diag, 
	    int mode, double factor, int nprint, int *nfev, int *njev, int *ipvt, 
	    double *qtf, double *wa1, double *wa2, double *wa3, double *wa4);

