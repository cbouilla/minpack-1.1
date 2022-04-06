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

