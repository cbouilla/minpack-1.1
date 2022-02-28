void initpt(int n, double * x, int nprob, double factor);
void vecfcn(int n, double * x, double * fvec, int nprob);
int vecjac(int n, double *x, double *fjac, int ldfjac, int nprob);

extern char * problem_name[];

struct test_case {
	int nprob;
	int n;
	double factor;
	int nfev;
	int info;
	double fnorm2;
};

extern struct test_case tests[];