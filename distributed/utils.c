#include <stdio.h>
#include <sys/time.h>

#include "pminpack.h"

/******************* utility functions ********************/

double wtime()
{
        struct timeval ts;
        gettimeofday(&ts, NULL);
        return (double) ts.tv_sec + ts.tv_usec / 1e6;
}

/* represent n in <= 6 char  */
void human_format(char * target, long n)
{
        if (n < 1000) {
                sprintf(target, "%ld", n);
                return;
        }
        if (n < 1000000) {
                sprintf(target, "%.1fK", n / 1e3);
                return;
        }
        if (n < 1000000000) {
                sprintf(target, "%.1fM", n / 1e6);
                return;
        }
        if (n < 1000000000000ll) {
                sprintf(target, "%.1fG", n / 1e9);
                return;
        }
        if (n < 1000000000000000ll) {
                sprintf(target, "%.1fT", n / 1e12);
                return;
        }
}


int scalapack_numroc(int n, int nb, int rank, int srcrank, int nprocs)
{
        return numroc_(&n, &nb, &rank, &srcrank, &nprocs);
}
