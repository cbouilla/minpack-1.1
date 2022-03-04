README for MINPACK
==================

Minpack includes software for solving nonlinear equations and
nonlinear least squares problems.  Five algorithmic paths each include
a core subroutine and an easy-to-use driver.  The algorithms proceed
either from an analytic specification of the Jacobian matrix or
directly from the problem functions.  The paths include facilities for
systems of equations with a banded Jacobian matrix, for least squares
problems with a large amount of data, and for checking the consistency
of the Jacobian matrix with the functions.

Jorge Moré, Burt Garbow, and Ken Hillstrom at Argonne National Laboratory.


References
==========

* M. J. D. Powell, *A Hybrid Method for Nonlinear Equations*.
  Numerical Methods for Nonlinear Algebraic Equations, P. Rabinowitz, editor. Gordon and Breach, 1970.

* Jorge J. Moré, *The Levenberg-Marquardt Algorithm, Implementation and Theory*. 
  Numerical Analysis, G. A. Watson, editor. Lecture Notes in Mathematics 630, Springer-Verlag, 1977.
  https://doi.org/10.1007/BFb0067700

* J. J. Moré, B. S. Garbow, and K. E. Hillstrom, *User Guide for MINPACK-1*, 
  Argonne National Laboratory Report ANL-80-74, Argonne, Ill., 1980.

* J. J. Moré, D. C. Sorensen, K. E. Hillstrom, and B. S. Garbow, *The MINPACK Project*, 
  in Sources and Development of Mathematical Software, W. J. Cowell, ed., Prentice-Hall, pages 88-111, 1984.

* J.J. Moré, B.S. Garbow and K.E. Hillstrom, *Testing unconstrained optimization software*. 
  ACM Trans. Math. Soft. 7, 1 (March 1981), 17-41. https://doi.org/10.1145/355934.355943


MINPACK versions
================

* The original MINPACK-1 Fortran version from 1980. https://www.netlib.org/minpack/
  
  Extensive original documentation. Extensive battery of tests (described in
  Algorithm 566 from ACM TOMS). The original test code processed all test
  problems and printed a summary that had to be manually inspected.


* MINPACK-2. https://ftp.mcs.anl.gov/pub/MINPACK-2/

  A MINPACK-2 project is announced. Contains some code and some new test problems,
  but is not backward compatible with the original MINPACK-1.


* CMINPACK from Manolis Lourakis.  http://www.netlib.org/minpack/cminpack.tar
  
  In 2002, Lourakis (lourakis at ics forth gr) published a C version derived
  from the fortran code using `f2c` and limited manual editing. This does not
  offer much more than directly linking against the original Fortran library.


* C/C++ MINPACK from Frédéric Devernay. http://devernay.free.fr/hacks/cminpack/
  
  Starting from 2007 and until 2021(at least), Devernay published much
  improved C versions (readable C code, friendlier C interface, no GOTOs,
  some tests, man pages, debian package, some CUDA support, CMAKE support,
  some LAPACK calls). A notable downside is that this does not run the
  original full test suite.
	
* Modernized Fortran code. https://github.com/fortran-lang/minpack

  Several people (including John Burkardt) have produced modernized Fortran
  versions of the original code (Fortran 90, free form, etc.).  Notably,
  there is an ongoing effort (as of 2022) to maintain a modern Fortran
  version of the original code.


In this repository
==================

This repository contains several versions of MINPACK-1 converted to in C and
potentially altered by Charles Bouillaguet.  All errors are mostly hiw own!

All versions should be drop-in replacement for the original Fortran code.

Tests
-----

The original 32 test cases have been ported to C, and modified to detect
success or failure. This uses the "Test Anything Protocol" that originates
from perl (using the `prove` test harness which is part of perl).  The new
set of tests is thus suitable for Continuous Integration.

It must be noted that the original Fortran code fails some tests. These
are "known failures".

Some limited benchmarking code has been added.

Testing MINPACK is difficult. A list of known issues is discussed below.



`fortran` version
-----------------

This is the original fortran code, along with the C tests.


`base` version
--------------

The code in the `base/` folder has been obtained by converting the original
fortran code to C (using f2c), restructuring the code to avoid most GOTOs,
and some manual editing.

The `lmstr` (non-linear least squares with limited storage) function has not
been included.

The main difference with the Fortran code is that the `dpmpar` function has
been removed, and replaced by constants from the `float.h` standard header.

There are tiny differences between the constants hardcoded in `dpmpar` and the
actual values (2.22044604926e-16 vs 2.2204460492503130808e-16 for machine epsilon). 
*These changes alter the behavior of the code*. It makes two tests fail
(does not converge).

The test suite has been updated to consider these as "known failures".

Using a different compiler alters the behavior of the code (e.g. tests fail
with `gcc`, succeed with `clang`). Changing compiler options may make some
test fails (e.g. pass with `-O1` but fail with `-O2`; actually observed with
`gcc` on POWER8 processors).


`lapackified` version
---------------------

The code in the `lapackified` folder differs from the `base/` code as follows:

- It uses LAPACK to compute QR factorizations and related operations. This
  reduces the amount of code (completely removes the `qrfac` function) and
  makes it much faster on large instances.

- It uses some BLAS functions when applicable, notably to reduce code size.

- The `mldif` and `lmder` functions are very similar; they have been merged.

- The `hybrid` and `hybrj` functions are very similar; they have been merged.

- Some functions have been converted to zero-indexed arrays.


All-in-all, the `lapackified` code is about 25% smaller and much faster than
the `base` code.


However, the actual numerical results will depend on the actual BLAS
implementation.  Changing the BLAS library may make some tests fail.  The
reference blas yields different numerical values than OpenBLAS. OpenBLAS is
multi-threaded; changing the number of threads alters the numerical values.
Running in `valgrind` using OpenBLAS changes the numerical values (this does
not happen with the reference BLAS). 
