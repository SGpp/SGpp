// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <cstddef>

namespace sgpp {
namespace datadriven {

const int kDescriptor_ = 9;
const int dtype_ = 0;
const int ctxt_ = 1;
const int m_ = 2;
const int n_ = 3;
const int mb_ = 4;
const int nb_ = 5;
const int rsrc_ = 6;
const int csrc_ = 7;
const int lld_ = 8;

const char *const pblasNoTranspose = "N";
const char *const pblasTranspose = "T";
const char *const pblasConjugate = "C";

const char *const lowerTriangular = "L";
const char *const upperTriangular = "U";

#ifdef USE_SCALAPACK

extern "C" {

// ScaLAPACK tools

int numroc_(const size_t &n, const size_t &nb, const int &iproc, const int &isrcproc,
            const int &nprocs);

// redistribute matrix
void pdgemr2d_(const size_t &m, const size_t &n, const double *a, const size_t &ia,
               const size_t &ja, const int *desca, double *b, const size_t &ib, const size_t &jb,
               const int *descb, const int &ictxt);

// linear equations

// solve Ax=b using Cholesky decomposition
void pdpotrs_(const char *uplo, const size_t &n, const size_t &nrhs, const double *a,
              const size_t &ia, const size_t &ja, const int *desca, double *b, const size_t &ib,
              const size_t &jb, const int *descb, const int &info);

// Level 1 BLAS

// sub(y) := sub(y) + a*sub(x)
void pdaxpy_(const size_t &n, const double &a, const double *x, const size_t &ix, const size_t &jx,
             const int *descx, const size_t &incx, double *y, const size_t &iy, const size_t &jy,
             const int *descy, const size_t &incy);

// dot = sub(x)'*sub(y)
void pddot_(const size_t &n, double &dot, const double *x, const size_t &ix, const size_t &jx,
            const int *descx, const size_t &incx, const double *y, const size_t &iy,
            const size_t &jy, const int *descy, const size_t &incy);

// sub(x) = a*sub(x)
void pdscal_(const size_t &n, const double &a, double *x, const size_t &ix, const size_t &jx,
             const int *descx, const size_t &incx);

// Level 2 BLAS
// sub(y) := alpha*sub(A)*sub(x) + beta*sub(y)
void pdgemv_(const char *trans, const size_t &m, const size_t &n, const double &alpha,
             const double *a, const size_t &ia, const size_t &ja, const int *desca, const double *x,
             const size_t &ix, const size_t &jx, const int *descx, const size_t &incx,
             const double &beta, double *y, const size_t &iy, const size_t &jy, const int *descy,
             const size_t &incy);
// sub(A) := alpha*sub(x)*sub(y)' + sub(A)
void pdger_(const size_t &m, const size_t &n, const double &alpha, const double *x,
            const size_t &ix, const size_t &jx, const int *descx, const size_t &incx,
            const double &beta, const double *y, const size_t &iy, const size_t &jy,
            const int *descy, const size_t &incy, double *a, const size_t &ia, const size_t &ja,
            const int *desca);

// Level 3 BLAS
// sub(C) := alpha*op(sub(A))*op(sub(B)) + beta*sub(C)
void pdgemm_(const char *transa, const char *transb, const size_t &m, const size_t &n,
             const size_t &k, const double &alpha, const double *a, const size_t &ia,
             const size_t &ja, const int *desca, const double *b, const size_t &ib,
             const size_t &jb, const int *descb, const double &beta, double *c, const size_t &ic,
             const size_t &jc, const int *descc);

// sub(C):=beta*sub(C) + alpha*op(sub(A))
void pdgeadd_(const char *trans, const size_t &m, const size_t &n, const double &alpha,
              const double *a, const size_t &ia, const size_t &ja, const int *desca,
              const double &beta, double *c, const size_t &ic, const size_t &jc, const int *descc);

// sub(C):=beta*sub(C) + alpha*sub(A)'
void pdtran_(const size_t &m, const size_t &n, const double &alpha, const double *a,
             const size_t &ia, const size_t &ja, const int *desca, const double &beta, double *c,
             const size_t &ic, const size_t &jc, const int *descc);
}

#endif /* USE_SCALAPACK */

}  // namespace datadriven
}  // namespace sgpp
