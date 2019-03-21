/*
 * Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * scalapack.hpp
 *
 * Created on: Jan 15, 2019
 *     Author: Jan Schopohl
 */
#pragma once

#include <cstddef>

namespace sgpp {
namespace datadriven {

const int dlen_ = 9;
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

extern "C" {

/**
 *  ictxt   (global output) int
 *          ICTXT specifies the BLACS context handle identifying the
 *          created process grid.  The context itself is global.
 *
 *  nprow   (global input) int
 *          NPROW specifies the number of process rows in the grid
 *          to be created.
 *
 *  npcol   (global input) int
 *          NPCOL specifies the number of process columns in the grid
 *          to be created.
 */
void sl_init_(int &ictxt, int &nprow, int &npcol);

void pdgsev(int *n, int *nrhs, double *a, int *ia, int *ja, int *desca, int *ipiv, double *b,
            int *ib, int *jb, int *descb, int *info);

void pdlaset_(char &uplo, int &m, int &n, double &alpha, double &beta, double *a, int &ia, int &ja,
              int *desca);

void pdlaprnt_(const int &m, const int &n, const double *a, const int &ia, const int &ja,
               const int *desca, const int &irprnt, const int &icprnt, const char *cmatnm,
               const int &nout, double *work);

void descinit_(int *desc, int m, int n, int mb, int nb, int irsrc, int icsrc, int ictxt, int lld,
               int info);

// ScaLAPACK tools

int numroc_(const int &n, const int &nb, const int &iproc, const int &isrcproc, const int &nprocs);

// redistribute matrix
void pdgemr2d_(const int &m, const int &n, const double *a, const int &ia, const int &ja,
               const int *desca, double *b, const int &ib, const int &jb, const int *descb,
               const int &ictxt);

void Cpdgemr2d(int m, int n, const double *a, int ia, int ja, const int *desca, double *b, int ib,
               int jb, const int *descb, int ictxt);

// Level 1 BLAS

// sub(y) := sub(y) + a*sub(x)
void pdaxpy_(const int &n, const double &a, const double *x, const int &ix, const int &jx,
             const int *descx, const int &incx, double *y, const int &iy, const int &jy,
             const int *descy, const int &incy);

// dot = sub(x)'*sub(y)
void pddot_(const int &n, double &dot, const double *x, const int &ix, const int &jx,
            const int *descx, const int &incx, const double *y, const int &iy, const int &jy,
            const int *descy, const int &incy);

// sub(x) = a*sub(x)
void pdscal_(const int &n, const double &a, double *x, const int &ix, const int &jx,
             const int *descx, const int &incx);

// Level 2 BLAS
// sub(y) := alpha*sub(A)*sub(x) + beta*sub(y)
void pdgemv_(const char *trans, const int &m, const int &n, const double &alpha, const double *a,
             const int &ia, const int &ja, const int *desca, const double *x, const int &ix,
             const int &jx, const int *descx, const int &incx, const double &beta, double *y,
             const int &iy, const int &jy, const int *descy, const int &incy);
// sub(A) := alpha*sub(x)*sub(y)' + sub(A)
void pdger_(const int &m, const int &n, const double &alpha, const double *x, const int &ix,
            const int &jx, const int *descx, const int &incx, const double &beta, const double *y,
            const int &iy, const int &jy, const int *descy, const int &incy, double *a,
            const int &ia, const int &ja, const int *desca);

// Level 3 BLAS
// sub(C) := alpha*op(sub(A))*op(sub(B)) + beta*sub(C)
void pdgemm_(const char *transa, const char *transb, const int &m, const int &n, const int &k,
             const double &alpha, const double *a, const int &ia, const int &ja, const int *desca,
             const double *b, const int &ib, const int &jb, const int *descb, const double &beta,
             double *c, const int &ic, const int &jc, const int *descc);

// sub(C):=beta*sub(C) + alpha*op(sub(A))
void pdgeadd_(const char *trans, const int &m, const int &n, const double &alpha, const double *a,
              const int &ia, const int &ja, const int *desca, const double &beta, double *c,
              const int &ic, const int &jc, const int *descc);

// sub(C):=beta*sub(C) + alpha*sub(A)'
void pdtran_(const int &m, const int &n, const double &alpha, const double *a, const int &ia,
             const int &ja, const int *desca, const double &beta, double *c, const int &ic,
             const int &jc, const int *descc);
}
}  // namespace datadriven
}  // namespace sgpp