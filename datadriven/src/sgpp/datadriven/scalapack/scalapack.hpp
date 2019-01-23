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

void pdgemm(const char *transa, const char *transb, const int *m, const int *n, const int *k,
            const double *alpha, const double *a, const int *ia, const int *ja, const int *desca,
            const double *b, const int *ib, const int *jb, const int *descb, const double *beta,
            double *c, const int *ic, const int *jc, const int *descc);

void pdlaset_(char &uplo, int &m, int &n, double &alpha, double &beta, double *a, int &ia, int &ja,
              int *desca);

void pdlaprnt_(int &m, int &n, double *a, int &ia, int &ja, int *desca, int &irprnt, int &icprnt,
               const char *cmatnm, int &nout, double *work);

void descinit_(int *desc, int m, int n, int mb, int nb, int irsrc, int icsrc, int ictxt, int lld,
               int info);

int numroc_(const int &n, const int &nb, const int &iproc, const int &isrcproc, const int &nprocs);

void Cpdgemr2d(const int &m, const int &n, const double *a, const int &ia, const int &ja,
               const int *desca, double *b, const int &ib, const int &jb, const int *descb,
               const int &ictxt);
}
}  // namespace datadriven
}  // namespace sgpp