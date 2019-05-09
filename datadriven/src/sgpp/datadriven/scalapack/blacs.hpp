/*
 * Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * blacs.hpp
 *
 * Created on: Jan 15, 2019
 *     Author: Jan Schopohl
 */
#pragma once

#ifdef USE_SCALAPACK

#include <cstddef>

namespace sgpp {
namespace datadriven {
extern "C" {

// support routines

void Cblacs_pinfo(int &mypnum, int &nprocs);

int Cblacs_pnum(int icontxt, int prow, int pcol);

void Cblacs_get(int icontxt, int what, int &val);

void Cblacs_gridinit(int &icontxt, const char *order, int nprow, int npcol);

void Cblacs_gridinfo(int icontxt, int &nprow, int &npcol, int &myprow, int &mypcol);

void Cblacs_gridexit(int icontxt);

void Cblacs_exit(int cont);

void Cblacs_barrier(int icontxt, const char *scope);

// broadcasts

// send
void Cdgebs2d(int icontxt, const char *scope, const char *top, size_t m, size_t n, const double *a,
              size_t lda);

// receive
void Cdgebr2d(int icontxt, const char *scope, const char *top, size_t m, size_t n, double *a,
              size_t lda, int rsrc, int csrc);

// p2p send/receive

void Cdgesd2d(int icontxt, size_t m, size_t n, const double *a, size_t lda, int rdest, int cdest);

void Cdgerv2d(int icontxt, size_t m, size_t n, double *a, size_t lda, int rsrc, int csrc);
}
}  // namespace datadriven
}  // namespace sgpp

#endif /* USE_SCALAPACK */
