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

namespace sgpp {
namespace datadriven {
extern "C" {

void blacs_pinfo_(int &mypnum, int &nprocs);

void blacs_get_(int &icontxt, const int &what, int &val);

void blacs_gridinit_(int &icontxt, const char *order, const int &nprow, int &npcol);

void blacs_gridinfo_(int &icontxt, int &nprow, int &npcol, int &myprow, int &mypcol);

void blacs_gridexit_(const int &icontxt);

void blacs_exit_(const int &cont);

void blacs_barrier(int &icontxt, const char *scope);
}
}  // namespace datadriven
}  // namespace sgpp