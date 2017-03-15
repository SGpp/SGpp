// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef MORTONORDERKERNEL_HPP
#define MORTONORDERKERNEL_HPP

#include <stdint.h>

/// Kernel wrapper for the Morton order generator
void zorder(double* pos, size_t* perm, size_t m, size_t DIM);

#endif  // MORTONORDERKERNEL_HPP
