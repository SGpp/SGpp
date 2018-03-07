// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef MORTONORDERKERNEL_HPP
#define MORTONORDERKERNEL_HPP

#include <stdint.h>

///@cond DOXY_IGNORE // NOLINT()
namespace sgpp {
namespace datadriven {
namespace OpMultiEvalCudaDetail {

/// Kernel wrapper for the Morton order generator
void zorder(double* pos, size_t* perm, size_t m, size_t DIM);

}  // namespace OpMultiEvalCudaDetail
}  // namespace datadriven
}  // namespace sgpp
///@endcond // NOLINT()

#endif  // MORTONORDERKERNEL_HPP
