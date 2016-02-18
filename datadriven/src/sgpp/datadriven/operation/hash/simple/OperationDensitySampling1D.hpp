// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef OPERATIONDENSITYSAMPLING1D_HPP
#define OPERATIONDENSITYSAMPLING1D_HPP

#include <sgpp/base/grid/Grid.hpp>

#include <sgpp/globaldef.hpp>
#include <cstring>

namespace SGPP {
namespace datadriven {

/**
*Sample 1D Probability Density Function
*/

class OperationDensitySampling1D {
 public:
  OperationDensitySampling1D() {}
  virtual ~OperationDensitySampling1D() {}

  /**
   * Sampling on 1D grid
   *
   * @param alpha Coefficient vector for current grid (1D grid)
   * @param num_samples # of samples to draw
   * @param samples Output DataVector
  * @param seedp seed
   */
  virtual void doSampling1D(base::DataVector* alpha, size_t num_samples, base::DataVector*& samples,
                            unsigned int* seedp) = 0;
};
}  // namespace datadriven
}  // namespace SGPP
#endif /* OPERATIONDENSITYSAMPLING1D_HPP */
