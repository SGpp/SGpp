// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef OPERATIONDENSITYMARGTO1DLINEAR_HPP
#define OPERATIONDENSITYMARGTO1DLINEAR_HPP

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/datadriven/operation/hash/simple/OperationDensityMargTo1D.hpp>

#include <sgpp/globaldef.hpp>
#include <cstring>

namespace SGPP {
namespace datadriven {

/**
 * keep applying marginalize to function until it's reduced to only 1 dimension
 */

class OperationDensityMargTo1DLinear : public OperationDensityMargTo1D {
 public:
  explicit OperationDensityMargTo1DLinear(base::Grid* grid) : grid(grid) {}
  virtual ~OperationDensityMargTo1DLinear() {}

  /**
   * Keep applying marginalizes to (Density) Functions, until it's reduced to 1 dimension (dim_x)
   *
   * @param alpha Coefficient vector for current grid
   * @param grid_x output 1D grid pointer
   * @param alpha_x Coefficient vector for new grid (grid_x). Will be initialized.
   * @param dim_x Target dimension, all other dimensions will be marginalized
   */
  void margToDimX(base::DataVector* alpha, base::Grid*& grid_x, base::DataVector*& alpha_x,
                  size_t dim_x);

 protected:
  base::Grid* grid;
  void marg_next_dim(base::Grid* g_in, base::DataVector* a_in, base::Grid*& g_out,
                     base::DataVector*& a_out, size_t dims, size_t dim_x, size_t& count);
};
}  // namespace datadriven
}  // namespace SGPP
#endif /* OPERATIONDENSITYMARGTO1DLINEAR_HPP */
