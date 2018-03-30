// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef OPERATIONDENSITYMARGINALIZE_HPP
#define OPERATIONDENSITYMARGINALIZE_HPP

#include <sgpp/base/grid/Grid.hpp>

#include <sgpp/globaldef.hpp>
#include <cstring>

namespace sgpp {
namespace datadriven {

/**
 * Marginalize Probability Density Function
 */

class OperationDensityMarginalize {
 public:
  explicit OperationDensityMarginalize(base::Grid* grid) : grid(grid) {}
  virtual ~OperationDensityMarginalize() {}

  /**
   * Marginalizes (Density) Functions
   *
   * @param alpha Coefficient vector for current grid
   * @param mg Referenz of grid pointer
   * @param malpha Coefficient vector for new grid (mg). Will be resized.
   * @param mdim Marginalize in dimension mdim
   */
  virtual void doMarginalize(base::DataVector& alpha, base::Grid*& mg, base::DataVector& malpha,
                             unsigned int mdim);

 protected:
  base::Grid* grid;
};
}  // namespace datadriven
}  // namespace sgpp
#endif /* OPERATIONDENSITYMARGINALIZE_HPP */
