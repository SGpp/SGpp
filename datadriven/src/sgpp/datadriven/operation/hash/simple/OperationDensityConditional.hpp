// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef OPERATIONDENSITYCONDITIONAL_HPP
#define OPERATIONDENSITYCONDITIONAL_HPP

#include <sgpp/base/grid/Grid.hpp>
#include <cstring>

#include <sgpp/globaldef.hpp>


namespace SGPP {
namespace datadriven {

/**
 * Conditionalize Probability Density Function
 */

class OperationDensityConditional {
 public:
  OperationDensityConditional() {}
  virtual ~OperationDensityConditional() {}

  /**
   * Conditional (Density) Functions
   *
   * @param alpha Coefficient vector for current grid
   * @param mg Referenz of grid pointer
   * @param malpha Coefficient vector for new grid (mg). Will be resized.
   * @param mdim Marginalize in dimension mdim
   * @param xbar Point at which to conditionalize
   */
  virtual void doConditional(base::DataVector& alpha, base::Grid*& mg,
                             base::DataVector& malpha, unsigned int mdim, float_t xbar) = 0;
};

}
}
#endif /* OPERATIONDENSITYCONDITIONAL_HPP */