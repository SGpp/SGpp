// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef OPERATIONSECONDMOMENTLINEAR_HPP
#define OPERATIONSECONDMOMENTLINEAR_HPP

#include <sgpp/base/operation/hash/OperationSecondMoment.hpp>
#include <sgpp/base/grid/Grid.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {
namespace base {

/**
 * SecondMomemnt of sparse grid function, linear grid without boundaries
 */
class OperationSecondMomentLinear : public OperationSecondMoment {
 public:
  /**
   * Constructor of OperationSecondMomentLinear
   *
   * @param storage Pointer to the grid's GridStorage object
   */
  OperationSecondMomentLinear(GridStorage* storage) : storage(storage) {}

  virtual ~OperationSecondMomentLinear() override {}

  /**
   * Compute second moment of the function
   * @f[ \int_{\Omega} x^2\cdot f(x) dx. @f]
   *
   * @param alpha Coefficient vector for current grid
   */
  float_t doQuadrature(DataVector& alpha) override;

 protected:
  // Pointer to the grid's GridStorage object
  GridStorage* storage;
};

}  // namespace base
}  // namespace SGPP

#endif /* OPERATIONSECONDMOMENTLINEAR_HPP */
