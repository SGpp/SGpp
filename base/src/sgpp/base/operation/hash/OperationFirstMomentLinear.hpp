// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef OPERATIONFIRSTMOMENTLINEAR_HPP
#define OPERATIONFIRSTMOMENTLINEAR_HPP

#include <sgpp/base/operation/hash/OperationFirstMoment.hpp>
#include <sgpp/base/grid/Grid.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {
namespace base {

/**
 * FirstMomemnt of sparse grid function, linear grid without boundaries
 */
class OperationFirstMomentLinear : public OperationFirstMoment {
 public:
  /**
   * Constructor of OperationFirstMomentLinear
   *
   * @param storage Pointer to the grid's GridStorage object
   */
  explicit OperationFirstMomentLinear(GridStorage* storage) :
  storage(storage) {}

  ~OperationFirstMomentLinear() override {}

  /**
   * Compute first moment of the function
   * @f[ \int_{\Omega} x\cdot f(x) dx. @f]
   *
   * @param alpha Coefficient vector for current grid
   */
  float_t doQuadrature(const DataVector& alpha) override;

 protected:
  // Pointer to the grid's GridStorage object
  GridStorage* storage;
};

}  // namespace base
}  // namespace SGPP

#endif /* OPERATIONFIRSTMOMENTLINEAR_HPP */
