// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef OPERATIONFIRSTMOMENTBSPLINEBOUNDARY_HPP
#define OPERATIONFIRSTMOMENTBSPLINEBOUNDARY_HPP

#include <sgpp/base/operation/hash/OperationFirstMoment.hpp>
#include <sgpp/base/grid/Grid.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace base {

/**
 * FirstMomemnt of sparse grid function, linear grid without boundaries
 */
class OperationFirstMomentBsplineBoundary : public OperationFirstMoment {
 public:
  /**
   * Constructor of OperationFirstMomentBsplineBoundary
   *
   * @param grid Pointer to a sparse grid object
   */
  explicit OperationFirstMomentBsplineBoundary(Grid* grid) : grid(grid) {}

  ~OperationFirstMomentBsplineBoundary() override {}

  /**
   * Compute first moment of the function
   * @f[ \int_{\Omega} x\cdot f(x) dx. @f]
   *
   * @param alpha Coefficient vector for current grid
   * @param bounds describes the boundaries of the hypercube of the original function
   */
  double doQuadrature(const DataVector& alpha, DataMatrix* bounds = nullptr) override;

 protected:
  // Pointer to the grid object (Grid needed for getDegree() function)
  sgpp::base::Grid* grid;
};

}  // namespace base
}  // namespace sgpp

#endif /* OPERATIONFIRSTMOMENTBSPLINEBOUNDARY_HPP */
