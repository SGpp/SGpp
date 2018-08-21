// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef OPERATIONSECONDMOMENTMODBSPLINECLENSHAWCURTIS_HPP
#define OPERATIONSECONDMOMENTMODBSPLINECLENSHAWCURTIS_HPP

#include <sgpp/base/operation/hash/OperationSecondMoment.hpp>
#include <sgpp/base/grid/Grid.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace base {

/**
 * FirstMomemnt of sparse grid function, linear grid without boundaries
 */
class OperationSecondMomentModBsplineClenshawCurtis : public OperationSecondMoment {
 public:
  /**
   * Constructor of OperationSecondMomentModBsplineClenshawCurtis
   *
   * @param grid Pointer to a sparse grid object
   */
  explicit OperationSecondMomentModBsplineClenshawCurtis(Grid* grid)
      : grid(grid), clenshawCurtisTable(base::ClenshawCurtisTable::getInstance()) {}

  ~OperationSecondMomentModBsplineClenshawCurtis() override {}

  /**
   * Compute first moment of the function
   * @f[ \int_{\Omega} x\cdot f(x) dx. @f]
   *
   * @param alpha Coefficient vector for current grid
   * @param bounds describes the boundaries of the hypercube of the original function
   */
  double doQuadrature(DataVector& alpha, DataMatrix* bounds = nullptr) override;

 protected:
  // Pointer to the grid object (Grid needed for getDegree() function)
  sgpp::base::Grid* grid;
  sgpp::base::ClenshawCurtisTable& clenshawCurtisTable;
};

}  // namespace base
}  // namespace sgpp

#endif /* OPERATIONSECONDMOMENTMODBSPLINECLENSHAWCURTIS_HPP */
