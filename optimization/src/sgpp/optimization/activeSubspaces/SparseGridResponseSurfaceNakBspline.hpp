// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/base/grid/type/NakBsplineBoundaryGrid.hpp>
#include <sgpp/base/grid/type/NakBsplineGrid.hpp>
#include <sgpp/base/grid/type/NakBsplineModifiedGrid.hpp>
#include <sgpp/base/operation/hash/common/basis/NakBsplineBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/NakBsplineBoundaryBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/NakBsplineModifiedBasis.hpp>
#include <sgpp/optimization/activeSubspaces/ResponseSurface.hpp>
#include <sgpp/optimization/sle/solver/Auto.hpp>
#include <sgpp/optimization/sle/system/HierarchisationSLE.hpp>
#include <sgpp/optimization/tools/Printer.hpp>

#include <iostream>

namespace sgpp {
namespace optimization {

class SparseGridResponseSurfaceNakBspline : public ResponseSurface {
 public:
  SparseGridResponseSurfaceNakBspline(sgpp::optimization::WrapperScalarFunction objectiveFunc,
                                      sgpp::base::GridType gridType, size_t degree = 3)
      : objectiveFunc(objectiveFunc), gridType(gridType), degree(degree) {
    initialize();
  };

  /**
   * Destructor
   */
  virtual ~SparseGridResponseSurfaceNakBspline() {}

  /**
   * sets numDim, grid and basis according to objectiveFunction and gridType
   */
  void initialize();

  void createRegularResponseSurface(size_t level);
  double eval(sgpp::base::DataVector v);
  double evalGradient(sgpp::base::DataVector v, sgpp::base::DataVector& gradient);

 private:
  sgpp::optimization::WrapperScalarFunction objectiveFunc;
  sgpp::base::GridType gridType;
  size_t degree;
  size_t numDim;
  std::unique_ptr<sgpp::base::Grid> grid;
  std::unique_ptr<sgpp::base::SBasis> basis;
};

}  // namespace optimization
}  // namespace sgpp
