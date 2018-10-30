// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/base/grid/generation/functors/SurplusRefinementFunctor.hpp>
#include <sgpp/base/grid/type/NakBsplineBoundaryGrid.hpp>
#include <sgpp/base/grid/type/NakBsplineGrid.hpp>
#include <sgpp/base/grid/type/NakBsplineModifiedGrid.hpp>
#include <sgpp/base/operation/hash/common/basis/NakBsplineBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/NakBsplineBoundaryBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/NakBsplineModifiedBasis.hpp>
#include <sgpp/optimization/activeSubspaces/ResponseSurface.hpp>
#include <sgpp/optimization/function/scalar/WrapperScalarFunction.hpp>
#include <sgpp/optimization/sle/solver/Auto.hpp>
#include <sgpp/optimization/sle/system/HierarchisationSLE.hpp>
#include <sgpp/optimization/tools/Printer.hpp>

#include <iostream>

namespace sgpp {
namespace optimization {

class SparseGridResponseSurfaceNakBspline : public ResponseSurface {
 public:
  SparseGridResponseSurfaceNakBspline(
      size_t dim, std::shared_ptr<sgpp::optimization::WrapperScalarFunction> objectiveFunc,
      sgpp::base::GridType gridType, size_t degree = 3)
      : ResponseSurface(dim), objectiveFunc(objectiveFunc), gridType(gridType), degree(degree) {
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
  void createSurplusAdaptiveResponseSurface(size_t maxNumGridPoints, size_t initialLevel);
  double eval(sgpp::base::DataVector v);
  double evalGradient(sgpp::base::DataVector v, sgpp::base::DataVector& gradient);

 private:
  std::shared_ptr<sgpp::optimization::WrapperScalarFunction> objectiveFunc;
  sgpp::base::GridType gridType;
  size_t degree;
  size_t numDim;
  std::unique_ptr<sgpp::base::Grid> grid;
  std::unique_ptr<sgpp::base::SBasis> basis;

  void refineSurplusAdaptive(size_t refinementsNum, sgpp::base::DataVector& alpha);
  sgpp::base::DataVector calculateInterpolationCoefficients();
};

}  // namespace optimization
}  // namespace sgpp
