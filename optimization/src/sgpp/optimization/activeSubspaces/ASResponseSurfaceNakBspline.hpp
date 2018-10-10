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
#include <sgpp/optimization/activeSubspaces/ASResponseSurface.hpp>

namespace sgpp {
namespace optimization {

class ASResponseSurfaceNakBspline : public ASResponseSurface {
 public:
  ASResponseSurfaceNakBspline(Eigen::MatrixXd W1, sgpp::base::GridType gridType, size_t degree = 3)
      : ASResponseSurface(W1), degree(degree), gridType(gridType){};

  void createRegularSurfaceFromDetectionPoints(sgpp::base::DataMatrix evaluationPoints,
                                               sgpp::base::DataVector functionValues, size_t level);
  double eval(sgpp::base::DataVector v);
  double evalGradient(sgpp::base::DataVector v, sgpp::base::DataVector& gradient);

 private:
  size_t degree = 3;
  sgpp::base::GridType gridType;
  std::unique_ptr<sgpp::base::Grid> grid;
};

}  // namespace optimization
}  // namespace sgpp
