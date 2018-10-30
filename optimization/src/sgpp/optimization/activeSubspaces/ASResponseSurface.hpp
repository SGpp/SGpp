// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once
//#ifdef USE_EIGEN

#include <sgpp/optimization/activeSubspaces/ResponseSurface.hpp>

#include <iostream>

namespace sgpp {
namespace optimization {

class ASResponseSurface : public ResponseSurface {
 public:
  ASResponseSurface(size_t dim, Eigen::MatrixXd W1) : ResponseSurface(dim), W1(W1){};

  /**
   * Destructor
   */
  virtual ~ASResponseSurface() {}

  virtual void createRegularReducedSurfaceFromDetectionPoints(
      sgpp::base::DataMatrix evaluationPoints, sgpp::base::DataVector functionValues,
      size_t level) = 0;

 protected:
  Eigen::MatrixXd W1;
};

}  // namespace optimization
}  // namespace sgpp

//#endif /* USE_EIGEN */
