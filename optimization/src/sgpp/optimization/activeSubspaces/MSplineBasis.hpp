// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/base/datatypes/DataVector.hpp>

#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Eigenvalues>
#include <iostream>

namespace sgpp {
namespace optimization {

/**
 * M-spline basis
 */
class MSplineBasis {
 public:
  /**
   * Constructor
   *
   * @param xi	knots defining the M-spline
   */
  MSplineBasis() {}

  /**
   * Destructor.
   */
  ~MSplineBasis() {}

  /*
   * This is the numerically more stable iterative way to evaluate M-Splines
   * https://projecteuclid.org/download/pdf_1/euclid.ss/1177012761
   *
   */
  double eval(size_t degree, size_t index, double x, sgpp::base::DataVector xi);

  double xpowplus(double x, size_t n);

  double w(size_t v, Eigen::VectorXd xi);

  /**
   * This is the original Schoenberg way to evaluate M-Splines https://doi.org/10.1007/BF02788653
   * It is numerically not as favourable as the iterative version and only here for testing.
   */
  double evalTruncated(double x, Eigen::VectorXd xi);

 protected:
};

}  // namespace optimization
}  // namespace sgpp
