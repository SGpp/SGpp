// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/datadriven/activeSubspaces/EigenFunctionalities.hpp>

#include <functional>
#include <iostream>
#include <map>
#include <tuple>
#include <unordered_map>
#include <utility>

namespace sgpp {
namespace datadriven {

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
  explicit MSplineBasis(sgpp::base::DataVector xi) : xi(xi) {}

  /**
   * Constructor
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
  double eval(size_t degree, size_t index, double x);

  double xpowplus(double x, size_t n);

  double w(size_t v, Eigen::VectorXd xi);

  /**
   * This is the original Schoenberg way to evaluate M-Splines https://doi.org/10.1007/BF02788653
   * It is numerically not as favourable as the iterative version and only here for testing.
   */
  double evalTruncated(double x, Eigen::VectorXd xi);

  void setXi(sgpp::base::DataVector xi) { this->xi = xi; }

 private:
  sgpp::base::DataVector xi;
  // tuple used as hash to store scalar products in innerProducts

  typedef std::tuple<size_t, size_t, double> mSplineHashType;

  // hash storage for scalar products. Holds all calculated scalar products s.t. they do not have to
  // calculated again if the same combination of indices, levels and dx is queried
  std::map<mSplineHashType, double> precalculatedValues;
  //  std::unordered_map<mSplineHashType, double, key_hash, key_equal> precalculatedValues;
};

}  // namespace datadriven
}  // namespace sgpp
