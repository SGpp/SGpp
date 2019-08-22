// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef LASSOFUNCTION_HPP
#define LASSOFUNCTION_HPP

#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/solver/sle/fista/RegularizationFunction.hpp>

#include <algorithm>

namespace sgpp {
namespace solver {

/**
 * @brief The LassoFunction class.
 * @details This class models the standard lasso regularization penalty
 * \f$ \Vert \boldsymbol{\alpha} \Vert_1 \f$.
 */
class LassoFunction : public RegularizationFunction {
 public:
  /**
   * @brief LassoFunction
   * @param lambda controls the regularization strength.
   */
  explicit LassoFunction(double lambda) : lambda(lambda) {}

  ~LassoFunction() override {}

  double eval(base::DataVector weights) override {
    weights.abs();
    return lambda * weights.sum();
  }

  base::DataVector prox(const base::DataVector& weights, double stepsize) override {
    stepsize = lambda * stepsize;
    auto proxVec = base::DataVector(weights.getSize());
#pragma omp parallel for
    for (size_t i = 0; i < weights.getSize(); ++i) {
      const double left = std::max(weights[i] - stepsize, 0.0);
      const double right = std::max(-weights[i] - stepsize, 0.0);
      proxVec[i] = left - right;
    }
    return proxVec;
  }

 private:
  const double lambda;
};

}  //  namespace solver
}  //  namespace sgpp

#endif  // LASSOFUNCTION_HPP
