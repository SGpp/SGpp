// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef RIDGEFUNCTION_HPP
#define RIDGEFUNCTION_HPP

#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/solver/sle/fista/RegularizationFunction.hpp>

namespace sgpp {
namespace solver {
/**
 * @brief The RidgeFunction class
 * @details Corresponds to the regularization function
 * \f$ \Vert \boldsymbol{\alpha} \Vert_2^2\f$.
 */
class RidgeFunction : public RegularizationFunction {
 public:
    /**
   * @brief RidgeFunction
   * @param lambda controls the regularization strength.
   */
  explicit RidgeFunction(double lambda) : lambda(lambda) {}

  ~RidgeFunction() override {}

  // (lambda |x|_2 )
  double eval(sgpp::base::DataVector weights) override {
    return lambda * weights.dotProduct(weights);
  }

  base::DataVector prox(const sgpp::base::DataVector& weights, double stepsize) override {
    auto proxVec = base::DataVector(weights.getSize());
#pragma omp parallel for
    for (size_t i = 0; i < weights.getSize(); ++i) {
      double scale = 1 + 2 * stepsize * lambda;
      proxVec[i] = weights[i] / scale;
    }
    return proxVec;
  }

 private:
  const double lambda;
};

}  //  namespace solver
}  //  namespace sgpp

#endif  // RIDGEFUNCTION_HPP
