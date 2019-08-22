// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef ELASTICNETFUNCTION_HPP
#define ELASTICNETFUNCTION_HPP

#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/solver/sle/fista/LassoFunction.hpp>
#include <sgpp/solver/sle/fista/RegularizationFunction.hpp>

#include <cmath>

namespace sgpp {
namespace solver {

/**
 * @brief The ElasticNetFunction class
 * @details Corresponds to the regularization functional
 * \f$ \left(1 - \gamma \right) \Vert \boldsymbol{\alpha} \Vert_2  + \gamma \Vert \boldsymbol{\alpha} \Vert _1 \f$.
 */
class ElasticNetFunction : public RegularizationFunction {
 public:
  // Accept lambda * [(1 - l1Ratio) ||x||_1 +  l1Ratio * |x|_2^2]
  // Internally we use ( lambda |x|_1 + gamma |x|_2^2)
    /**
   * @brief ElasticNetFunction
   * @param lambda controls the regularization strength.
   * @param l1Ratio (called \f$ \gamma \f$  above) controls the amount of \f$ l_1 \f$ regularization.
   * A value of one corresponds to the lasso, and a value of zero to the ridge regularization.
   */
  ElasticNetFunction(double lambda, double l1Ratio)
      : lambda(lambda * l1Ratio), gamma(lambda * (1 - l1Ratio)), lasso(this->lambda) {}

  ~ElasticNetFunction() override {}

  double eval(sgpp::base::DataVector weights) override {
    double l1 = 0.0;
    double l2 = 0.0;
#pragma omp parallel for reduction(+ : l1, l2)
    for (size_t i = 0; i < weights.getSize(); ++i) {
      l1 += std::abs(weights[i]);
      l2 += weights[i] * weights[i];
    }
    l1 *= lambda;
    l2 *= gamma;
    const double en = (l1 + l2);
    return en;
  }

  base::DataVector prox(const sgpp::base::DataVector& weights, double stepsize) override {
    base::DataVector proxVec = lasso.prox(weights, stepsize);
    const double multiplicator = 1 / (1 + 2 * stepsize * gamma);
    proxVec.mult(multiplicator);
    return proxVec;
  }

 private:
  const double lambda;
  const double gamma;
  LassoFunction lasso;
};

}  //  namespace solver
}  //  namespace sgpp

#endif  // ELASTICNETFUNCTION_HPP
