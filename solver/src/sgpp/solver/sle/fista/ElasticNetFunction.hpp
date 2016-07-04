// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef ELASTICNETFUNCTION_HPP
#define ELASTICNETFUNCTION_HPP

#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/solver/sle/fista/RegularizationFunction.hpp>
#include <sgpp/solver/sle/fista/LassoFunction.hpp>

#include <cmath>

namespace sgpp {
namespace solver {

class ElasticNetFunction : public RegularizationFunction {
 public:
  ElasticNetFunction(double lambda, double gamma) : lambda(lambda), gamma(gamma), lasso(lambda) {}

  // ( lambda |x|_1 + (gamma/2) |x|_2)
  double eval(sgpp::base::DataVector weights) const override {
    double l1 = 0.0;
    double l2 = 0.0;
#pragma omp parallel for reduction(+ : l1, l2)
    for (size_t i = 0; i < weights.getSize(); ++i) {
      l1 += std::abs(weights[i]);
      l2 += weights[i] * weights[i];
    }
    l1 *= lambda;
    l2 *= (gamma / 2.0);
    const double en = (l1 + l2);
    return en;
  }

  base::DataVector prox(const sgpp::base::DataVector& weights, double stepsize) const override {
    base::DataVector proxVec = lasso.prox(weights, stepsize);
    const double multiplicator = 1 / (1 + stepsize * gamma);
    proxVec.mult(multiplicator);
    return proxVec;
  }

 private:
  const double lambda;
  const double gamma;
  const LassoFunction lasso;
};

}  //  namespace solver
}  //  namepsace sgpp

#endif  // ELASTICNETFUNCTION_HPP
