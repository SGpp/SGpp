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

class LassoFunction : public RegularizationFunction {
 public:
  explicit LassoFunction(double lambda) : lambda(lambda) {}

  double eval(base::DataVector weights) const override {
    weights.abs();
    return lambda * weights.sum();
  }

  base::DataVector prox(const base::DataVector& weights, double stepsize) const override {
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
}  //  namepsace sgpp

#endif  // LASSOFUNCTION_HPP
