// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef ZEROFUNCTION_HPP
#define ZEROFUNCTION_HPP

#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/solver/sle/fista/RegularizationFunction.hpp>

namespace sgpp {
namespace solver {
/**
 * @brief The ZeroFunction class
 * @details Corresponds to the trivial regularization function, which
 * is zero in all cases.
 */
class ZeroFunction : public RegularizationFunction {
 public:
  ~ZeroFunction() override {}
  double eval(base::DataVector weights) override { return 0.0; }

  base::DataVector prox(const base::DataVector& weights, double stepsize) override {
    return weights;
  }
};

}  //  namespace solver
}  //  namespace sgpp

#endif  // ZEROFUNCTION_HPP
