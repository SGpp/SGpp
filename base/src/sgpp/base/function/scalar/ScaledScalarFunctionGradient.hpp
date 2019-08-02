// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/base/function/scalar/ScalarFunctionGradient.hpp>
#include <sgpp/globaldef.hpp>

#include <cstddef>
#include <memory>
#include <functional>

namespace sgpp {
namespace base {

class ScaledScalarFunctionGradient : public ScalarFunctionGradient {
 public:
  explicit ScaledScalarFunctionGradient(const ScalarFunctionGradient& fGradientOrig) :
    ScaledScalarFunctionGradient(
        fGradientOrig,
        DataVector(fGradientOrig.getNumberOfParameters()),
        DataVector(fGradientOrig.getNumberOfParameters()),
        1.0) {}

  ScaledScalarFunctionGradient(const ScalarFunctionGradient& fGradientOrig,
                               const DataVector& lowerBounds,
                               const DataVector& upperBounds,
                               double valueFactor) :
    ScalarFunctionGradient(fGradientOrig.getNumberOfParameters()),
    lowerBounds(lowerBounds),
    upperBounds(upperBounds),
    valueFactor(valueFactor),
    xScaled(DataVector(d)) {
    fGradientOrig.clone(this->fGradientOrig);
  }

  ~ScaledScalarFunctionGradient() override {}

  inline double eval(const DataVector& x, DataVector& gradient) override {
    // scale x from restricted domain
    for (size_t t = 0; t < d; t++) {
      xScaled[t] = lowerBounds[t] + x[t] * (upperBounds[t] - lowerBounds[t]);
    }

    // multiply with valueFactor
    const double y = valueFactor * fGradientOrig->eval(xScaled, gradient);

    for (size_t t = 0; t < d; t++) {
      // scale gradient
      gradient[t] *= valueFactor * (upperBounds[t] - lowerBounds[t]);
    }

    return y;
  }

  void clone(std::unique_ptr<ScalarFunctionGradient>& clone) const override {
    clone = std::unique_ptr<ScalarFunctionGradient>(
        new ScaledScalarFunctionGradient(*fGradientOrig, lowerBounds, upperBounds, valueFactor));
  }

  const DataVector& getLowerBounds() const { return lowerBounds; }
  void setLowerBounds(const DataVector& lowerBounds) { this->lowerBounds = lowerBounds; }

  const DataVector& getUpperBounds() const { return upperBounds; }
  void setUpperBounds(const DataVector& upperBounds) { this->upperBounds = upperBounds; }

  double getValueFactor(double valueFactor) const { return valueFactor; }
  void setValueFactor(double valueFactor) { this->valueFactor = valueFactor; }

 protected:
  std::unique_ptr<ScalarFunctionGradient> fGradientOrig;
  DataVector lowerBounds;
  DataVector upperBounds;
  double valueFactor;
  DataVector xScaled;
};
}  // namespace base
}  // namespace sgpp
