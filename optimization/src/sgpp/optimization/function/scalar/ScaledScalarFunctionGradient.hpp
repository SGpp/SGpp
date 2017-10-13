// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef SGPP_OPTIMIZATION_FUNCTION_SCALAR_SCALEDSCALARFUNCTIONGRADIENT_HPP
#define SGPP_OPTIMIZATION_FUNCTION_SCALAR_SCALEDSCALARFUNCTIONGRADIENT_HPP

#include <sgpp/globaldef.hpp>
#include <sgpp/optimization/function/scalar/ScalarFunctionGradient.hpp>

#include <cstddef>
#include <memory>
#include <functional>

namespace sgpp {
namespace optimization {

class ScaledScalarFunctionGradient : public ScalarFunctionGradient {
 public:
  explicit ScaledScalarFunctionGradient(const ScalarFunctionGradient& fGradientOrig) :
    ScaledScalarFunctionGradient(
        fGradientOrig,
        base::DataVector(fGradientOrig.getNumberOfParameters()),
        base::DataVector(fGradientOrig.getNumberOfParameters()),
        1.0) {}

  ScaledScalarFunctionGradient(const ScalarFunctionGradient& fGradientOrig,
                               const base::DataVector& lowerBounds,
                               const base::DataVector& upperBounds,
                               double valueFactor) :
    ScalarFunctionGradient(fGradientOrig.getNumberOfParameters()),
    lowerBounds(lowerBounds),
    upperBounds(upperBounds),
    valueFactor(valueFactor),
    xScaled(base::DataVector(d)) {
    fGradientOrig.clone(this->fGradientOrig);
  }

  ~ScaledScalarFunctionGradient() override {}

  inline double eval(const base::DataVector& x, base::DataVector& gradient) override {
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

  const base::DataVector& getLowerBounds() const { return lowerBounds; }
  void setLowerBounds(const base::DataVector& lowerBounds) { this->lowerBounds = lowerBounds; }

  const base::DataVector& getUpperBounds() const { return upperBounds; }
  void setUpperBounds(const base::DataVector& upperBounds) { this->upperBounds = upperBounds; }

  double getValueFactor(double valueFactor) const { return valueFactor; }
  void setValueFactor(double valueFactor) { this->valueFactor = valueFactor; }

 protected:
  std::unique_ptr<ScalarFunctionGradient> fGradientOrig;
  base::DataVector lowerBounds;
  base::DataVector upperBounds;
  double valueFactor;
  base::DataVector xScaled;
};
}  // namespace optimization
}  // namespace sgpp

#endif /* SGPP_OPTIMIZATION_FUNCTION_SCALAR_SCALEDSCALARFUNCTIONGRADIENT_HPP */
