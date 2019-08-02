// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/base/function/scalar/ScalarFunction.hpp>
#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace base {

class ScaledScalarFunction : public ScalarFunction {
 public:
  explicit ScaledScalarFunction(const ScalarFunction& fOrig) :
    ScaledScalarFunction(
        fOrig,
        DataVector(fOrig.getNumberOfParameters()),
        DataVector(fOrig.getNumberOfParameters()),
        1.0) {}

  ScaledScalarFunction(const ScalarFunction& fOrig,
                       const DataVector& lowerBounds,
                       const DataVector& upperBounds,
                       double valueFactor) :
    ScalarFunction(fOrig.getNumberOfParameters()),
    lowerBounds(lowerBounds),
    upperBounds(upperBounds),
    valueFactor(valueFactor),
    xScaled(DataVector(d)) {
    fOrig.clone(this->fOrig);
  }

  ~ScaledScalarFunction() override {}

  inline double eval(const DataVector& x) override {
    // scale x from restricted domain
    for (size_t t = 0; t < d; t++) {
      xScaled[t] = lowerBounds[t] + x[t] * (upperBounds[t] - lowerBounds[t]);
    }

    // multiply with valueFactor
    return valueFactor * fOrig->eval(xScaled);
  }

  void clone(std::unique_ptr<ScalarFunction>& clone) const override {
    clone = std::unique_ptr<ScalarFunction>(
        new ScaledScalarFunction(*fOrig, lowerBounds, upperBounds, valueFactor));
  }

  const DataVector& getLowerBounds() const { return lowerBounds; }
  void setLowerBounds(const DataVector& lowerBounds) { this->lowerBounds = lowerBounds; }

  const DataVector& getUpperBounds() const { return upperBounds; }
  void setUpperBounds(const DataVector& upperBounds) { this->upperBounds = upperBounds; }

  double getValueFactor(double valueFactor) const { return valueFactor; }
  void setValueFactor(double valueFactor) { this->valueFactor = valueFactor; }

 protected:
  std::unique_ptr<ScalarFunction> fOrig;
  DataVector lowerBounds;
  DataVector upperBounds;
  double valueFactor;
  DataVector xScaled;
};
}  // namespace base
}  // namespace sgpp
