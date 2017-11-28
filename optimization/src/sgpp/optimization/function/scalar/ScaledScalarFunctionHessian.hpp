// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef SGPP_OPTIMIZATION_FUNCTION_SCALAR_SCALEDSCALARFUNCTIONHESSIAN_HPP
#define SGPP_OPTIMIZATION_FUNCTION_SCALAR_SCALEDSCALARFUNCTIONHESSIAN_HPP

#include <sgpp/globaldef.hpp>
#include <sgpp/optimization/function/scalar/ScalarFunctionHessian.hpp>

#include <cstddef>
#include <memory>
#include <functional>

namespace sgpp {
namespace optimization {

class ScaledScalarFunctionHessian : public ScalarFunctionHessian {
 public:
  explicit ScaledScalarFunctionHessian(const ScalarFunctionHessian& fHessianOrig) :
    ScaledScalarFunctionHessian(
        fHessianOrig,
        base::DataVector(fHessianOrig.getNumberOfParameters()),
        base::DataVector(fHessianOrig.getNumberOfParameters()),
        1.0) {}

  ScaledScalarFunctionHessian(const ScalarFunctionHessian& fHessianOrig,
                              const base::DataVector& lowerBounds,
                              const base::DataVector& upperBounds,
                              double valueFactor) :
    ScalarFunctionHessian(fHessianOrig.getNumberOfParameters()),
    lowerBounds(lowerBounds),
    upperBounds(upperBounds),
    valueFactor(valueFactor),
    xScaled(base::DataVector(d)) {
    fHessianOrig.clone(this->fHessianOrig);
  }

  ~ScaledScalarFunctionHessian() override {}

  inline double eval(const base::DataVector& x,
                     base::DataVector& gradient,
                     base::DataMatrix& hessian) override {
    // scale x from restricted domain
    for (size_t t = 0; t < d; t++) {
      xScaled[t] = lowerBounds[t] + x[t] * (upperBounds[t] - lowerBounds[t]);
    }

    // multiply with valueFactor
    const double y = valueFactor * fHessianOrig->eval(xScaled, gradient, hessian);

    for (size_t t = 0; t < d; t++) {
      // scale gradient
      gradient[t] *= valueFactor * (upperBounds[t] - lowerBounds[t]);

      for (size_t t2 = 0; t2 < d; t2++) {
        // scale Hessian
        hessian(t, t2) *=
            valueFactor * (upperBounds[t] - lowerBounds[t]) * (upperBounds[t2] - lowerBounds[t2]);
      }
    }

    return y;
  }

  void clone(std::unique_ptr<ScalarFunctionHessian>& clone) const override {
    clone = std::unique_ptr<ScalarFunctionHessian>(
        new ScaledScalarFunctionHessian(*fHessianOrig, lowerBounds, upperBounds, valueFactor));
  }

  const base::DataVector& getLowerBounds() const { return lowerBounds; }
  void setLowerBounds(const base::DataVector& lowerBounds) { this->lowerBounds = lowerBounds; }

  const base::DataVector& getUpperBounds() const { return upperBounds; }
  void setUpperBounds(const base::DataVector& upperBounds) { this->upperBounds = upperBounds; }

  double getValueFactor(double valueFactor) const { return valueFactor; }
  void setValueFactor(double valueFactor) { this->valueFactor = valueFactor; }

 protected:
  std::unique_ptr<ScalarFunctionHessian> fHessianOrig;
  base::DataVector lowerBounds;
  base::DataVector upperBounds;
  double valueFactor;
  base::DataVector xScaled;
};
}  // namespace optimization
}  // namespace sgpp

#endif /* SGPP_OPTIMIZATION_FUNCTION_SCALAR_SCALEDSCALARFUNCTIONHESSIAN_HPP */
