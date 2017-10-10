// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/globaldef.hpp>

#include <sgpp/optimization/fuzzy/FuzzyExtensionPrincipleViaOptimization.hpp>
#include <sgpp/optimization/fuzzy/InterpolatedFuzzyInterval.hpp>

#include <limits>
#include <vector>

namespace sgpp {
namespace optimization {

namespace {
class ScaledScalarFunction : public ScalarFunction {
 public:
  explicit ScaledScalarFunction(const ScalarFunction& fOrig) :
    ScaledScalarFunction(
        fOrig,
        base::DataVector(fOrig.getNumberOfParameters()),
        base::DataVector(fOrig.getNumberOfParameters()),
        1.0) {}

  ScaledScalarFunction(const ScalarFunction& fOrig,
                       const base::DataVector& lowerBounds,
                       const base::DataVector& upperBounds,
                       double valueFactor) :
    ScalarFunction(fOrig.getNumberOfParameters()),
    lowerBounds(lowerBounds),
    upperBounds(upperBounds),
    valueFactor(valueFactor),
    xScaled(base::DataVector(d)) {
    fOrig.clone(this->fOrig);
  }

  ~ScaledScalarFunction() override {}

  inline double eval(const base::DataVector& x) override {
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

  base::DataVector& getLowerBounds() { return lowerBounds; }
  base::DataVector& getUpperBounds() { return upperBounds; }
  void setValueFactor(double valueFactor) { this->valueFactor = valueFactor; }

 protected:
  std::unique_ptr<ScalarFunction> fOrig;
  base::DataVector lowerBounds;
  base::DataVector upperBounds;
  double valueFactor;
  base::DataVector xScaled;
};

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

  base::DataVector& getLowerBounds() { return lowerBounds; }
  base::DataVector& getUpperBounds() { return upperBounds; }
  void setValueFactor(double valueFactor) { this->valueFactor = valueFactor; }

 protected:
  std::unique_ptr<ScalarFunctionGradient> fGradientOrig;
  base::DataVector lowerBounds;
  base::DataVector upperBounds;
  double valueFactor;
  base::DataVector xScaled;
};

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
            valueFactor * (upperBounds[t] - lowerBounds[t]) *
            valueFactor * (upperBounds[t2] - lowerBounds[t2]);
      }
    }

    return y;
  }

  void clone(std::unique_ptr<ScalarFunctionHessian>& clone) const override {
    clone = std::unique_ptr<ScalarFunctionHessian>(
        new ScaledScalarFunctionHessian(*fHessianOrig, lowerBounds, upperBounds, valueFactor));
  }

  base::DataVector& getLowerBounds() { return lowerBounds; }
  base::DataVector& getUpperBounds() { return upperBounds; }
  void setValueFactor(double valueFactor) { this->valueFactor = valueFactor; }

 protected:
  std::unique_ptr<ScalarFunctionHessian> fHessianOrig;
  base::DataVector lowerBounds;
  base::DataVector upperBounds;
  double valueFactor;
  base::DataVector xScaled;
};
}  // namespace

FuzzyExtensionPrincipleViaOptimization::FuzzyExtensionPrincipleViaOptimization(
    const ScalarFunction& f,
    size_t numberOfAlphaSegments) :
        FuzzyExtensionPrinciple(f),
        defaultOptimizer(optimizer::MultiStart(f)),
        m(numberOfAlphaSegments) {
  defaultOptimizer.clone(optimizer);
}

FuzzyExtensionPrincipleViaOptimization::FuzzyExtensionPrincipleViaOptimization(
    const optimizer::UnconstrainedOptimizer& optimizer,
    size_t numberOfAlphaSegments) :
      FuzzyExtensionPrinciple(optimizer.getObjectiveFunction()),
      defaultOptimizer(optimizer::MultiStart(optimizer.getObjectiveFunction())),
      m(numberOfAlphaSegments) {
  optimizer.clone(this->optimizer);

  if (optimizer.getObjectiveGradient() != nullptr) {
    optimizer.getObjectiveGradient()->clone(fGradient);
  }

  if (optimizer.getObjectiveHessian() != nullptr) {
    optimizer.getObjectiveHessian()->clone(fHessian);
  }
}

FuzzyExtensionPrincipleViaOptimization::FuzzyExtensionPrincipleViaOptimization(
    const FuzzyExtensionPrincipleViaOptimization& other) :
    FuzzyExtensionPrinciple(*other.f),
    defaultOptimizer(optimizer::MultiStart(*other.f)),
    m(other.m) {
  other.optimizer->clone(optimizer);

  if (other.fGradient.get() != nullptr) {
    other.fGradient->clone(fGradient);
  }

  if (fHessian.get() != nullptr) {
    other.fHessian->clone(fHessian);
  }
}

FuzzyExtensionPrincipleViaOptimization::~FuzzyExtensionPrincipleViaOptimization() {}

FuzzyInterval* FuzzyExtensionPrincipleViaOptimization::apply(
    const std::vector<const FuzzyInterval*>& x) const {
  const size_t d = f->getNumberOfParameters();
  ScaledScalarFunction fScaled(*f);
  std::unique_ptr<ScaledScalarFunctionGradient> fGradientScaled;
  std::unique_ptr<ScaledScalarFunctionHessian> fHessianScaled;

  // create scaled gradient function if gradient is given
  if (fGradient.get() != nullptr) {
    fGradientScaled.reset(new ScaledScalarFunctionGradient(*fGradient));
  }

  // create scaled Hessian function if Hessian is given
  if (fHessian.get() != nullptr) {
    fHessianScaled.reset(new ScaledScalarFunctionHessian(*fHessian));
  }

  // write restricted optimization domain directly to member of fScaled
  base::DataVector& lowerBounds(fScaled.getLowerBounds());
  base::DataVector& upperBounds(fScaled.getUpperBounds());

  // result data
  base::DataVector xData(2 * m + 2);
  base::DataVector alphaData(2 * m + 2);

  // save last optimal min/max value and corresponding argmin/argmax point
  // to make sure that the optimum for a smaller alpha (hence for a larger optimization domain)
  // is not worse that for a larger alpha
  double lastOptimalValueMin = std::numeric_limits<double>::infinity();
  double lastOptimalValueMax = -std::numeric_limits<double>::infinity();
  base::DataVector lastOptimalPointMin(d);
  base::DataVector lastOptimalPointMax(d);

  // iterate through alphas from 1 to 0
  for (size_t j = m + 1; j-- > 0;) {
    const double alpha = static_cast<double>(j) / static_cast<double>(m);

    // determine input parameter confidence interval,
    // directly changing the optimization domain in fScaled
    for (size_t t = 0; t < d; t++) {
      lowerBounds[t] = x[t]->evaluateConfidenceIntervalLowerBound(alpha);
      upperBounds[t] = x[t]->evaluateConfidenceIntervalUpperBound(alpha);
    }

    // set optimization domain in fGradientScaled if available
    if (fGradientScaled.get() != nullptr) {
      fGradientScaled->getLowerBounds() = lowerBounds;
      fGradientScaled->getUpperBounds() = upperBounds;
    }

    // set optimization domain in fHessianScaled if available
    if (fHessianScaled.get() != nullptr) {
      fHessianScaled->getLowerBounds() = lowerBounds;
      fHessianScaled->getUpperBounds() = upperBounds;
    }

    // compute minimum (lower bound of output confidence interval)
    {
      // value factor of 1 means minimization
      fScaled.setValueFactor(1.0);
      optimizer->setObjectiveFunction(fScaled);

      if (fGradientScaled.get() != nullptr) {
        fGradientScaled->setValueFactor(1.0);
        optimizer->setObjectiveGradient(fGradientScaled.get());
      }

      if (fHessianScaled.get() != nullptr) {
        fHessianScaled->setValueFactor(1.0);
        optimizer->setObjectiveHessian(fHessianScaled.get());
      }

      if (j < m) {
        // set starting point of optimizer to optimal point of the larger alpha
        // (optimizer is allowed to ignore the starting point though)
        optimizer->setStartingPoint(lastOptimalPointMin);
      }

      optimizer->optimize();

      // check if optimal value is indeed smaller than the last minimum
      if (optimizer->getOptimalValue() < lastOptimalValueMin) {
        xData[j] = optimizer->getOptimalValue();
        lastOptimalValueMin = optimizer->getOptimalValue();
        lastOptimalPointMin = optimizer->getOptimalPoint();
      } else {
        xData[j] = lastOptimalValueMin;
      }

      alphaData[j] = alpha;
    }

    // compute maximum (upper bound of output confidence interval)
    {
      // value factor of -1 means maximization
      fScaled.setValueFactor(-1.0);
      optimizer->setObjectiveFunction(fScaled);

      if (fGradientScaled.get() != nullptr) {
        fGradientScaled->setValueFactor(-1.0);
        optimizer->setObjectiveGradient(fGradientScaled.get());
      }

      if (fHessianScaled.get() != nullptr) {
        fHessianScaled->setValueFactor(-1.0);
        optimizer->setObjectiveHessian(fHessianScaled.get());
      }

      if (j < m) {
        // set starting point of optimizer to optimal point of the larger alpha
        // (optimizer is allowed to ignore the starting point though)
        optimizer->setStartingPoint(lastOptimalPointMax);
      }

      optimizer->optimize();

      // check if optimal value is indeed larger than the last maximum
      // (don't forget the minus because the value factor is not incorporated in getOptimalValue)
      if (-optimizer->getOptimalValue() > lastOptimalValueMax) {
        xData[2*m+1-j] = -optimizer->getOptimalValue();
        lastOptimalValueMax = -optimizer->getOptimalValue();
        lastOptimalPointMax = optimizer->getOptimalPoint();
      } else {
        xData[2*m+1-j] = lastOptimalValueMax;
      }

      alphaData[2*m+1-j] = alpha;
    }
  }

  // interpolate between alpha data points
  return new InterpolatedFuzzyInterval(xData, alphaData);
}

size_t FuzzyExtensionPrincipleViaOptimization::getNumberOfAlphaSegments() const {
  return m;
}

void FuzzyExtensionPrincipleViaOptimization::setNumberOfAlphaSegments(
    size_t numberOfAlphaSegments) {
  m = numberOfAlphaSegments;
}

}  // namespace optimization
}  // namespace sgpp
