// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/globaldef.hpp>

#include <sgpp/optimization/fuzzy/FuzzyExtensionPrinciple.hpp>
#include <sgpp/optimization/fuzzy/InterpolatedFuzzyInterval.hpp>

#include <vector>

namespace sgpp {
namespace optimization {

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
    for (size_t t = 0; t < d; t++) {
      xScaled[t] = lowerBounds[t] + x[t] * (upperBounds[t] - lowerBounds[t]);
    }

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
    for (size_t t = 0; t < d; t++) {
      xScaled[t] = lowerBounds[t] + x[t] * (upperBounds[t] - lowerBounds[t]);
    }

    const double y = valueFactor * fGradientOrig->eval(xScaled, gradient);

    for (size_t t = 0; t < d; t++) {
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
    for (size_t t = 0; t < d; t++) {
      xScaled[t] = lowerBounds[t] + x[t] * (upperBounds[t] - lowerBounds[t]);
    }

    const double y = valueFactor * fHessianOrig->eval(xScaled, gradient, hessian);

    for (size_t t = 0; t < d; t++) {
      gradient[t] *= valueFactor * (upperBounds[t] - lowerBounds[t]);

      for (size_t t2 = 0; t2 < d; t2++) {
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

FuzzyExtensionPrinciple::FuzzyExtensionPrinciple(
    const ScalarFunction& f,
    size_t numberOfAlphaSegments) :
        defaultOptimizer(optimizer::MultiStart(f)),
        m(numberOfAlphaSegments) {
  defaultOptimizer.clone(optimizer);
  f.clone(this->f);
}

FuzzyExtensionPrinciple::FuzzyExtensionPrinciple(
    const optimizer::UnconstrainedOptimizer& optimizer,
    size_t numberOfAlphaSegments) :
      defaultOptimizer(optimizer::MultiStart(optimizer.getObjectiveFunction())),
      m(numberOfAlphaSegments) {
  optimizer.clone(this->optimizer);
  optimizer.getObjectiveFunction().clone(this->f);

  if (optimizer.getObjectiveGradient() != nullptr) {
    optimizer.getObjectiveGradient()->clone(fGradient);
  }

  if (optimizer.getObjectiveHessian() != nullptr) {
    optimizer.getObjectiveHessian()->clone(fHessian);
  }
}

FuzzyExtensionPrinciple::FuzzyExtensionPrinciple(const FuzzyExtensionPrinciple& other) :
    defaultOptimizer(optimizer::MultiStart(*other.f)),
    m(other.m) {
  other.optimizer->clone(optimizer);
  other.f->clone(f);

  if (other.fGradient.get() != nullptr) {
    other.fGradient->clone(fGradient);
  }

  if (fHessian.get() != nullptr) {
    other.fHessian->clone(fHessian);
  }
}

void FuzzyExtensionPrinciple::apply(const std::vector<const FuzzyInterval*>& x,
                                    std::unique_ptr<FuzzyInterval>& y) const {
  const size_t d = f->getNumberOfParameters();
  ScaledScalarFunction fScaled(*f);
  std::unique_ptr<ScaledScalarFunctionGradient> fGradientScaled;
  std::unique_ptr<ScaledScalarFunctionHessian> fHessianScaled;

  if (fGradient.get() != nullptr) {
    fGradientScaled.reset(new ScaledScalarFunctionGradient(*fGradient));
  }

  if (fHessian.get() != nullptr) {
    fHessianScaled.reset(new ScaledScalarFunctionHessian(*fHessian));
  }

  base::DataVector& lowerBounds(fScaled.getLowerBounds());
  base::DataVector& upperBounds(fScaled.getUpperBounds());

  base::DataVector xData(2 * m + 2);
  base::DataVector alphaData(2 * m + 2);

  for (size_t j = 0; j <= m; j++) {
    const double alpha = static_cast<double>(j) / static_cast<double>(m);

    for (size_t t = 0; t < d; t++) {
      lowerBounds[t] = x[t]->evaluateConfidenceIntervalLowerBound(alpha);
      upperBounds[t] = x[t]->evaluateConfidenceIntervalUpperBound(alpha);
    }

    if (fGradientScaled.get() != nullptr) {
      fGradientScaled->getLowerBounds() = lowerBounds;
      fGradientScaled->getUpperBounds() = upperBounds;
    }

    if (fHessianScaled.get() != nullptr) {
      fHessianScaled->getLowerBounds() = lowerBounds;
      fHessianScaled->getUpperBounds() = upperBounds;
    }

    {
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

      optimizer->optimize();

      xData[j] = optimizer->getOptimalValue();
      alphaData[j] = alpha;
    }

    {
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

      optimizer->optimize();

      xData[2*m+1-j] = -optimizer->getOptimalValue();
      alphaData[2*m+1-j] = alpha;
    }
  }

  y.reset(new InterpolatedFuzzyInterval(xData, alphaData));
}

size_t FuzzyExtensionPrinciple::getNumberOfAlphaSegments() const {
  return m;
}

void FuzzyExtensionPrinciple::setNumberOfAlphaSegments(size_t numberOfAlphaSegments) {
  m = numberOfAlphaSegments;
}

}  // namespace optimization
}  // namespace sgpp
