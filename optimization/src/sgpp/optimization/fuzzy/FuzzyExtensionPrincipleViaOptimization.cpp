// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/globaldef.hpp>

#include <sgpp/optimization/fuzzy/FuzzyExtensionPrincipleViaOptimization.hpp>

#include <limits>
#include <vector>

namespace sgpp {
namespace optimization {

FuzzyExtensionPrincipleViaOptimization::FuzzyExtensionPrincipleViaOptimization(
    const ScalarFunction& f,
    size_t numberOfAlphaSegments) :
        FuzzyExtensionPrinciple(f, numberOfAlphaSegments),
        defaultOptimizer(optimizer::MultiStart(f)) {
  defaultOptimizer.clone(optimizer);
}

FuzzyExtensionPrincipleViaOptimization::FuzzyExtensionPrincipleViaOptimization(
    const optimizer::UnconstrainedOptimizer& optimizer,
    size_t numberOfAlphaSegments) :
      FuzzyExtensionPrinciple(optimizer.getObjectiveFunction(), numberOfAlphaSegments),
      defaultOptimizer(optimizer::MultiStart(optimizer.getObjectiveFunction())) {
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
    FuzzyExtensionPrinciple(*other.f, other.m),
    defaultOptimizer(optimizer::MultiStart(*other.f)) {
  other.optimizer->clone(optimizer);

  if (other.fGradient.get() != nullptr) {
    other.fGradient->clone(fGradient);
  }

  if (fHessian.get() != nullptr) {
    other.fHessian->clone(fHessian);
  }
}

FuzzyExtensionPrincipleViaOptimization::~FuzzyExtensionPrincipleViaOptimization() {}

void FuzzyExtensionPrincipleViaOptimization::prepareApply() {
  fScaled.reset(new ScaledScalarFunction(*f));

  // create scaled gradient function if gradient is given
  if (fGradient.get() != nullptr) {
    fGradientScaled.reset(new ScaledScalarFunctionGradient(*fGradient));
  }

  // create scaled Hessian function if Hessian is given
  if (fHessian.get() != nullptr) {
    fHessianScaled.reset(new ScaledScalarFunctionHessian(*fHessian));
  }
}

void FuzzyExtensionPrincipleViaOptimization::optimizeForSingleAlphaLevel(
    size_t j, base::DataVector& minimumPoint, double& minimumValue,
    base::DataVector& maximumPoint, double& maximumValue) {
  const size_t d = f->getNumberOfParameters();
  const base::DataVector& lowerBounds = optimizationDomainsLowerBounds[j];
  const base::DataVector& upperBounds = optimizationDomainsUpperBounds[j];

  fScaled->setLowerBounds(lowerBounds);
  fScaled->setUpperBounds(upperBounds);

  if (fGradientScaled.get() != nullptr) {
    fGradientScaled->setLowerBounds(lowerBounds);
    fGradientScaled->setUpperBounds(upperBounds);
  }

  if (fHessianScaled.get() != nullptr) {
    fHessianScaled->setLowerBounds(lowerBounds);
    fHessianScaled->setUpperBounds(upperBounds);
  }

  // compute minimum (lower bound of output confidence interval)
  {
    // value factor of 1 means minimization
    fScaled->setValueFactor(1.0);
    optimizer->setObjectiveFunction(*fScaled);

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
      optimizer->setStartingPoint(minimumPoints[j + 1]);
    }

    optimizer->optimize();

    // save minimum points
    minimumPoint = optimizer->getOptimalPoint();
    minimumValue = optimizer->getOptimalValue();

    for (size_t t = 0; t < d; t++) {
      minimumPoint[t] = lowerBounds[t] + minimumPoint[t] * (upperBounds[t] - lowerBounds[t]);
    }
  }

  // compute maximum (upper bound of output confidence interval)
  {
    // value factor of -1 means maximization
    fScaled->setValueFactor(-1.0);
    optimizer->setObjectiveFunction(*fScaled);

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
      optimizer->setStartingPoint(maximumPoints[j + 1]);
    }

    optimizer->optimize();

    // save minimum points
    // (don't forget the minus because the value factor is not incorporated in getOptimalValue)
    maximumPoint = optimizer->getOptimalPoint();
    maximumValue = -optimizer->getOptimalValue();

    for (size_t t = 0; t < d; t++) {
      maximumPoint[t] = lowerBounds[t] + maximumPoint[t] * (upperBounds[t] - lowerBounds[t]);
    }
  }
}

void FuzzyExtensionPrincipleViaOptimization::clone(
    std::unique_ptr<FuzzyExtensionPrinciple>& clone) const {
  clone = std::unique_ptr<FuzzyExtensionPrinciple>(
      new FuzzyExtensionPrincipleViaOptimization(*this));
}

}  // namespace optimization
}  // namespace sgpp
