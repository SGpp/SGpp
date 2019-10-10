// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/globaldef.hpp>

#include <sgpp/base/tools/Printer.hpp>
#include <sgpp/base/tools/RandomNumberGenerator.hpp>
#include <sgpp/optimization/optimizer/unconstrained/MultiStart.hpp>

#include <algorithm>
#include <cstdlib>
#include <limits>
#include <string>
#include <vector>

namespace sgpp {
namespace optimization {
namespace optimizer {

MultiStart::MultiStart(const base::ScalarFunction& f, size_t maxFcnEvalCount, size_t populationSize)
    : UnconstrainedOptimizer(f, nullptr, nullptr, maxFcnEvalCount),
      defaultOptimizer(NelderMead(f)) {
  defaultOptimizer.clone(optimizer);
  initialize(populationSize);
}

MultiStart::MultiStart(const UnconstrainedOptimizer& optimizer, size_t maxFcnEvalCount,
                       size_t populationSize)
    : UnconstrainedOptimizer(optimizer.getObjectiveFunction(),
                             optimizer.getObjectiveGradient(),
                             optimizer.getObjectiveHessian(),
                             maxFcnEvalCount),
      defaultOptimizer(NelderMead(*f)) {
  optimizer.clone(this->optimizer);
  initialize(populationSize);
}

MultiStart::MultiStart(const MultiStart& other)
    : UnconstrainedOptimizer(other),
      defaultOptimizer(NelderMead(*f)) {
  other.optimizer->clone(optimizer);
  initialize(other.populationSize);
}

MultiStart::~MultiStart() {}

void MultiStart::initialize(size_t populationSize) {
  this->populationSize = (populationSize > 0)
                             ? populationSize
                             : std::min(10 * f->getNumberOfParameters(), static_cast<size_t>(100));
}

void MultiStart::setObjectiveFunction(const base::ScalarFunction& f) {
  UnconstrainedOptimizer::setObjectiveFunction(f);
  optimizer->setObjectiveFunction(f);
}

void MultiStart::setObjectiveGradient(const base::ScalarFunctionGradient* fGradient) {
  UnconstrainedOptimizer::setObjectiveGradient(fGradient);
  optimizer->setObjectiveGradient(fGradient);
}

void MultiStart::setObjectiveHessian(const base::ScalarFunctionHessian* fHessian) {
  UnconstrainedOptimizer::setObjectiveHessian(fHessian);
  optimizer->setObjectiveHessian(fHessian);
}

void MultiStart::optimize() {
  base::Printer::getInstance().printStatusBegin("Optimizing (multi-start)...");

  const size_t d = f->getNumberOfParameters();

  xOpt.resize(0);
  fOpt = std::numeric_limits<double>::quiet_NaN();
  xHist.resize(0, d);
  fHist.resize(0);
  kHist.clear();

  std::vector<base::DataVector> x0(populationSize, base::DataVector(d));
  std::vector<size_t> roundN(populationSize, 0);
  size_t remainingN = N;

  // split the number of function evaluations evenly up for all points,
  // generate pseudorandom starting points
  for (size_t k = 0; k < populationSize; k++) {
    roundN[k] = static_cast<size_t>(
        std::ceil(static_cast<double>(remainingN) / static_cast<double>(populationSize - k)));
    remainingN -= roundN[k];

    for (size_t t = 0; t < d; t++) {
      x0[k][t] = base::RandomNumberGenerator::getInstance().getUniformRN();
    }
  }

  base::DataVector xCurrentOpt(d);
  double fCurrentOpt = std::numeric_limits<double>::infinity();

  // temporarily save x0 and N (will be overwritten by the loop)
  const base::DataVector tmpX0(optimizer->getStartingPoint());
  const size_t tmpN = optimizer->getN();
  const bool statusPrintingEnabled = base::Printer::getInstance().isStatusPrintingEnabled();

  if (statusPrintingEnabled) {
    base::Printer::getInstance().disableStatusPrinting();
  }

  size_t pointsDone = 0;

#pragma omp parallel shared(x0, roundN, xCurrentOpt, fCurrentOpt, pointsDone)
  {
    UnconstrainedOptimizer* curOptimizerPtr = optimizer.get();
#ifdef _OPENMP
    std::unique_ptr<UnconstrainedOptimizer> curOptimizer;

    if (omp_get_max_threads() > 1) {
      optimizer->clone(curOptimizer);
      curOptimizerPtr = curOptimizer.get();
    }

#endif /* _OPENMP */

    base::DataVector xLocalOpt(d);
    double fLocalOpt;

#pragma omp for ordered schedule(dynamic)

    for (size_t k = 0; k < populationSize; k++) {
      // optimize with k-th starting point
      curOptimizerPtr->setStartingPoint(x0[k]);
      curOptimizerPtr->setN(roundN[k]);
      curOptimizerPtr->optimize();

      xLocalOpt = curOptimizerPtr->getOptimalPoint();
      fLocalOpt = curOptimizerPtr->getOptimalValue();

#pragma omp critical
      {
        if (fLocalOpt < fCurrentOpt) {
          // this point is the best so far
          xCurrentOpt = xLocalOpt;
          fCurrentOpt = fLocalOpt;
        }
      }

#pragma omp atomic
      pointsDone++;

// status printing
      if (statusPrintingEnabled) {
        char str[10];
        snprintf(str, sizeof(str), "%.1f%%",
                 static_cast<double>(k) / static_cast<double>(populationSize) * 100.0);
        base::Printer::getInstance().getMutex().lock();
        base::Printer::getInstance().enableStatusPrinting();
        base::Printer::getInstance().printStatusUpdate(std::string(str) +
                                                       ", f(x) = " + std::to_string(fCurrentOpt));
        base::Printer::getInstance().disableStatusPrinting();
        base::Printer::getInstance().getMutex().unlock();
      }

#pragma omp critical
      {
        xHist.appendRow(xCurrentOpt);
        fHist.append(fCurrentOpt);
        kHist.push_back(curOptimizerPtr->getHistoryOfOptimalPoints().getNrows());
      }
    }
  }

  // set x0 and N to initial values
  optimizer->setStartingPoint(tmpX0);
  optimizer->setN(tmpN);

  xOpt.resize(d);
  xOpt = xCurrentOpt;
  fOpt = fCurrentOpt;

  if (statusPrintingEnabled) {
    base::Printer::getInstance().enableStatusPrinting();
  }

  base::Printer::getInstance().printStatusUpdate("100.0%, f(x) = " + std::to_string(fOpt));
  base::Printer::getInstance().printStatusEnd();
}

size_t MultiStart::getPopulationSize() const { return populationSize; }

void MultiStart::setPopulationSize(size_t populationSize) { this->populationSize = populationSize; }

const std::vector<size_t>& MultiStart::getHistoryOfInnerIterations() const { return kHist; }

void MultiStart::clone(std::unique_ptr<UnconstrainedOptimizer>& clone) const {
  clone = std::unique_ptr<UnconstrainedOptimizer>(new MultiStart(*this));
}
}  // namespace optimizer
}  // namespace optimization
}  // namespace sgpp
