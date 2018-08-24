// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/globaldef.hpp>

#include <sgpp/optimization/fuzzy/FuzzyExtensionPrinciple.hpp>
#include <sgpp/optimization/fuzzy/InterpolatedFuzzyInterval.hpp>
#include <sgpp/optimization/tools/Printer.hpp>

#include <limits>
#include <string>
#include <vector>

namespace sgpp {
namespace optimization {

FuzzyExtensionPrinciple::FuzzyExtensionPrinciple(
    const ScalarFunction& f,
    size_t numberOfAlphaSegments) :
      m(numberOfAlphaSegments) {
  f.clone(this->f);
}

FuzzyExtensionPrinciple::FuzzyExtensionPrinciple(const FuzzyExtensionPrinciple& other) :
  m(other.m),
  alphaLevels(other.alphaLevels),
  optimizationDomainsLowerBounds(other.optimizationDomainsLowerBounds),
  optimizationDomainsUpperBounds(other.optimizationDomainsUpperBounds),
  minimumPoints(other.minimumPoints),
  minimumValues(other.minimumValues),
  maximumPoints(other.maximumPoints),
  maximumValues(other.maximumValues) {
  other.f->clone(f);
}

FuzzyExtensionPrinciple::~FuzzyExtensionPrinciple() {}

void FuzzyExtensionPrinciple::prepareApply() {}

FuzzyInterval* FuzzyExtensionPrinciple::apply(
    const std::vector<const FuzzyInterval*>& xFuzzy) {
  Printer::getInstance().printStatusBegin("Applying fuzzy extension principle...");

  const size_t d = f->getNumberOfParameters();

  // result data
  base::DataVector xData(2 * m + 2);
  base::DataVector alphaData(2 * m + 2);
  alphaLevels.resize(m + 1);
  optimizationDomainsLowerBounds.resize(m + 1, base::DataVector(d));
  optimizationDomainsUpperBounds.resize(m + 1, base::DataVector(d));
  minimumPoints.resize(m + 1, base::DataVector(d));
  minimumValues.resize(m + 1);
  maximumPoints.resize(m + 1, base::DataVector(d));
  maximumValues.resize(m + 1);

  Printer::getInstance().printStatusUpdate("calculating confidence intervals");
  Printer::getInstance().printStatusNewLine();

  // calculate alpha levels and corresponding confidence intervals
  for (size_t j = 0; j <= m; j++) {
    alphaLevels[j] = static_cast<double>(j) / static_cast<double>(m);

    // determine input parameter confidence interval,
    // directly changing the optimization domain in fScaled
    for (size_t t = 0; t < d; t++) {
      optimizationDomainsLowerBounds[j][t] =
          xFuzzy[t]->evaluateConfidenceIntervalLowerBound(alphaLevels[j]);
      optimizationDomainsUpperBounds[j][t] =
          xFuzzy[t]->evaluateConfidenceIntervalUpperBound(alphaLevels[j]);
    }
  }

  // current optimal points/values
  std::vector<std::unique_ptr<base::DataVector>> curMinimumPoints;
  base::DataVector curMinimumValues(m + 1);
  std::vector<std::unique_ptr<base::DataVector>> curMaximumPoints;
  base::DataVector curMaximumValues(m + 1);

  for (size_t j = 0; j <= m; j++) {
    curMinimumPoints.emplace_back(new base::DataVector(d));
    curMaximumPoints.emplace_back(new base::DataVector(d));
  }

  // save last optimal min/max value and corresponding argmin/argmax point
  // to make sure that the optimum for a smaller alpha (hence for a larger optimization domain)
  // is not worse that for a larger alpha
  double lastOptimalValueMin = std::numeric_limits<double>::infinity();
  double lastOptimalValueMax = -std::numeric_limits<double>::infinity();
  base::DataVector lastOptimalPointMin(d);
  base::DataVector lastOptimalPointMax(d);

  const bool statusPrintingEnabled = Printer::getInstance().isStatusPrintingEnabled();

  if (statusPrintingEnabled) {
    Printer::getInstance().disableStatusPrinting();
  }

  size_t alphaLevelsDone = 0;

#pragma omp parallel shared(curMinimumPoints, curMinimumValues, \
    alphaLevelsDone, curMaximumPoints, curMaximumValues) default(none)
  {
    std::unique_ptr<FuzzyExtensionPrinciple> curFuzzyExtensionPrinciple;
    clone(curFuzzyExtensionPrinciple);

    // call custom preparation method (may be required and implemented by sub-class)
    curFuzzyExtensionPrinciple->prepareApply();

#pragma omp for schedule(static)
    for (size_t j = 0; j <= m; j++) {
      curFuzzyExtensionPrinciple->optimizeForSingleAlphaLevel(
          j, *curMinimumPoints[j], curMinimumValues[j],
          *curMaximumPoints[j], curMaximumValues[j]);

#pragma omp atomic
      alphaLevelsDone++;

      // status message
      if (statusPrintingEnabled) {
        char str[10];
        snprintf(str, sizeof(str), "%.1f%%",
                 static_cast<double>(alphaLevelsDone) / static_cast<double>(m + 1) * 100.0);
        Printer::getInstance().getMutex().lock();
        Printer::getInstance().enableStatusPrinting();
        Printer::getInstance().printStatusUpdate("optimizing (" + std::string(str) +
                                                 ")");
        Printer::getInstance().disableStatusPrinting();
        Printer::getInstance().getMutex().unlock();
      }
    }
  }

  if (statusPrintingEnabled) {
    Printer::getInstance().enableStatusPrinting();
  }

  Printer::getInstance().printStatusUpdate("optimizing (100.0%)");
  Printer::getInstance().printStatusNewLine();

  // iterate through alphas from 1 to 0
  for (size_t j = m + 1; j-- > 0;) {
    // check if optimal value is indeed smaller than the last minimum
    if (curMinimumValues[j] < lastOptimalValueMin) {
      xData[j] = curMinimumValues[j];
      lastOptimalValueMin = curMinimumValues[j];
      lastOptimalPointMin = *curMinimumPoints[j];
    } else {
      xData[j] = lastOptimalValueMin;
    }

    minimumPoints[j] = lastOptimalPointMin;
    minimumValues[j] = lastOptimalValueMin;
    alphaData[j] = alphaLevels[j];

    // check if optimal value is indeed larger than the last maximum
    if (curMaximumValues[j] > lastOptimalValueMax) {
      xData[2*m+1-j] = curMaximumValues[j];
      lastOptimalValueMax = curMaximumValues[j];
      lastOptimalPointMax = *curMaximumPoints[j];
    } else {
      xData[2*m+1-j] = lastOptimalValueMax;
    }

    maximumPoints[j] = lastOptimalPointMax;
    maximumValues[j] = lastOptimalValueMax;
    alphaData[2*m+1-j] = alphaLevels[j];
  }

  Printer::getInstance().printStatusEnd();

  // interpolate between alpha data points
  return new InterpolatedFuzzyInterval(xData, alphaData);
}

size_t FuzzyExtensionPrinciple::getNumberOfAlphaSegments() const {
  return m;
}

void FuzzyExtensionPrinciple::setNumberOfAlphaSegments(
    size_t numberOfAlphaSegments) {
  m = numberOfAlphaSegments;
}

const base::DataVector& FuzzyExtensionPrinciple::getAlphaLevels() const {
  return alphaLevels;
}

const std::vector<base::DataVector>&
FuzzyExtensionPrinciple::getOptimizationDomainsLowerBounds() const {
  return optimizationDomainsLowerBounds;
}

const std::vector<base::DataVector>&
FuzzyExtensionPrinciple::getOptimizationDomainsUpperBounds() const {
  return optimizationDomainsUpperBounds;
}

const std::vector<base::DataVector>& FuzzyExtensionPrinciple::getMinimumPoints() const {
  return minimumPoints;
}

const base::DataVector& FuzzyExtensionPrinciple::getMinimumValues() const {
  return minimumValues;
}

const std::vector<base::DataVector>& FuzzyExtensionPrinciple::getMaximumPoints() const {
  return maximumPoints;
}

const base::DataVector& FuzzyExtensionPrinciple::getMaximumValues() const {
  return maximumValues;
}

}  // namespace optimization
}  // namespace sgpp
