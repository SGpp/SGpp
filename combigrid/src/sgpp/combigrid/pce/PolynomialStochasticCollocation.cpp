// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/combigrid/functions/AbstractInfiniteFunctionBasis1D.hpp>
#include <sgpp/combigrid/functions/OrthogonalPolynomialBasis1D.hpp>
#include <sgpp/combigrid/integration/GaussLegendreQuadrature.hpp>
#include <sgpp/combigrid/operation/CombigridMultiOperation.hpp>
#include <sgpp/combigrid/operation/CombigridOperation.hpp>
#include <sgpp/combigrid/operation/CombigridTensorOperation.hpp>

#include <sgpp/combigrid/operation/multidim/AveragingLevelManager.hpp>

#include <sgpp/combigrid/pce/CombigridSurrogateModel.hpp>
#include <sgpp/combigrid/pce/PolynomialStochasticCollocation.hpp>
#include <sgpp/combigrid/pce/SGppToDakota.hpp>

#include <sgpp/combigrid/algebraic/FirstMomentNormStrategy.hpp>
#include <sgpp/combigrid/algebraic/VarianceNormStrategy.hpp>

#include <sgpp/base/exception/application_exception.hpp>
#include <sgpp/base/tools/GaussLegendreQuadRule1D.hpp>

#include <algorithm>
#include <cmath>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace sgpp {
namespace combigrid {

PolynomialStochasticCollocation::PolynomialStochasticCollocation(
    sgpp::combigrid::CombigridSurrogateModelConfiguration& config)
    : CombigridSurrogateModel(config),
      weightFunctions(0),
      computedMeanFlag(false),
      ev(0.0),
      computedVarianceFlag(false),
      var(0.0) {
  // create vector of function bases
  if (config.weightFunctions.size() == 0) {
    if (!config.weightFunction) {
      throw sgpp::base::application_exception(
          "PolynomialStochasticCollocation: no weight function available");
    }
    sgpp::combigrid::SingleFunction weightFunction = *config.weightFunction.get();
    for (size_t idim = 0; idim < numDims; idim++) {
      weightFunctions.push_back(weightFunction);
    }
  } else if (numDims == config.weightFunctions.size()) {
    for (size_t idim = 0; idim < numDims; idim++) {
      weightFunctions.push_back(config.weightFunctions[idim]);
    }
  } else {
    throw sgpp::base::application_exception(
        "PolynomialStochasticCollocation: number of weight functions do not match with the number "
        "of dimensions of the operation");
  }

  // initialize legendre basis
  sgpp::combigrid::OrthogonalPolynomialBasis1DConfiguration basisConfig;
  basisConfig.polyParameters.type_ = sgpp::combigrid::OrthogonalPolynomialBasisType::LEGENDRE;
  basisConfig.polyParameters.lowerBound_ = 0.0;
  basisConfig.polyParameters.upperBound_ = 1.0;
  legendreBasis = std::make_shared<sgpp::combigrid::OrthogonalPolynomialBasis1D>(basisConfig);

  updateConfig(config);

  initializeBounds();
  initializeNormStrategies();
}

PolynomialStochasticCollocation::~PolynomialStochasticCollocation() {}

// --------------------------------------------------------------------------------------

void PolynomialStochasticCollocation::initializeBounds() {
  if (config.bounds.size() == 0) {
    config.bounds.resize(2 * numDims);
    for (size_t idim = 0; idim < numDims; idim++) {
      config.bounds[2 * idim] = 0.0;
      config.bounds[2 * idim + 1] = 1.0;
    }
  } else {
    if (config.bounds.size() != 2 * numDims) {
      throw sgpp::base::application_exception(
          "PolynomialStochasticCollocation::initializeBounds - not enough arguments for bounds "
          "specified");
    }
  }
}

void PolynomialStochasticCollocation::initializeNormStrategies() {
  firstMomentNormstrategy.reset(
      new FirstMomentNormStrategy(legendreBasis, weightFunctions, false, config.bounds));
  varianceNormStrategy.reset(
      new VarianceNormStrategy(legendreBasis, weightFunctions, false, config.bounds));
}

double PolynomialStochasticCollocation::eval(sgpp::base::DataVector& x) {
  double ans = 0.0;
  for (auto it = expansionCoefficients.getValues()->getStoredDataIterator(); it->isValid();
       it->moveToNext()) {
    MultiIndex ix = it->getMultiIndex();
    double coeff = it->value().value();

    // evaluate tensor product
    double poly = 1.0;
    for (size_t k = 0; k < ix.size(); k++) {
      poly *= legendreBasis->evaluate(ix[k], x[k]);
    }
    ans += coeff * poly;
  }
  return ans;
}

void PolynomialStochasticCollocation::eval(sgpp::base::DataMatrix& xs,
                                           sgpp::base::DataVector& res) {
  size_t numSamples = xs.getNrows();
  res.resize(numSamples);
  res.setAll(0.0);
  for (auto it = expansionCoefficients.getValues()->getStoredDataIterator(); it->isValid();
       it->moveToNext()) {
    MultiIndex ix = it->getMultiIndex();
    double coeff = it->value().value();

    for (size_t i = 0; i < numSamples; i++) {
      // evaluate tensor product
      double poly = 1.0;
      for (size_t k = 0; k < ix.size(); k++) {
        poly *= legendreBasis->evaluate(ix[k], xs.get(i, k));
      }

      res[i] += coeff * poly;
    }
  }
}

double PolynomialStochasticCollocation::computeMean() {
  return firstMomentNormstrategy->norm(expansionCoefficients);
}

double PolynomialStochasticCollocation::mean() {
  if (!computedMeanFlag) {
    ev = computeMean();
    computedMeanFlag = true;
  }
  return ev;
}

double PolynomialStochasticCollocation::computeVariance() {
  return varianceNormStrategy->norm(expansionCoefficients);
}

double PolynomialStochasticCollocation::variance() {
  if (!computedVarianceFlag) {
    var = computeVariance();
    computedVarianceFlag = true;
  }
  return var;
}

void PolynomialStochasticCollocation::getComponentSobolIndices(
    sgpp::base::DataVector& componentSsobolIndices, bool normalized) {
  throw sgpp::base::application_exception(
      "PolynomialStochasticCollocation::getComponentSobolIndices - not implemented.");
}

void PolynomialStochasticCollocation::getTotalSobolIndices(
    sgpp::base::DataVector& totalSobolIndices, bool normalized) {
  throw sgpp::base::application_exception(
      "PolynomialStochasticCollocation::getTotalSobolIndices - not implemented.");
}

void PolynomialStochasticCollocation::updateConfig(
    sgpp::combigrid::CombigridSurrogateModelConfiguration config) {
  // update the tensor operation
  if (config.pointHierarchies.size() == numDims && config.storage) {
    combigridTensorOperation =
        sgpp::combigrid::CombigridTensorOperation::createOperationTensorPolynomialInterpolation(
            config.pointHierarchies, config.storage, legendreBasis);
  }

  if (!combigridTensorOperation) {
    throw sgpp::base::application_exception(
        "PolynomialStochasticCollocation:updateConfig - tensor operation is null -> something is "
        "seriously wrong");
  }

  if (config.levelManager) {
    combigridTensorOperation->setLevelManager(config.levelManager);
  }

  if (config.levelStructure) {
    combigridTensorOperation->getLevelManager()->addLevelsFromStructure(config.levelStructure);

    // update status
    computedMeanFlag = false;
    computedVarianceFlag = false;
  }

  if (config.enableLevelManagerStatsCollection) {
    combigridTensorOperation->getLevelManager()->enableStatsCollection();
  }

  expansionCoefficients = combigridTensorOperation->getResult();
}

size_t PolynomialStochasticCollocation::numGridPoints() {
  return combigridTensorOperation->numGridPoints();
}

std::shared_ptr<LevelInfos> PolynomialStochasticCollocation::getInfoOnAddedLevels() {
  return combigridTensorOperation->getLevelManager()->getInfoOnAddedLevels();
}

} /* namespace combigrid */
} /* namespace sgpp */
