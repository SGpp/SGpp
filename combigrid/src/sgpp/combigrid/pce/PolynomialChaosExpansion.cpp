// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/combigrid/functions/AbstractInfiniteFunctionBasis1D.hpp>
#include <sgpp/combigrid/functions/OrthogonalBasisFunctionsCollection.hpp>
#include <sgpp/combigrid/operation/CombigridMultiOperation.hpp>
#include <sgpp/combigrid/operation/CombigridOperation.hpp>
#include <sgpp/combigrid/operation/CombigridTensorOperation.hpp>

#include <sgpp/combigrid/operation/multidim/AveragingLevelManager.hpp>

#include <sgpp/combigrid/pce/PolynomialChaosExpansion.hpp>
#include <sgpp/combigrid/pce/SGppToDakota.hpp>

#include <sgpp/base/exception/application_exception.hpp>

#ifdef USE_DAKOTA
#include <pecos_data_types.hpp>
#endif

#include <vector>

namespace sgpp {
namespace combigrid {

PolynomialChaosExpansion::PolynomialChaosExpansion(
    sgpp::combigrid::CombigridSurrogateModelConfiguration& config)
    : CombigridSurrogateModel(config), basisFunctions(0), computedSobolIndicesFlag(false) {
  if (config.basisFunctions.size() == 0 && config.basisFunction) {
    for (size_t idim = 0; idim < numDims; idim++) {
      basisFunctions.push_back(config.basisFunction);
    }
  } else if (config.basisFunctions.size() == numDims) {
    for (size_t idim = 0; idim < numDims; idim++) {
      basisFunctions.push_back(config.basisFunctions[idim]);
    }
  } else {
    throw sgpp::base::application_exception(
        "PolynomialChaosExpansion: number of basis function do not match with the number of "
        "dimensions of the operation");
  }

  updateConfig(config);
}

PolynomialChaosExpansion::~PolynomialChaosExpansion() {}

double PolynomialChaosExpansion::eval(sgpp::base::DataVector& x) {
  double ans = 0.0;
  for (auto it = expansionCoefficients.getValues()->getStoredDataIterator(); it->isValid();
       it->moveToNext()) {
    MultiIndex ix = it->getMultiIndex();
    double coeff = it->value().value();

    // evaluate tensor product
    double poly = 1.0;
    for (size_t k = 0; k < ix.size(); k++) {
      poly *= basisFunctions[k]->evaluate(ix[k], x[k]);
    }
    ans += coeff * poly;
  }
  return ans;
}

void PolynomialChaosExpansion::eval(sgpp::base::DataMatrix& xs, sgpp::base::DataVector& res) {
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
        poly *= basisFunctions[k]->evaluate(ix[k], xs.get(i, k));
      }

      res[i] += coeff * poly;
    }
  }
}

double PolynomialChaosExpansion::mean() {
  return expansionCoefficients.get(MultiIndex(numDims, 0)).getValue();
}

double PolynomialChaosExpansion::variance() {
  double var = 0.0;
  auto it = expansionCoefficients.getValues()->getStoredDataIterator();
  if (!it->isValid()) {
    return 0.0;
  }
  it->moveToNext();  // ignore first entry (belonging to mean)

  for (; it->isValid(); it->moveToNext()) {
    double coeff = it->value().value();
    var += coeff * coeff;
  }
  return var;
}

void PolynomialChaosExpansion::computeComponentSobolIndices() {
  // compute the component sobol indices
  if (computedSobolIndicesFlag) {
    return;
  }

  size_t numSobolIndices = static_cast<size_t>(std::pow(2, numDims) - 1);
  sobolIndices.resizeZero(numSobolIndices);

  // load index vectors
  auto it_coeffs = expansionCoefficients.getValues()->getStoredDataIterator();
  std::vector<MultiIndex> indexList;
  for (; it_coeffs->isValid(); it_coeffs->moveToNext()) {
    indexList.push_back(it_coeffs->getMultiIndex());
  }

  for (size_t i = 0; i < numSobolIndices; i++) {
    // loop over all the remaining basis functions
    for (std::vector<MultiIndex>::iterator it_ixlist = indexList.begin();
         it_ixlist != indexList.end();) {
      MultiIndex multiIndex = *it_ixlist;

      // check if all the current dimensions are set
      bool dimsAreSet = true;
      size_t idim = 0;
      while (dimsAreSet && idim < numDims) {
        // get kth bit from permuation mask
        size_t isSet = ((i + 1) & (1 << idim)) >> idim;  // in {0, 1}
        // check if the degree of the basis in the current dimension
        // is > 0 if set and equal to 0 if not set
        if (isSet == 0) {
          dimsAreSet &= multiIndex[idim] == 0;
        } else {
          dimsAreSet &= multiIndex[idim] > 0;
        }
        idim++;
      }

      // add the squared coefficient of the current term to the
      // sobol index matrix if all dimensions are present
      // in the current term
      if (dimsAreSet) {
        // load index and delete it from list since every term
        // can just contribute to one sobol index
        indexList.erase(it_ixlist);
        double coefficient = expansionCoefficients.get(multiIndex).value();
        sobolIndices[i] += coefficient * coefficient;
      } else {
        it_ixlist++;
      }
    }
  }

  // update the computed flag
  computedSobolIndicesFlag = true;
}

void PolynomialChaosExpansion::getComponentSobolIndices(
    sgpp::base::DataVector& componentSobolIndices, bool normalized) {
  computeComponentSobolIndices();
  // copy sobol indices to output vector
  componentSobolIndices.resize(sobolIndices.getSize());
  componentSobolIndices.copyFrom(sobolIndices);

  // divide all the entries by the variance to obtain the Sobol indices
  if (normalized) {
    double var = variance();
    if (var > 1e-14) {
      componentSobolIndices.mult(1. / var);
    }
  }
}

void PolynomialChaosExpansion::getTotalSobolIndices(sgpp::base::DataVector& totalSobolIndices,
                                                    bool normalized) {
  computeComponentSobolIndices();
  totalSobolIndices.resizeZero(numDims);
  for (size_t idim = 0; idim < numDims; idim++) {
    for (size_t iperm = 0; iperm < sobolIndices.getSize(); iperm++) {
      // check if the current dimension is set in the current key of the sobol index
      size_t isSet = ((iperm + 1) & (1 << idim)) >> idim;  // in {0, 1}
      if (isSet == 1) {
        totalSobolIndices[idim] += sobolIndices[iperm];
      }
    }
  }

  // divide all the entries by the variance to obtain the Sobol indices
  if (normalized) {
    double var = variance();
    if (var > 1e-14) {
      totalSobolIndices.mult(1. / variance());
    }
  }
}

void PolynomialChaosExpansion::updateConfig(
    sgpp::combigrid::CombigridSurrogateModelConfiguration config) {
  if (config.tensorOperation) {
    config.loadFromCombigridOperation(config.tensorOperation);
    combigridTensorOperation = config.tensorOperation;
  } else if (config.pointHierarchies.size() == numDims && config.storage) {
    combigridTensorOperation =
        sgpp::combigrid::CombigridTensorOperation::createOperationTensorPolynomialInterpolation(
            config.pointHierarchies, config.storage, basisFunctions);
    config.tensorOperation = combigridTensorOperation;
  }

  if (!combigridTensorOperation) {
    throw sgpp::base::application_exception(
        "PolynomialChaosExpansion:updateConfig - tensor operation is null -> something is "
        "seriously wrong");
  }

  if (config.levelManager) {
    combigridTensorOperation->setLevelManager(config.levelManager);
  }

  if (config.levelStructure) {
    combigridTensorOperation->getLevelManager()->addLevelsFromStructure(config.levelStructure);
    computedSobolIndicesFlag = false;
  }

  if (config.enableLevelManagerStatsCollection) {
    combigridTensorOperation->getLevelManager()->enableStatsCollection();
  }

  expansionCoefficients = combigridTensorOperation->getResult();
}

size_t PolynomialChaosExpansion::numGridPoints() {
  return combigridTensorOperation->numGridPoints();
}

std::shared_ptr<LevelInfos> PolynomialChaosExpansion::getInfoOnAddedLevels() {
  return combigridTensorOperation->getLevelManager()->getInfoOnAddedLevels();
}

} /* namespace combigrid */
} /* namespace sgpp */
