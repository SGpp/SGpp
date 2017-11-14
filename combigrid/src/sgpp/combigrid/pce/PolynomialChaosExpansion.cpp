// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/combigrid/pce/PolynomialChaosExpansion.hpp>
#include <sgpp/combigrid/pce/SGppToDakota.hpp>
#include <sgpp/combigrid/operation/CombigridOperation.hpp>
#include <sgpp/combigrid/operation/CombigridMultiOperation.hpp>
#include <sgpp/combigrid/operation/CombigridTensorOperation.hpp>
#include <sgpp/combigrid/functions/AbstractInfiniteFunctionBasis1D.hpp>
#include <pecos_data_types.hpp>

#include <vector>

namespace sgpp {
namespace combigrid {

PolynomialChaosExpansion::PolynomialChaosExpansion(
    std::shared_ptr<sgpp::combigrid::CombigridOperation> combigridOperation,
    std::shared_ptr<sgpp::combigrid::AbstractInfiniteFunctionBasis1D> functionBasis)
    : numDims(combigridOperation->numDims()),
      combigridOperation(combigridOperation),
      combigridMultiOperation(nullptr),
      combigridTensorOperation(nullptr),
      numGridPoints(combigridOperation->numGridPoints()) {
  // create tensor operation for pce transformation
  this->combigridTensorOperation =
      sgpp::combigrid::CombigridTensorOperation::createOperationTensorPolynomialInterpolation(
          combigridOperation->getPointHierarchies(), combigridOperation->getStorage(),
          combigridOperation->getLevelManager(), functionBasis);
}

PolynomialChaosExpansion::PolynomialChaosExpansion(
    std::shared_ptr<sgpp::combigrid::CombigridMultiOperation> combigridMultiOperation,
    std::shared_ptr<sgpp::combigrid::AbstractInfiniteFunctionBasis1D> functionBasis)
    : numDims(combigridMultiOperation->numDims()),
      combigridOperation(nullptr),
      combigridMultiOperation(combigridMultiOperation),
      combigridTensorOperation(nullptr),
      numGridPoints(combigridMultiOperation->numGridPoints()) {
  // create tensor operation for pce transformation
  this->combigridTensorOperation =
      sgpp::combigrid::CombigridTensorOperation::createOperationTensorPolynomialInterpolation(
          combigridMultiOperation->getPointHierarchies(), combigridMultiOperation->getStorage(),
          combigridMultiOperation->getLevelManager(), functionBasis);
}

PolynomialChaosExpansion::PolynomialChaosExpansion(
    std::shared_ptr<sgpp::combigrid::CombigridTensorOperation> combigridTensorOperation,
    std::shared_ptr<sgpp::combigrid::AbstractInfiniteFunctionBasis1D> functionBasis)
    : numDims(combigridTensorOperation->numDims()),
      combigridOperation(nullptr),
      combigridMultiOperation(nullptr),
      combigridTensorOperation(combigridTensorOperation),
      numGridPoints(combigridTensorOperation->numGridPoints()) {
  // create tensor operation for pce transformation
  this->combigridTensorOperation =
      CombigridTensorOperation::createOperationTensorPolynomialInterpolation(
          combigridTensorOperation->getPointHierarchies(), combigridTensorOperation->getStorage(),
          combigridTensorOperation->getLevelManager(), functionBasis);
}

PolynomialChaosExpansion::PolynomialChaosExpansion(
    std::shared_ptr<sgpp::combigrid::CombigridOperation> combigridOperation,
    std::vector<std::shared_ptr<sgpp::combigrid::AbstractInfiniteFunctionBasis1D>>& functionBases)
    : numDims(combigridMultiOperation->numDims()),
      combigridOperation(combigridOperation),
      combigridMultiOperation(nullptr),
      combigridTensorOperation(nullptr),
      numGridPoints(combigridMultiOperation->numGridPoints()) {
  // create tensor operation for pce transformation
  this->combigridTensorOperation =
      sgpp::combigrid::CombigridTensorOperation::createOperationTensorPolynomialInterpolation(
          combigridOperation->getPointHierarchies(), combigridOperation->getStorage(),
          combigridOperation->getLevelManager(), functionBases);
}

PolynomialChaosExpansion::PolynomialChaosExpansion(
    std::shared_ptr<sgpp::combigrid::CombigridMultiOperation> combigridMultiOperation,
    std::vector<std::shared_ptr<sgpp::combigrid::AbstractInfiniteFunctionBasis1D>>& functionBases)
    : numDims(combigridMultiOperation->numDims()),
      combigridOperation(nullptr),
      combigridMultiOperation(combigridMultiOperation),
      combigridTensorOperation(nullptr),
      numGridPoints(combigridMultiOperation->numGridPoints()) {
  // create tensor operation for pce transformation
  this->combigridTensorOperation =
      sgpp::combigrid::CombigridTensorOperation::createOperationTensorPolynomialInterpolation(
          combigridMultiOperation->getPointHierarchies(), combigridMultiOperation->getStorage(),
          combigridMultiOperation->getLevelManager(), functionBases);
}

PolynomialChaosExpansion::PolynomialChaosExpansion(
    std::shared_ptr<sgpp::combigrid::CombigridTensorOperation> combigridTensorOperation,
    std::vector<std::shared_ptr<sgpp::combigrid::AbstractInfiniteFunctionBasis1D>>& functionBases)
    : numDims(combigridTensorOperation->numDims()),
      combigridOperation(nullptr),
      combigridMultiOperation(nullptr),
      combigridTensorOperation(nullptr),
      numGridPoints(combigridTensorOperation->numGridPoints()) {
  // create tensor operation for pce transformation
  this->combigridTensorOperation =
      sgpp::combigrid::CombigridTensorOperation::createOperationTensorPolynomialInterpolation(
          combigridTensorOperation->getPointHierarchies(), combigridTensorOperation->getStorage(),
          combigridTensorOperation->getLevelManager(), functionBases);
}

PolynomialChaosExpansion::~PolynomialChaosExpansion() {}

double PolynomialChaosExpansion::mean() {
  //  if (orthogPoly != nullptr) {
  // #ifdef USE_DAKOTA
  //    return orthogPoly->mean();
  // #endif
  //  }
  if (numGridPoints < combigridTensorOperation->numGridPoints()) {
    expansionCoefficients = combigridTensorOperation->getResult();
    numGridPoints = combigridTensorOperation->numGridPoints();
  }
  return expansionCoefficients.get(MultiIndex(numDims, 0)).getValue();
}

double PolynomialChaosExpansion::variance() {
  //  if (orthogPoly != nullptr) {
  // #ifdef USE_DAKOTA
  //    return orthogPoly->variance();
  // #endif
  //  }
  // PCE norm, i.e. variance (if using appropriate normalized orthogonal polynomials)
  if (numGridPoints < combigridTensorOperation->numGridPoints()) {
    expansionCoefficients = combigridTensorOperation->getResult();
    numGridPoints = combigridTensorOperation->numGridPoints();
  }

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

void PolynomialChaosExpansion::getComponentSobolIndices(
    sgpp::base::DataVector& componentSobolIndices, bool normalized) {
  if (numGridPoints < combigridTensorOperation->numGridPoints()) {
    computeComponentSobolIndices();
    numGridPoints = combigridTensorOperation->numGridPoints();
  }
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

void PolynomialChaosExpansion::computeComponentSobolIndices() {
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
}

void PolynomialChaosExpansion::getTotalSobolIndices(sgpp::base::DataVector& totalSobolIndices,
                                                    bool normalized) {
  // compute the component sobol indices
  if (numGridPoints < combigridTensorOperation->numGridPoints()) {
    computeComponentSobolIndices();
    numGridPoints = combigridTensorOperation->numGridPoints();
  }

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

std::shared_ptr<sgpp::combigrid::CombigridTensorOperation>
PolynomialChaosExpansion::getCombigridTensorOperation() {
  return combigridTensorOperation;
}

} /* namespace combigrid */
} /* namespace sgpp */
