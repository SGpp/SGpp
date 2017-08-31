// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/combigrid/operation/multidim/fullgrid/FunctionValuesCoefficientsStorage.hpp>
#include <sgpp/combigrid/operation/multidim/fullgrid/AbstractBasisCoefficientsStorage.hpp>
#include <sgpp/combigrid/storage/AbstractMultiStorageIterator.hpp>

#include <vector>
#include <memory>

namespace sgpp {
namespace combigrid {

FunctionValuesCoefficientsStorage::FunctionValuesCoefficientsStorage()
    : AbstractBasisCoefficientsStorage() {}

FunctionValuesCoefficientsStorage::~FunctionValuesCoefficientsStorage() {}

void FunctionValuesCoefficientsStorage::computeCoefficients(
    MultiIndex const& level, std::shared_ptr<AbstractCombigridStorage>& storage,
    MultiIndex& multiBounds, std::vector<bool>& orderingConfiguration) {
  MultiIndexIterator it(multiBounds);
  auto funcIter = storage->getGuidedIterator(level, it, orderingConfiguration);
  coefficients[level] = std::make_shared<std::vector<double>>();
  while (true) {
    // get function value and partial product and multiply them together with the last basis
    // coefficient, then add the resulting value to the total sum
    coefficients[level]->push_back(funcIter->value());

    // increment iterator
    int h = funcIter->moveToNext();

    if (h < 0) {
      break;  // all indices have been traversed, stop iteration
    }
  }
}

} /* namespace combigrid */
} /* namespace sgpp */
