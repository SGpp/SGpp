// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/combigrid/operation/multidim/fullgrid/AbstractBasisCoefficientsStorage.hpp>
#include <sgpp/combigrid/storage/AbstractMultiStorageIterator.hpp>

#include <vector>
#include <memory>

namespace sgpp {
namespace combigrid {

AbstractBasisCoefficientsStorage::AbstractBasisCoefficientsStorage() {}

AbstractBasisCoefficientsStorage::~AbstractBasisCoefficientsStorage() {}

std::shared_ptr<std::vector<double>> AbstractBasisCoefficientsStorage::getCoefficients(
    MultiIndex const& level, std::shared_ptr<AbstractCombigridStorage>& storage,
    MultiIndex& multiBounds, std::vector<bool>& orderingConfiguration) {
  if (coefficients.find(level) == coefficients.end()) {
    computeCoefficients(level, storage, multiBounds, orderingConfiguration);
  }
  return coefficients[level];
}

} /* namespace combigrid */
} /* namespace sgpp */
