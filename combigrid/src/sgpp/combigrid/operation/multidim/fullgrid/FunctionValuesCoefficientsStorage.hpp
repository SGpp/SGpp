// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/combigrid/operation/multidim/fullgrid/AbstractBasisCoefficientsStorage.hpp>
#include <sgpp/combigrid/storage/AbstractCombigridStorage.hpp>
#include <sgpp/combigrid/storage/AbstractMultiStorage.hpp>

#include <memory>
#include <vector>

namespace sgpp {
namespace combigrid {

template <class V>
class FunctionValuesCoefficientsStorage : public AbstractBasisCoefficientsStorage<V> {
 public:
  FunctionValuesCoefficientsStorage();

  virtual ~FunctionValuesCoefficientsStorage();

  void computeCoefficients(MultiIndex const& level,
                           std::shared_ptr<AbstractCombigridStorage>& storage,
                           MultiIndex& multiBounds, std::vector<bool>& orderingConfiguration,
                           std::vector<std::vector<V>> const& basisValues) override;
};

template <class V>
FunctionValuesCoefficientsStorage<V>::FunctionValuesCoefficientsStorage()
    : AbstractBasisCoefficientsStorage<V>() {}

template <class V>
FunctionValuesCoefficientsStorage<V>::~FunctionValuesCoefficientsStorage() {}

template <class V>
void FunctionValuesCoefficientsStorage<V>::computeCoefficients(
    MultiIndex const& level, std::shared_ptr<AbstractCombigridStorage>& storage,
    MultiIndex& multiBounds, std::vector<bool>& orderingConfiguration,
    std::vector<std::vector<V>> const& basisValues) {
  MultiIndexIterator it(multiBounds);
  auto funcIter = storage->getGuidedIterator(level, it, orderingConfiguration);
  this->coefficients[level] = std::make_shared<std::vector<double>>();
  while (true) {
    // get function value and partial product and multiply them together with the last basis
    // coefficient, then add the resulting value to the total sum
    this->coefficients[level]->push_back(funcIter->value());

    // increment iterator
    int h = funcIter->moveToNext();

    if (h < 0) {
      break;  // all indices have been traversed, stop iteration
    }
  }
}

} /* namespace combigrid */
} /* namespace sgpp */
