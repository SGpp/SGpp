// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/combigrid/algebraic/FloatScalarVector.hpp>
#include <sgpp/combigrid/storage/AbstractCombigridStorage.hpp>
#include <sgpp/combigrid/storage/AbstractMultiStorage.hpp>

#include <map>
#include <memory>
#include <vector>

namespace sgpp {
namespace combigrid {

template <class V>
class AbstractBasisCoefficientsStorage {
 public:
  AbstractBasisCoefficientsStorage();

  virtual ~AbstractBasisCoefficientsStorage();
  virtual void computeCoefficients(MultiIndex const& level,
                                   std::shared_ptr<AbstractCombigridStorage>& storage,
                                   MultiIndex& multiBounds,
                                   std::vector<bool>& orderingConfiguration,
                                   std::vector<std::vector<V>> const& basisValues) = 0;
  std::shared_ptr<std::vector<double>> getCoefficients(
      MultiIndex const& level, std::shared_ptr<AbstractCombigridStorage>& storage,
      MultiIndex& multiBounds, std::vector<bool>& orderingConfiguration,
      std::vector<std::vector<V>> const& basisValues);

 protected:
  std::map<MultiIndex, std::shared_ptr<std::vector<double>>> coefficients;
  base::DataVector functionValues;
};

template <class V>
AbstractBasisCoefficientsStorage<V>::AbstractBasisCoefficientsStorage() {}

template <class V>
AbstractBasisCoefficientsStorage<V>::~AbstractBasisCoefficientsStorage() {}

template <class V>
std::shared_ptr<std::vector<double>> AbstractBasisCoefficientsStorage<V>::getCoefficients(
    MultiIndex const& level, std::shared_ptr<AbstractCombigridStorage>& storage,
    MultiIndex& multiBounds, std::vector<bool>& orderingConfiguration,
    std::vector<std::vector<V>> const& basisValues) {
  if (coefficients.find(level) == coefficients.end()) {
    computeCoefficients(level, storage, multiBounds, orderingConfiguration, basisValues);
  }
  return coefficients[level];
}

} /* namespace combigrid */
} /* namespace sgpp */
