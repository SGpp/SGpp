// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include "CombigridTreeStorage.hpp"

#include <sgpp/combigrid/serialization/FloatSerializationStrategy.hpp>
#include <sgpp/combigrid/serialization/TreeStorageSerializationStrategy.hpp>
#include <sgpp/combigrid/threading/PtrGuard.hpp>

#include <mutex>
#include <string>
#include <vector>

namespace sgpp {
namespace combigrid {

class CombigridTreeStorageImpl {
 public:
  CombigridTreeStorageImpl(
      std::vector<std::shared_ptr<AbstractPointHierarchy>> const &p_pointHierarchies,
      MultiFunction p_func)
      : func(p_func),
        pointHierarchies(p_pointHierarchies),
        storage(new TreeStorage<std::shared_ptr<TreeStorage<double>>>(
            p_pointHierarchies.size(),
            [](MultiIndex const &level) { return std::shared_ptr<TreeStorage<double>>(nullptr); })),
        mutexPtr(nullptr) {
    setFunctions();
  }

  /**
   * Sets the computation functions for the storage and the storages it contains
   */
  void setFunctions() {
    auto innerLambda = [this](MultiIndex const &index, MultiIndex const &level) -> double {
      std::shared_ptr<base::DataVector> coordinates;
      {
        PtrGuard guard(this->mutexPtr);
        size_t numDimensions = pointHierarchies.size();
        coordinates = std::make_shared<base::DataVector>(numDimensions);
        for (size_t d = 0; d < numDimensions; ++d) {
          (*coordinates)[d] = pointHierarchies[d]->getPoint(level[d], index[d]);
        }
      }

      return func(*coordinates);
    };

    auto outerLambda = [innerLambda, this](MultiIndex const &level) {
      // capture level by copy because the reference might not be valid anymore at the time the
      // lambda is called
      return std::shared_ptr<TreeStorage<double>>(new TreeStorage<double>(
          pointHierarchies.size(), [innerLambda, level, this](MultiIndex const &index) -> double {
            return innerLambda(index, level);
          }));
    };

    storage->setFunc(outerLambda);

    for (auto it = storage->getStoredDataIterator(); it->isValid(); it->moveToNext()) {
      auto level = it->getMultiIndex();
      it->value()->setFunc([innerLambda, level, this](MultiIndex const &index) -> double {
        return innerLambda(index, level);
      });
    }
  }

  std::function<double(base::DataVector const &)> func;
  std::vector<std::shared_ptr<AbstractPointHierarchy>> pointHierarchies;
  std::shared_ptr<TreeStorage<std::shared_ptr<TreeStorage<double>>>> storage;
  std::shared_ptr<std::mutex> mutexPtr;
};

CombigridTreeStorage::CombigridTreeStorage(
    std::vector<std::shared_ptr<AbstractPointHierarchy>> const &p_pointHierarchies,
    MultiFunction p_func)
    : impl(new CombigridTreeStorageImpl(p_pointHierarchies, p_func)) {}

CombigridTreeStorage::~CombigridTreeStorage() {}

std::shared_ptr<AbstractMultiStorageIterator<double>> CombigridTreeStorage::getGuidedIterator(
    const MultiIndex &level, MultiIndexIterator &iterator,
    std::vector<bool> orderingConfiguration) {
  // set level to zero for all nested hierarchies
  MultiIndex reducedLevel = level;
  size_t numDimensions = impl->pointHierarchies.size();
  IterationPolicy policy;

  for (size_t d = 0; d < numDimensions; ++d) {
    if (impl->pointHierarchies[d]->isNested()) {
      reducedLevel[d] = 0;
    }

    if (orderingConfiguration[d]) {
      policy.setIterator(d, impl->pointHierarchies[d]->getSortedPermutationIterator(level[d]));
    }
  }

  return impl->storage->get(reducedLevel)->getGuidedIterator(iterator, policy);
}

size_t CombigridTreeStorage::getNumEntries() {
  size_t result = 0;

  auto it = impl->storage->getStoredDataIterator();

  while (it->isValid()) {
    auto innerIt = it->value()->getStoredDataIterator();

    while (innerIt->isValid()) {
      ++result;
      innerIt->moveToNext();
    }

    it->moveToNext();
  }

  return result;
}

std::string CombigridTreeStorage::serialize() {
  std::shared_ptr<AbstractSerializationStrategy<double>> floatSerializationStrategy(
      new FloatSerializationStrategy<double>());

  std::shared_ptr<AbstractSerializationStrategy<std::shared_ptr<TreeStorage<double>>>>
  innerSerializationStrategy(new TreeStorageSerializationStrategy<double>(
      impl->pointHierarchies.size(), floatSerializationStrategy));

  TreeStorageSerializationStrategy<std::shared_ptr<TreeStorage<double>>> outerSerializationStrategy(
      impl->pointHierarchies.size(), innerSerializationStrategy);

  return outerSerializationStrategy.serialize(impl->storage);  // TODO(holzmudd)
}

void CombigridTreeStorage::deserialize(const std::string &str) {
  std::shared_ptr<AbstractSerializationStrategy<double>> floatSerializationStrategy(
      new FloatSerializationStrategy<double>());

  std::shared_ptr<AbstractSerializationStrategy<std::shared_ptr<TreeStorage<double>>>>
  innerSerializationStrategy(new TreeStorageSerializationStrategy<double>(
      impl->pointHierarchies.size(), floatSerializationStrategy));

  TreeStorageSerializationStrategy<std::shared_ptr<TreeStorage<double>>> outerSerializationStrategy(
      impl->pointHierarchies.size(), innerSerializationStrategy);

  impl->storage = outerSerializationStrategy.deserialize(str);

  impl->setFunctions();
}

void CombigridTreeStorage::set(const MultiIndex &level, const MultiIndex &index, double value) {
  impl->storage->get(level)->set(index, value);
}
}  // namespace combigrid
} /* namespace sgpp*/

void sgpp::combigrid::CombigridTreeStorage::setMutex(std::shared_ptr<std::mutex> mutexPtr) {
  impl->mutexPtr = mutexPtr;
}
