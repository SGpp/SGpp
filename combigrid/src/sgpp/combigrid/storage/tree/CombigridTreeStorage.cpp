// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/combigrid/storage/tree/CombigridTreeStorage.hpp>

#include <sgpp/combigrid/serialization/FloatSerializationStrategy.hpp>
#include <sgpp/combigrid/serialization/TreeStorageSerializationStrategy.hpp>
#include <sgpp/combigrid/threading/PtrGuard.hpp>

#include <iostream>
#include <mutex>
#include <string>
#include <vector>

namespace sgpp {
namespace combigrid {

class CombigridTreeStorageImpl {
 public:
  CombigridTreeStorageImpl(
      std::vector<std::shared_ptr<AbstractPointHierarchy>> const &p_pointHierarchies,
      MultiFunction p_func, bool exploitNesting)
      : func(p_func),
        pointHierarchies(p_pointHierarchies),
        mutexPtr(nullptr),
        exploitNesting(exploitNesting) {
    storage = std::make_shared<TreeStorage<std::shared_ptr<TreeStorage<double>>>>(
        p_pointHierarchies.size(),
        [](MultiIndex const &level) { return std::shared_ptr<TreeStorage<double>>(nullptr); });
    setFunctions();
  }

  /**
   * Sets the computation functions for the storage and the storages it contains
   */
  void setFunctions() {
    auto innerLambda = [this](MultiIndex const &index, MultiIndex const &level) -> double {
      std::shared_ptr<base::DataVector> coordinates;
      {
        CGLOG_SURROUND(PtrGuard guard(this->mutexPtr));
        size_t numDimensions = pointHierarchies.size();
        coordinates = std::make_shared<base::DataVector>(numDimensions);
        for (size_t d = 0; d < numDimensions; ++d) {
          (*coordinates)[d] = pointHierarchies[d]->getPoint(level[d], index[d]);
        }
        CGLOG("leave guard(this->mutexPtr) in CGStorage");
      }

      return func(*coordinates);
    };

    auto outerLambda = [innerLambda, this](MultiIndex const &level) {
      // capture level by copy because the reference might not be valid anymore at the time the
      // lambda is called
      return std::make_shared<TreeStorage<double>>(
          pointHierarchies.size(), [innerLambda, level, this](MultiIndex const &index) -> double {
            return innerLambda(index, level);
          });
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
  std::shared_ptr<std::recursive_mutex> mutexPtr;
  bool exploitNesting;
};

CombigridTreeStorage::CombigridTreeStorage(
    std::vector<std::shared_ptr<AbstractPointHierarchy>> const &p_pointHierarchies,
    MultiFunction p_func) {
  impl = std::make_unique<CombigridTreeStorageImpl>(p_pointHierarchies, p_func, true);
}

CombigridTreeStorage::CombigridTreeStorage(
    std::vector<std::shared_ptr<AbstractPointHierarchy>> const &p_pointHierarchies,
    bool exploitNesting, MultiFunction p_func) {
  impl = std::make_unique<CombigridTreeStorageImpl>(p_pointHierarchies, p_func, exploitNesting);
}

CombigridTreeStorage::~CombigridTreeStorage() {}

std::shared_ptr<AbstractMultiStorageIterator<double>> CombigridTreeStorage::getGuidedIterator(
    const MultiIndex &level, MultiIndexIterator &iterator,
    std::vector<bool> orderingConfiguration) {
  // set level to zero for all nested hierarchies
  MultiIndex reducedLevel = level;
  size_t numDimensions = impl->pointHierarchies.size();
  IterationPolicy policy;

  for (size_t d = 0; d < numDimensions; ++d) {
    if (impl->pointHierarchies[d]->isNested() && impl->exploitNesting) {
      reducedLevel[d] = 0;
    }

    if (orderingConfiguration[d]) {
      policy.setIterator(d, impl->pointHierarchies[d]->getSortedPermutationIterator(level[d]));
    }
  }

  return impl->storage->get(reducedLevel)->getGuidedIterator(iterator, policy);
}

std::shared_ptr<TreeStorage<double>> CombigridTreeStorage::getStorage(const MultiIndex &level) {
  // set level to zero for all nested hierarchies
  MultiIndex reducedLevel = level;
  size_t numDimensions = impl->pointHierarchies.size();
  IterationPolicy policy;

  for (size_t d = 0; d < numDimensions; ++d) {
    if (impl->pointHierarchies[d]->isNested() && impl->exploitNesting) {
      reducedLevel[d] = 0;
    }
  }

  return impl->storage->get(reducedLevel);
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
  std::shared_ptr<AbstractSerializationStrategy<double>> floatSerializationStrategy =
      std::make_shared<FloatSerializationStrategy<double>>();

  std::shared_ptr<AbstractSerializationStrategy<std::shared_ptr<TreeStorage<double>>>>
      innerSerializationStrategy = std::make_shared<TreeStorageSerializationStrategy<double>>(
          impl->pointHierarchies.size(), floatSerializationStrategy);

  TreeStorageSerializationStrategy<std::shared_ptr<TreeStorage<double>>> outerSerializationStrategy(
      impl->pointHierarchies.size(), innerSerializationStrategy);

  return outerSerializationStrategy.serialize(impl->storage);
}

void CombigridTreeStorage::deserialize(const std::string &str) {
  std::shared_ptr<AbstractSerializationStrategy<double>> floatSerializationStrategy =
      std::make_shared<FloatSerializationStrategy<double>>();

  std::shared_ptr<AbstractSerializationStrategy<std::shared_ptr<TreeStorage<double>>>>
      innerSerializationStrategy = std::make_shared<TreeStorageSerializationStrategy<double>>(
          impl->pointHierarchies.size(), floatSerializationStrategy);

  TreeStorageSerializationStrategy<std::shared_ptr<TreeStorage<double>>> outerSerializationStrategy(
      impl->pointHierarchies.size(), innerSerializationStrategy);

  impl->storage = outerSerializationStrategy.deserialize(str);

  impl->setFunctions();
}

void CombigridTreeStorage::set(const MultiIndex &level, const MultiIndex &index, double value) {
  MultiIndex reducedLevel = level;
  size_t numDimensions = impl->pointHierarchies.size();
  for (size_t d = 0; d < numDimensions; ++d) {
    if (impl->pointHierarchies[d]->isNested() && impl->exploitNesting) {
      reducedLevel[d] = 0;
    }
  }
  impl->storage->get(reducedLevel)->set(index, value);
}

double CombigridTreeStorage::get(MultiIndex const &level, MultiIndex const &index) {
  MultiIndex reducedLevel = level;
  size_t numDimensions = impl->pointHierarchies.size();
  for (size_t d = 0; d < numDimensions; ++d) {
    if (impl->pointHierarchies[d]->isNested() && impl->exploitNesting) {
      reducedLevel[d] = 0;
    }
  }
  return impl->storage->get(reducedLevel)->get(index);
}

void CombigridTreeStorage::setMutex(std::shared_ptr<std::recursive_mutex> mutexPtr) {
  impl->mutexPtr = mutexPtr;
}

}  // namespace combigrid
} /* namespace sgpp*/
