// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef COMBIGRID_SRC_SGPP_COMBIGRID_STORAGE_TREE_LOWESTTREESTORAGENODE_HPP_
#define COMBIGRID_SRC_SGPP_COMBIGRID_STORAGE_TREE_LOWESTTREESTORAGENODE_HPP_

#include <sgpp/combigrid/definitions.hpp>
#include <sgpp/combigrid/storage/tree/AbstractTreeStorageNode.hpp>
#include <sgpp/combigrid/storage/tree/TreeStorageContext.hpp>

#include <algorithm>
#include <vector>

namespace sgpp {
namespace combigrid {

/**
 * The Tree does not contain single elements of T in its leafs. For efficiency reasons, the lowest
 * Node contains a vector of T's instead.
 * See also InternalTreeStorageNode.
 */
template <typename T>
class LowestTreeStorageNode : public AbstractTreeStorageNode<T> {
 public:
  std::vector<StorageStatus> statusVector;
  std::vector<T> elements;
  TreeStorageContext<T> &context;

  LowestTreeStorageNode(TreeStorageContext<T> &context)
      : statusVector(), elements(), context(context) {}

  virtual ~LowestTreeStorageNode() {}

  virtual size_t numChildren() const {
    return std::count_if(statusVector.begin(), statusVector.end(),
                         [](StorageStatus s) { return s == StorageStatus::STORED; });
  }

  virtual bool isLeaf() const { return true; }

  void ensureVectorEntry(size_t currentIndex) {
    while (currentIndex >= elements.size()) {
      elements.push_back(T());
      statusVector.push_back(StorageStatus::NOT_STORED);
    }
  }

  virtual T &get(MultiIndex const &index, size_t depth = 0) {
    size_t currentIndex = index[depth];
    ensureVectorEntry(currentIndex);
    if (statusVector[currentIndex] != StorageStatus::STORED) {
      elements[currentIndex] = context.func(index);
      statusVector[currentIndex] = StorageStatus::STORED;
    }
    return elements[currentIndex];
  }

  virtual void set(MultiIndex const &index, T const &value, size_t depth = 0) {
    size_t currentIndex = index[depth];
    ensureVectorEntry(currentIndex);
    elements[currentIndex] = value;
    statusVector[currentIndex] = StorageStatus::STORED;
  }

  virtual bool containsIndex(MultiIndex const &index, size_t depth = 0) const {
    return index[depth] < statusVector.size() &&
           statusVector[index[depth]] == StorageStatus::STORED;
  }
};

} /* namespace combigrid */
} /* namespace sgpp*/

#endif /* COMBIGRID_SRC_SGPP_COMBIGRID_STORAGE_TREE_LOWESTTREESTORAGENODE_HPP_ */
