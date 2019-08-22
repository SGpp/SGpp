// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef COMBIGRID_SRC_SGPP_COMBIGRID_STORAGE_TREE_TREESTORAGESTOREDDATAITERATOR_HPP_
#define COMBIGRID_SRC_SGPP_COMBIGRID_STORAGE_TREE_TREESTORAGESTOREDDATAITERATOR_HPP_

#include <sgpp/combigrid/storage/AbstractMultiStorageIterator.hpp>
#include <sgpp/combigrid/storage/IterationPolicy.hpp>
#include <sgpp/combigrid/storage/tree/InternalTreeStorageNode.hpp>
#include <sgpp/combigrid/storage/tree/LowestTreeStorageNode.hpp>

#include <cstddef>
#include <memory>
#include <vector>

namespace sgpp {
namespace combigrid {

/**
 * Iterator for the TreeStorage class that only traverses entries stored in the storage.
 * The entries are traversed in lexicographical order (lexicographically ascending multi-indices).
 * For a detailed method description, see AbstractMultiStorageIterator.
 */
template <typename T>
class TreeStorageStoredDataIterator : public AbstractMultiStorageIterator<T> {
  MultiIndex internalIndex;
  std::vector<InternalTreeStorageNode<T> *> internalNodes;
  LowestTreeStorageNode<T> *lowestNode;

  /**
   * This is the index for the lowestNode
   */
  size_t lowestIndex;

  bool valid;

  /**
   * Helper function that moves to the next storage entry. This entry might not have been computed
   * yet and so this function might have to be called multiple times until a valid entry is reached.
   */
  int moveToNextImpl() {
    size_t numInternalDimensions = internalIndex.size();
    size_t lastInternalDim = numInternalDimensions - 1;

    size_t nextLowestIndex = lowestIndex + 1;
    if (nextLowestIndex < lowestNode->elements.size()) {
      lowestIndex = nextLowestIndex;

      return 0;
    }

    // in all other cases, we are at the end in the last dimension, so we have to advance in another
    // dimension
    lowestIndex = 0;

    int d = static_cast<int>(numInternalDimensions) - 1;
    while (d >= 0 && internalIndex[d] + 1 >= internalNodes[d]->children.size()) {
      internalIndex[d] = 0;
      --d;
    }

    if (d < 0) {
      valid = false;
      return -1;
    }

    size_t firstValue = internalIndex[d] + 1;
    internalIndex[d] = firstValue;
    size_t dim = d + 1;
    if (dim < numInternalDimensions) {
      internalNodes[dim] = dynamic_cast<InternalTreeStorageNode<T>*>(
          internalNodes[d]->children[firstValue].get());

      ++dim;
      while (dim < numInternalDimensions) {
        internalNodes[dim] = dynamic_cast<InternalTreeStorageNode<T>*>(
            internalNodes[dim - 1]->children[0].get());
        ++dim;
      }

      firstValue = 0;
    }

    lowestNode = dynamic_cast<LowestTreeStorageNode<T>*>(
        internalNodes[lastInternalDim]->children[firstValue].get());

    return static_cast<int>(numInternalDimensions - d);
  }

 public:
  TreeStorageStoredDataIterator(AbstractTreeStorageNode<T> *root, size_t numDimensions)
      : internalIndex(numDimensions - 1, 0),
        internalNodes(),
        lowestNode(nullptr),
        lowestIndex(0),
        valid(true) {
    auto current = root;

    for (; numDimensions >= 2; --numDimensions) {
      auto internal = dynamic_cast<InternalTreeStorageNode<T>*>(current);

      if (internal->children.size() == 0) {
        // the storage is empty
        valid = false;
        return;
      }

      internalNodes.push_back(internal);
      current = internal->children[0].get();
    }
    lowestNode = dynamic_cast<LowestTreeStorageNode<T>*>(current);

    // valid = lowestNode->elements.size() > 0;
    if (lowestNode->elements.size() == 0) {
      moveToNext();
    }
  }

  virtual ~TreeStorageStoredDataIterator() {}

  /**
   * @return Returns the difference of the greatest dimension and the lowest updated dimension,
   * i. e. the lowest dimension where the corresponding multi-index changed.
   * This will be mostly be zero, if the last dimension has many points.
   * Returns -1 if the next entry is invalid.
   */
  virtual int moveToNext() {
    // return moveToNextImpl();

    int h = 0;
    do {
      int newH = moveToNextImpl();
      if (newH == -1) {
        return -1;
      }
      if (newH > h) {
        h = newH;
      }
    } while (lowestIndex >= lowestNode->statusVector.size() ||
             lowestNode->statusVector[lowestIndex] != StorageStatus::STORED);

    return h;
  }

  virtual T &value() { return lowestNode->elements[lowestIndex]; }

  virtual void setValue(T const &input) { lowestNode->elements[lowestIndex] = input; }

  /**
   * @return returns true if the iterator points to a valid position.
   */
  virtual bool isValid() { return valid; }

  virtual size_t indexAt(size_t d) const {
    if (d < internalIndex.size()) {
      return internalIndex[d];
    }
    return lowestIndex;
  }

  virtual MultiIndex getMultiIndex() const {
    MultiIndex index = internalIndex;
    index.push_back(lowestIndex);
    return index;
  }

  /**
   * Returns true because the iterator only iterates over already stored data.
   */
  virtual bool computationRequested() { return true; }

  /**
   * Returns a dummy function because the values pointed to are already stored.
   */
  virtual std::function<T()> requestComputationTask() {
    return []() { return T(); };  // no computation necessary
  }
};

} /* namespace combigrid */
} /* namespace sgpp*/

#endif /* COMBIGRID_SRC_SGPP_COMBIGRID_STORAGE_TREE_TREESTORAGESTOREDDATAITERATOR_HPP_ */
