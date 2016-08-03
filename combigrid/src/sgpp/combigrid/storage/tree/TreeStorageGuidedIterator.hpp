// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef COMBIGRID_SRC_SGPP_COMBIGRID_STORAGE_TREE_TREESTORAGEGUIDEDITERATOR_HPP_
#define COMBIGRID_SRC_SGPP_COMBIGRID_STORAGE_TREE_TREESTORAGEGUIDEDITERATOR_HPP_

#include <sgpp/combigrid/storage/AbstractMultiStorageIterator.hpp>
#include <sgpp/combigrid/storage/tree/InternalTreeStorageNode.hpp>
#include <sgpp/combigrid/storage/tree/LowestTreeStorageNode.hpp>
#include <sgpp/combigrid/storage/IterationPolicy.hpp>
#include <sgpp/combigrid/common/MultiIndexIterator.hpp>

#include <sgpp/combigrid/definitions.hpp>
#include <vector>
#include <memory>
#include <cstddef>

namespace sgpp {
namespace combigrid {

/**
 * Iterator class that travels "along" a MultiIndexIterator through a TreeStorage.
 * If entries are not already contained, they are created during iteration.
 */
template <typename T>
class TreeStorageGuidedIterator : public AbstractMultiStorageIterator<T> {
  MultiIndexIterator &iterator;
  MultiIndex permutedIndex;

  std::vector<InternalTreeStorageNode<T> *> internalNodes;
  LowestTreeStorageNode<T> *lowestNode;

  bool valid;

  IterationPolicy policy;

  AbstractTreeStorageNode<T> *getChild(InternalTreeStorageNode<T> *node, size_t depth, size_t index,
                                       MultiIndex const &multiIndex) {
    if (index >= node->children.size()) {
      size_t numDimensions = permutedIndex.size();
      size_t remainingDimensions = numDimensions - depth - 1;

      if (remainingDimensions >= 2) {
        while (index >= node->children.size()) {
          CGLOG("TreeStorageGuidedIterator::getChild(): create internal tree storage node at depth "
                << depth);
          node->children.emplace_back(new InternalTreeStorageNode<T>(
              node->context, remainingDimensions - 1, multiIndex, depth + 1, depth, true));
        }
      } else {
        while (index >= node->children.size()) {
          CGLOG("TreeStorageGuidedIterator::getChild(): create lowest tree storage node at depth "
                << depth);
          node->children.emplace_back(
              new LowestTreeStorageNode<T>(node->context, multiIndex, depth + 1, depth, true));
        }
      }
    }

    return node->children[index].get();
  }

 public:
  /**
   * The number of dimensions is also contained in the iterator, but it is passed separately,
   * so the casts below can't fail because of an iterator with a wrong number of dimensions.
   */
  TreeStorageGuidedIterator(IterationPolicy const &policy, AbstractTreeStorageNode<T> *root,
                            size_t numDimensions, MultiIndexIterator &iterator)
      : iterator(iterator),
        permutedIndex(numDimensions, 0),
        internalNodes(),
        lowestNode(nullptr),
        valid(true),
        policy(policy) {
    size_t lastDim = numDimensions - 1;
    auto current = root;
    for (size_t d = 0; d < lastDim; ++d) {
      auto internal = (InternalTreeStorageNode<T> *)current;
      internalNodes.push_back(internal);
      size_t currentIndex = this->policy.value(d, 0);
      CGLOG("TreeStorageGuidedIterator(): currentIndex == " << currentIndex);
      permutedIndex[d] = currentIndex;
      current = getChild(internal, d, currentIndex, permutedIndex);
    }
    lowestNode = (LowestTreeStorageNode<T> *)current;
  }

  virtual ~TreeStorageGuidedIterator() {}

  /**
   * @return Returns the difference of the greatest dimension and the lowest updated dimension,
   * i. e. the lowest dimension where the corresponding multi-index changed.
   * This will be mostly be zero, if the last dimension has many points.
   * Returns -1 if the next entry is invalid.
   */
  virtual int moveToNext() {
    size_t numInternalDimensions = internalNodes.size();
    size_t lastDim = numInternalDimensions;  // for access to the permutation iterator in the policy
    size_t lastInternalDim = numInternalDimensions - 1;
    size_t numDimensions = numInternalDimensions + 1;

    int h = iterator.moveToNext();

    if (h == 0) {
      policy.moveToNext(lastDim);
      return 0;
    } else if (h < 0) {
      return h;
    }

    policy.reset(lastDim);

    // update all pointers in dimensions dim >= numDimensions - h
    // we assume dim >= 1 since h < numDimensions (otherwise the iteration would stop and then h ==
    // -1)
    size_t dim = numDimensions - h;
    size_t d = dim - 1;
    size_t nextValue = policy.moveAndGetValue(d, iterator.indexAt(d));
    permutedIndex[d] = nextValue;

    if (dim < numInternalDimensions) {
      internalNodes[dim] =
          (InternalTreeStorageNode<T> *)getChild(internalNodes[d], d, nextValue, permutedIndex);
      nextValue = policy.resetAndGetValue(dim, 0);
      permutedIndex[dim] = nextValue;

      ++dim;
      while (dim < numInternalDimensions) {
        internalNodes[dim] = (InternalTreeStorageNode<T> *)getChild(internalNodes[dim - 1], dim - 1,
                                                                    nextValue, permutedIndex);
        nextValue = policy.resetAndGetValue(dim, 0);
        permutedIndex[dim] = nextValue;

        ++dim;
      }
    }
    lowestNode = (LowestTreeStorageNode<T> *)getChild(internalNodes[lastInternalDim],
                                                      lastInternalDim, nextValue, permutedIndex);

    return h;
  }

  virtual T &value() {
    size_t lastDim = internalNodes.size();
    size_t permutedLowestIndex = policy.value(lastDim, iterator.indexAt(lastDim));
    permutedIndex[lastDim] = permutedLowestIndex;
    return lowestNode->get(permutedIndex, lastDim);
  }

  /**
   * @return returns true if the iterator points to a valid position.
   */
  virtual bool isValid() { return iterator.isValid(); }

  virtual size_t indexAt(size_t d) const { return iterator.indexAt(d); }

  virtual MultiIndex getMultiIndex() const { return iterator.getMultiIndex(); }

  virtual bool computationRequested() {
    size_t lastDim = internalNodes.size();
    size_t permutedLowestIndex = policy.value(lastDim, iterator.indexAt(lastDim));
    return permutedLowestIndex < lowestNode->statusVector.size() &&
           lowestNode->statusVector[permutedLowestIndex] >= StorageStatus::REQUESTED;
  }

  virtual std::function<T()> requestComputationTask() {
    size_t lastDim = internalNodes.size();
    size_t permutedLowestIndex = policy.value(lastDim, iterator.indexAt(lastDim));
    lowestNode->ensureVectorEntry(permutedLowestIndex);
    permutedIndex[lastDim] = permutedLowestIndex;
    lowestNode->statusVector[permutedLowestIndex] = StorageStatus::REQUESTED;

    auto myLowestNode = lowestNode;
    auto myPermutedIndex = permutedIndex;
    return
        [myLowestNode, myPermutedIndex]() { return myLowestNode->context.func(myPermutedIndex); };
  }
};

} /* namespace combigrid */
} /* namespace sgpp*/

#endif /* COMBIGRID_SRC_SGPP_COMBIGRID_STORAGE_TREE_TREESTORAGEGUIDEDITERATOR_HPP_ */
