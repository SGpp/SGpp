// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef COMBIGRID_SRC_SGPP_COMBIGRID_STORAGE_TREE_INTERNALTREESTORAGENODE_HPP_
#define COMBIGRID_SRC_SGPP_COMBIGRID_STORAGE_TREE_INTERNALTREESTORAGENODE_HPP_

#include <sgpp/combigrid/storage/tree/AbstractTreeStorageNode.hpp>
#include <sgpp/combigrid/storage/tree/LowestTreeStorageNode.hpp>

#include <vector>

namespace sgpp {
namespace combigrid {

/**
 * Each inner node of the storage contains a vector of other nodes.
 * For a d-dimensional multi-index, to find a storage entry, the components of the multi-index are
 * iterated over.
 * Each component defines the index of the next child to be traversed.
 * Because traversal would be slower if each leaf only contained one number,
 * the nodes for the last dimension contain a vector of numbers (or Ts) instead of having children
 * which themselves have one number.
 * These nodes are represented by LowestTreeStorageNode objects, all other nodes are
 * InternalTreeStorageNode objects.
 * Each InternalTreeStorageNode contains at least one child, even if the final LowestTreeStorageNode
 * contains no value. This is guaranteed by the constructor.
 */
template <typename T>
class InternalTreeStorageNode : public AbstractTreeStorageNode<T> {
 public:
  std::vector<std::unique_ptr<AbstractTreeStorageNode<T>>> children;
  TreeStorageContext<T> &context;

  /**
   * @param context TreeStorageContext<T>-object containing information about the storage.
   * @param remainingDimensions Number of dimensions that come after this node (at least 1, since
   * this is an internal node).
   */
  InternalTreeStorageNode(TreeStorageContext<T> &context, size_t remainingDimensions)
      : children(), context(context) {
    if (remainingDimensions > 1) {
      children.emplace_back(new InternalTreeStorageNode<T>(context, remainingDimensions - 1));
    } else {
      children.emplace_back(new LowestTreeStorageNode<T>(context));
    }
  }

  virtual ~InternalTreeStorageNode() {}

  virtual size_t numChildren() const { return children.size(); }

  virtual bool isLeaf() const { return false; }

  virtual T &get(MultiIndex const &index, size_t depth = 0) {
    size_t currentIndex = index[depth];
    size_t remainingDimensions = context.numDimensions - depth - 1;

    while (currentIndex >= children.size()) {
      if (remainingDimensions >= 2) {
        children.emplace_back(new InternalTreeStorageNode<T>(context, remainingDimensions - 1));
      } else {
        children.emplace_back(new LowestTreeStorageNode<T>(context));
      }
    }
    return children[currentIndex]->get(index, depth + 1);
  }

  virtual void set(MultiIndex const &index, T const &value, size_t depth = 0) {
    size_t currentIndex = index[depth];
    size_t remainingDimensions = context.numDimensions - depth - 1;

    while (currentIndex >= children.size()) {
      if (remainingDimensions >= 2) {
        children.emplace_back(new InternalTreeStorageNode<T>(context, remainingDimensions - 1));
      } else {
        children.emplace_back(new LowestTreeStorageNode<T>(context));
      }
    }

    children[currentIndex]->set(index, value, depth + 1);
  }

  virtual bool containsIndex(MultiIndex const &index, size_t depth = 0) const {
    return index[depth] < children.size() &&
           children[index[depth]]->containsIndex(index, depth + 1);
  }
};
}  // namespace combigrid
} /* namespace sgpp*/

#endif /* COMBIGRID_SRC_SGPP_COMBIGRID_STORAGE_TREE_INTERNALTREESTORAGENODE_HPP_ */
