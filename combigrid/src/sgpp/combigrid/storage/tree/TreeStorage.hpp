// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef COMBIGRID_SRC_SGPP_COMBIGRID_STORAGE_TREE_TREESTORAGE_HPP_
#define COMBIGRID_SRC_SGPP_COMBIGRID_STORAGE_TREE_TREESTORAGE_HPP_

#include <sgpp/combigrid/storage/AbstractMultiStorage.hpp>
#include <sgpp/combigrid/storage/tree/AbstractTreeStorageNode.hpp>
#include <sgpp/combigrid/storage/tree/LowestTreeStorageNode.hpp>
#include <sgpp/combigrid/storage/tree/InternalTreeStorageNode.hpp>
#include <sgpp/combigrid/storage/tree/TreeStorageContext.hpp>
#include <sgpp/combigrid/storage/tree/TreeStorageStoredDataIterator.hpp>
#include <sgpp/combigrid/storage/tree/TreeStorageGuidedIterator.hpp>

#include <memory>
#include <stdexcept>
#include <thread>

namespace sgpp {
namespace combigrid {

/**
 * The TreeStorage class is an optimized implementation of a storage that stores values addressed by
 * MultiIndex-objects.
 * It uses a tree with as many levels as dimensions in the multi-indices.
 * TreeStorage can be configured with a function that computes a value for a multi-index such that
 * TreeStorage caches its values.
 * Multi-Indices start from 0 in each dimension.
 */
template <typename T>
class TreeStorage : public AbstractMultiStorage<T> {
  TreeStorageContext<T>
      context;  // the context stores general information and is referenced by all nodes

  std::unique_ptr<AbstractTreeStorageNode<T>> root;

  TreeStorage(TreeStorage<T> const &) = delete;

 public:
  typedef std::function<T(MultiIndex const &)> function_type;

  /**
   * Constructor.
   * @param numDimensions number of dimensions of the multi-indices that the storage is addressed
   * with
   * @param func "Default-value-function" that is called to compute entries that are not already
   * stored.
   * The storage might serve as a function value cache such that these values only have to be
   * computed once.
   * If no function is specified, the default constructor of the stored type is used.
   */
  TreeStorage(size_t numDimensions, function_type func = multiIndexToDefaultValue<T>())
      : context(numDimensions, func), root(nullptr) {
    MultiIndex index;
    if (numDimensions <= 1) {
      root.reset(new LowestTreeStorageNode<T>(context, index, 0, 0, false));
    } else {
      root.reset(new InternalTreeStorageNode<T>(context, numDimensions - 1, index, 0, 0, false));
    }
  }

  virtual ~TreeStorage() {}

  virtual T &get(MultiIndex const &index) {
    if (index.size() != context.numDimensions) {
      throw std::runtime_error("TreeStorage::get(): index.size() != context.numDimensions");
    }
    return root->get(index);
  }

  virtual void setFunc(function_type func) { context.func = func; }

  /**
   * Unlike get(), this function does not activate a computation if there is no entry for the given
   * multi-index in the storage.
   * Instead, it directly sets the given value.
   */
  virtual void set(MultiIndex const &index, T const &value) { root->set(index, value); }

  virtual bool containsIndex(MultiIndex const &index) const { return root->containsIndex(index); }

  virtual std::shared_ptr<AbstractMultiStorageIterator<T>> getStoredDataIterator() {
    return std::shared_ptr<AbstractMultiStorageIterator<T>>(
        new TreeStorageStoredDataIterator<T>(root.get(), context.numDimensions));
  }

  virtual std::shared_ptr<AbstractMultiStorageIterator<T>> getGuidedIterator(
      MultiIndexIterator &indexIter, IterationPolicy const &policy = IterationPolicy::Default) {
    return std::shared_ptr<AbstractMultiStorageIterator<T>>(
        new TreeStorageGuidedIterator<T>(policy, root.get(), context.numDimensions, indexIter));
  }
};

} /* namespace combigrid */
} /* namespace sgpp*/

#endif /* COMBIGRID_SRC_SGPP_COMBIGRID_STORAGE_TREE_TREESTORAGE_HPP_ */
