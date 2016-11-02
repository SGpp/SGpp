// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef COMBIGRID_SRC_SGPP_COMBIGRID_STORAGE_ABSTRACTMULTISTORAGEITERATOR_HPP_
#define COMBIGRID_SRC_SGPP_COMBIGRID_STORAGE_ABSTRACTMULTISTORAGEITERATOR_HPP_

#include <sgpp/combigrid/definitions.hpp>
#include <sgpp/combigrid/threading/ThreadPool.hpp>
#include <sgpp/globaldef.hpp>

#include <cstddef>

namespace sgpp {
namespace combigrid {

/**
 * Abstract base class for iterators iterating over (existing or maybe also non-existing) entries of
 * an AbstractMultiStorage.
 */
template <typename T>
class AbstractMultiStorageIterator {
 public:
  virtual ~AbstractMultiStorageIterator() {}

  /**
   * @return Returns the difference of the greatest dimension and the lowest updated dimension,
   * i. e. the lowest dimension where the corresponding multi-index changed.
   * This will be mostly be zero, if the last dimension has many points.
   * Returns -1 if the next entry is invalid.
   */
  virtual int moveToNext() = 0;

  /**
   * @return Returns a reference to the value the iterator points to
   */
  virtual T &value() = 0;

  /**
   * @return Returns true iff the iterator points to a valid position.
   */
  virtual bool isValid() = 0;

  /**
   * @return Returns the index in the given dimension that the iterator points to.
   * If this iterator was configured with a certain policy to permute the indices it accesses, this
   * does not affect the index returned.
   * For example, in the beginning, indexAt(0) will always return 0, even if the index is 5 after
   * the permutation.
   */
  virtual size_t indexAt(size_t d) const = 0;

  /**
   * @return Returns the MultiIndex containing the components that can be retrieved using indexAt().
   */
  virtual MultiIndex getMultiIndex() const = 0;

  /**
   * @return True iff the value the iterator points to has already been computed or the computation
   * has been requested (in case of parallel evaluation, using requestComputationTask()).
   */
  virtual bool computationRequested() = 0;

  /**
   * @return Returns a function that computes the value at the current entry (if it is not already
   * stored). This is useful for parallel function evaluation.
   */
  virtual std::function<T()> requestComputationTask() = 0;
};

} /* namespace combigrid */
} /* namespace sgpp*/

#endif /* COMBIGRID_SRC_SGPP_COMBIGRID_STORAGE_ABSTRACTMULTISTORAGEITERATOR_HPP_ */
