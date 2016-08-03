// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef COMBIGRID_SRC_SGPP_COMBIGRID_STORAGE_ABSTRACTMULTISTORAGE_HPP_
#define COMBIGRID_SRC_SGPP_COMBIGRID_STORAGE_ABSTRACTMULTISTORAGE_HPP_

#include <sgpp/combigrid/definitions.hpp>
#include <sgpp/combigrid/storage/IterationPolicy.hpp>
#include <sgpp/combigrid/storage/AbstractMultiStorageIterator.hpp>
#include <sgpp/combigrid/common/MultiIndexIterator.hpp>

#include <memory>
#include <string>

namespace sgpp {
namespace combigrid {

/**
 * Interface for storage classes which store values of type T. These values can be accessed via a
 * multi-index.
 * The (efficiency-optimized) TreeStorage class poses certain restrictions on which values are
 * automatically included into the storage if other values are included.
 */
template <typename T>
class AbstractMultiStorage {
 public:
  virtual ~AbstractMultiStorage() {}

  // virtual std::string serialize() const = 0; // TODO(holzmudd): add SerializationStrategy
  // virtual void deserialize(std::string const &) = 0;

  /**
   * @return Returns a reference to the element stored at the given MultiIndex. If the element was
   * not yet in the storage, it is created.
   */
  virtual T &get(MultiIndex const &index) = 0;

  /**
   * Unlike get(), this function does not activate a computation if there is no entry for the given
   * multi-index in the storage.
   * Instead, it directly sets the given value.
   */
  virtual void set(MultiIndex const &index, T const &value) = 0;

  /**
   * @return Returns true iff the storage contains a value for the given MultiIndex.
   */
  virtual bool containsIndex(MultiIndex const &index) const = 0;

  /**
   * @return Returns an iterator that iterates over all values that are already stored in the
   * storage.
   */
  virtual std::shared_ptr<AbstractMultiStorageIterator<T> > getStoredDataIterator() = 0;

  /**
   * @param indexIter Iterator that defines which values should be iterated over.
   * @param policy Defines for each dimension the order that the values should be traversed in.
   * @return Returns an iterator that iterates over all values that indexIter iterates over,
   * possibly creating them on-the-fly.
   * Depending on the given policy, the values are not traversed in their natural order but in an
   * order given by the iterators in the policy.
   */
  virtual std::shared_ptr<AbstractMultiStorageIterator<T> > getGuidedIterator(
      MultiIndexIterator &indexIter, IterationPolicy const &policy = IterationPolicy::Default) = 0;
};

} /* namespace combigrid */
} /* namespace sgpp*/

#endif /* COMBIGRID_SRC_SGPP_COMBIGRID_STORAGE_ABSTRACTMULTISTORAGE_HPP_ */
