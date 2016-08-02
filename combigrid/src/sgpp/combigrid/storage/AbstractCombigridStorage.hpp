// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef COMBIGRID_SRC_SGPP_COMBIGRID_STORAGE_ABSTRACTCOMBIGRIDSTORAGE_HPP_
#define COMBIGRID_SRC_SGPP_COMBIGRID_STORAGE_ABSTRACTCOMBIGRIDSTORAGE_HPP_

#include <sgpp/combigrid/definitions.hpp>
#include <sgpp/combigrid/common/MultiIndexIterator.hpp>
#include <sgpp/combigrid/storage/AbstractMultiStorageIterator.hpp>

#include <memory>
#include <string>
#include <vector>
#include <mutex>

namespace sgpp{
namespace combigrid {

class AbstractCombigridStorage {
 public:
  virtual ~AbstractCombigridStorage();

  /**
   * @param level Level to traverse
   * @param iterator Iterator that defines which values should be iterated over
   * @param orderingConfiguration Defines for each dimension whether the grid points in that
   * dimension should be traversed in sorted order.
   * @return Returns an iterator that iterates along the multi-indices for a given level. These
   * multi-indices are given by the parameter iterator.
   * When moveToNext() is called on the returned iterator, it also calls moveToNext() on the
   * underlying MultiIndexIterator.
   * If the values are not already stored, they are created during iteration.
   */
  virtual std::shared_ptr<AbstractMultiStorageIterator<double> > getGuidedIterator(
      MultiIndex const &level, MultiIndexIterator &iterator,
      std::vector<bool> orderingConfiguration) = 0;

  /**
   * @return Returns the number of stored values (cumulated over all levels).
   */
  virtual size_t getNumEntries() = 0;

  virtual std::string serialize() = 0;
  virtual void deserialize(std::string const &str) = 0;
  virtual void set(MultiIndex const &level, MultiIndex const &index, double value) = 0;
  virtual void setMutex(std::shared_ptr<std::mutex> mutexPtr) = 0;
};

} /* namespace combigrid */
} /* namespace sgpp*/

#endif /* COMBIGRID_SRC_SGPP_COMBIGRID_STORAGE_ABSTRACTCOMBIGRIDSTORAGE_HPP_ */
