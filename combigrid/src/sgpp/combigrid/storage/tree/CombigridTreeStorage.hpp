// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef COMBIGRID_SRC_SGPP_COMBIGRID_STORAGE_TREE_COMBIGRIDTREESTORAGE_HPP_
#define COMBIGRID_SRC_SGPP_COMBIGRID_STORAGE_TREE_COMBIGRIDTREESTORAGE_HPP_

#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/combigrid/MultiFunction.hpp>
#include <sgpp/combigrid/definitions.hpp>
#include <sgpp/combigrid/grid/hierarchy/AbstractPointHierarchy.hpp>
#include <sgpp/combigrid/storage/AbstractCombigridStorage.hpp>
#include <sgpp/combigrid/storage/tree/TreeStorage.hpp>

#include <memory>
#include <string>
#include <vector>

namespace sgpp {
namespace combigrid {

class CombigridTreeStorageImpl;

/**
 * Implementation of the AbstractCombigridStorage using TreeStorage<TreeStorage<float_t>>. For
 * further information, refer to AbstractCombigridStorage.
 */
class CombigridTreeStorage : public AbstractCombigridStorage {
  std::unique_ptr<CombigridTreeStorageImpl> impl;

 public:
  /**
   * @param p_pointHierarchies Point hierarchies generating the points at which the function should
   * be evaluated.
   * @param p_func Function generating the values that are stored in the storage.
   */
  CombigridTreeStorage(
      std::vector<std::shared_ptr<AbstractPointHierarchy>> const &p_pointHierarchies,
      MultiFunction p_func);
  virtual ~CombigridTreeStorage();

  virtual std::shared_ptr<AbstractMultiStorageIterator<double>> getGuidedIterator(
      MultiIndex const &level, MultiIndexIterator &iterator,
      std::vector<bool> orderingConfiguration);

  /**
   * Returns the number of entries (all level-index pairs) in the storage, which indicates the
   * number of function evaluations that have been done. Currently, this is an O(n) method.
   */
  virtual size_t getNumEntries();

  virtual std::string serialize();
  virtual void deserialize(std::string const &str);

  virtual void set(MultiIndex const &level, MultiIndex const &index, double value);
  virtual void setMutex(std::shared_ptr<std::mutex> mutexPtr);
};
}  // namespace combigrid
} /* namespace sgpp*/

#endif /* COMBIGRID_SRC_SGPP_COMBIGRID_STORAGE_TREE_COMBIGRIDTREESTORAGE_HPP_ */
