// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef COMBIGRID_SRC_SGPP_COMBIGRID_OPERATION_MULTIDIM_ABSTRACTLEVELEVALUATOR_HPP_
#define COMBIGRID_SRC_SGPP_COMBIGRID_OPERATION_MULTIDIM_ABSTRACTLEVELEVALUATOR_HPP_

#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/combigrid/definitions.hpp>
#include <sgpp/combigrid/storage/tree/TreeStorage.hpp>
#include <sgpp/combigrid/threading/ThreadPool.hpp>

#include <memory>
#include <mutex>
#include <vector>

namespace sgpp {
namespace combigrid {

/**
 * This class abstracts a lot of methods of CombigridEvaluator for classes that do not want to carry
 * around its template parameter, e. g. LevelManager. AbstractLevelEvaluator provides the central
 * method addLevel() that does an evaluation on a given level to make its numerical approximation
 * more accurate.
 *
 * For function descriptions, refer to CombigridEvaluator();
 */
class AbstractLevelEvaluator {
 public:
  virtual ~AbstractLevelEvaluator();

  virtual bool addLevel(MultiIndex const &level) = 0;
  virtual std::vector<ThreadPool::Task> getLevelTasks(MultiIndex const &level,
                                                      ThreadPool::Task callback) = 0;
  virtual void setMutex(std::shared_ptr<std::mutex> mutexPtr) = 0;
  virtual bool containsLevel(MultiIndex const &level) = 0;
  virtual size_t maxNewPoints(MultiIndex const &level) = 0;
  virtual size_t numPoints(MultiIndex const &level) = 0;
  virtual double getDifferenceNorm(MultiIndex const &level) = 0;
  virtual size_t dim() const = 0;
  virtual void clear() = 0;
  virtual std::shared_ptr<TreeStorage<uint8_t>> getLevelStructure() = 0;
  virtual std::vector<base::DataVector> getAllGridPoints() = 0;
  virtual base::DataMatrix getGridPointMatrix() = 0;
};
}  // namespace combigrid
}  // namespace sgpp

#endif /* COMBIGRID_SRC_SGPP_COMBIGRID_OPERATION_MULTIDIM_ABSTRACTLEVELEVALUATOR_HPP_ */
