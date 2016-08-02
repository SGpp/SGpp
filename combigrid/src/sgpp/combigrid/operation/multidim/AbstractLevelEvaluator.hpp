/*
 * AbstractLevelEvaluator.hpp
 *
 *  Created on: 02.08.2016
 *      Author: david
 */

#ifndef COMBIGRID_SRC_SGPP_COMBIGRID_OPERATION_MULTIDIM_ABSTRACTLEVELEVALUATOR_HPP_
#define COMBIGRID_SRC_SGPP_COMBIGRID_OPERATION_MULTIDIM_ABSTRACTLEVELEVALUATOR_HPP_

#include <memory>
#include <mutex>
#include <sgpp/combigrid/definitions.hpp>
#include <sgpp/combigrid/storage/tree/TreeStorage.hpp>
#include <sgpp/combigrid/threading/ThreadPool.hpp>
#include <vector>

namespace sgpp {
namespace combigrid {

class AbstractLevelEvaluator {
 public:
  virtual ~AbstractLevelEvaluator();

  virtual bool addLevel(MultiIndex const &level) = 0;
  virtual std::vector<ThreadPool::Task> getLevelTasks(MultiIndex const &level,
                                                      ThreadPool::Task callback) = 0;
  virtual void setMutex(std::shared_ptr<std::mutex> mutexPtr) = 0;
  virtual bool containsLevel(MultiIndex const &level) = 0;
  virtual size_t maxNewPoints(MultiIndex const &level) = 0;
  virtual double getDifferenceNorm(MultiIndex const &level) = 0;
  virtual size_t dim() const = 0;
  virtual void clear() = 0;
  virtual std::shared_ptr<TreeStorage<uint8_t>> getLevelStructure() = 0;
};
}
}

#endif /* COMBIGRID_SRC_SGPP_COMBIGRID_OPERATION_MULTIDIM_ABSTRACTLEVELEVALUATOR_HPP_ */
