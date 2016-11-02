// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef COMBIGRID_SRC_SGPP_COMBIGRID_OPERATION_MULTIDIM_ABSTRACTFULLGRIDEVALUATOR_HPP_
#define COMBIGRID_SRC_SGPP_COMBIGRID_OPERATION_MULTIDIM_ABSTRACTFULLGRIDEVALUATOR_HPP_

#include <sgpp/combigrid/definitions.hpp>
#include <sgpp/combigrid/storage/AbstractCombigridStorage.hpp>

#include <memory>
#include <vector>

namespace sgpp {
namespace combigrid {

/**
 * Abstract class for doing an evaluation on a full grid, yielding a value of the template type V.
 * For this type, confer the comment in algebraic/FloatArrayVector.hpp
 * This class is used inside CombigridEvaluator to do evaluations on different grids.
 */
template <typename V>
class AbstractFullGridEvaluator {
 public:
  virtual ~AbstractFullGridEvaluator() {}

  virtual V eval(MultiIndex const &level) = 0;

  virtual size_t maxNewPoints(MultiIndex const &level) = 0;

  virtual size_t numPoints(MultiIndex const &level) = 0;

  virtual std::shared_ptr<AbstractCombigridStorage> getStorage() = 0;

  virtual std::vector<ThreadPool::Task> getLevelTasks(MultiIndex const &level,
                                                      ThreadPool::Task callback) = 0;

  virtual void setMutex(std::shared_ptr<std::mutex> mutexPtr) = 0;
};

} /* namespace combigrid */
} /* namespace sgpp*/

#endif /* COMBIGRID_SRC_SGPP_COMBIGRID_OPERATION_MULTIDIM_ABSTRACTFULLGRIDEVALUATOR_HPP_ */
