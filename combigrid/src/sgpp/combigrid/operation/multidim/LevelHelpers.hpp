// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef COMBIGRID_SRC_SGPP_COMBIGRID_OPERATION_MULTIDIM_LEVELHELPERS_HPP_
#define COMBIGRID_SRC_SGPP_COMBIGRID_OPERATION_MULTIDIM_LEVELHELPERS_HPP_

#include <sgpp/combigrid/definitions.hpp>
#include <sgpp/combigrid/storage/AbstractMultiStorage.hpp>
#include <sgpp/combigrid/storage/tree/TreeStorage.hpp>
#include <sgpp/combigrid/threading/ThreadPool.hpp>
#include <sgpp/combigrid/utils/BinaryHeap.hpp>

#include <functional>
#include <memory>
#include <queue>
#include <unordered_set>
#include <vector>

namespace sgpp {
namespace combigrid {

class QueueEntry {
 public:
  QueueEntry(MultiIndex const &level, double priority, size_t maxNewPoints)
      : level(level), priority(priority), maxNewPoints(maxNewPoints) {}

  MultiIndex level;
  double priority;
  size_t maxNewPoints;
};

class QueueComparator {
 public:
  int operator()(QueueEntry first, QueueEntry second) const {
    return std::less<double>()(first.priority, second.priority);
  }
};

// typedef std::priority_queue<QueueEntry, std::vector<QueueEntry>, QueueComparator>
// MultiIndexQueue;
/*typedef boost::heap::binomial_heap<QueueEntry, boost::heap::compare<QueueComparator>>
    MultiIndexQueue;*/
typedef sgpp::combigrid::BinaryHeap<QueueEntry, QueueComparator> MultiIndexQueue;

/**
 * Started: the computation of function values has been started
 * Terminated: the computation of function values has terminated
 * Completed: addLevel() has been called (can be later than termination if previous levels are not
 * terminated yet).
 */
enum class ComputationStage { NOT_STARTED, STARTED, TERMINATED, COMPLETED };

class LevelInfo {
 public:
  size_t numNotStartedPredecessors;
  size_t numNotCompletedPredecessors;
  std::vector<ThreadPool::Task> terminationListeners;
  ComputationStage computationStage;
  std::shared_ptr<MultiIndexQueue::handle_type> handle;
  double norm;
  size_t maxNewPoints;
  size_t numPoints;

  // TODO(holzmudd): add priority/num points

  /**
   * Creates a new level in its earliest stage. This means that it is not ready for computation yet.
   */
  explicit LevelInfo(size_t numPredecessors, size_t maxNewPoints, size_t numPoints)
      : numNotStartedPredecessors(numPredecessors),
        numNotCompletedPredecessors(numPredecessors),
        terminationListeners(),
        computationStage(ComputationStage::NOT_STARTED),
        handle(nullptr),
        norm(0.0),
        maxNewPoints(maxNewPoints),
        numPoints(numPoints) {}

  /**
   * Creates a new level in its latest stage, where everything has already been computed.
   */
  explicit LevelInfo(double norm, size_t maxNewPoints, size_t numPoints)
      : numNotStartedPredecessors(0),
        numNotCompletedPredecessors(0),
        terminationListeners(),
        computationStage(ComputationStage::COMPLETED),
        handle(nullptr),
        norm(norm),
        maxNewPoints(maxNewPoints),
        numPoints(numPoints) {}

  void setPriority(MultiIndexQueue &queue, double priority) {
    auto entry = *(*handle);
    entry.priority = priority;
    queue.update(*handle, entry);
  }
};

} /* namespace combigrid */
} /* namespace sgpp*/

#endif /* COMBIGRID_SRC_SGPP_COMBIGRID_OPERATION_MULTIDIM_LEVELHELPERS_HPP_ */
