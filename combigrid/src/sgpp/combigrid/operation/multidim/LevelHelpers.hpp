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
/**
 * This is a header containing helper classes for the implementation of LevelManager.
 */

/**
 * This class can be stored in a priority queue.
 */
class QueueEntry {
 public:
  QueueEntry(MultiIndex const &level, double priority, size_t maxNewPoints)
      : level(level), priority(priority), maxNewPoints(maxNewPoints) {}

  MultiIndex level;
  double priority;
  size_t maxNewPoints;
};

/**
 * This class is a comparator to compare objects of type QueueEntry.
 */
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

/**
 * use custom binary heap class because std::priority_queue does not provide methods to change an
 * element's priority.
 * We could use boost, but that would introduce an additional dependency
 */
typedef sgpp::combigrid::BinaryHeap<QueueEntry, QueueComparator> MultiIndexQueue;

/**
 * Started: the computation of function values has been started
 * Terminated: the computation of function values has terminated
 * Completed: addLevel() has been called (can be later than termination if previous levels are not
 * terminated yet).
 */
enum class ComputationStage { NOT_STARTED, STARTED, TERMINATED, COMPLETED };

/**
 * Stores necessary information about the computation status of a level. This is necessary for
 * adaptive level generation.
 */
class LevelInfo {
 public:
  /**
   * Number of predecessor levels which are in status ComputationStage::NOT_STARTED.
   * If this is > 0, then the computation for this level may not be started (to avoid evaluating the
   * function multiple times at the same point).
   */
  size_t numNotStartedPredecessors;

  /**
   * Number of predecessor levels which are not in status ComputationStage::COMPLETED
   */
  size_t numNotCompletedPredecessors;

  /**
   * The stage of computation of this level.
   */
  ComputationStage computationStage;

  /**
   * Handle to the queue entry if there is any. Via this handle, the QueueEntry can be retrieved and
   * its priority can be updated.
   */
  std::shared_ptr<MultiIndexQueue::handle_type> handle;

  /**
   * Norm of the difference (Delta) of this level (if already computed, zero otherwise).
   */
  double norm;

  /**
   * Number of new function evaluations (grid points) that have to be done to add this level.
   */
  size_t maxNewPoints;

  /**
   * Number of grid points in this level.
   */
  size_t numPoints;

  // TODO(holzmudd): add priority

  /**
   * Creates a new level in its earliest stage. This means that it is not ready for computation yet.
   */
  explicit LevelInfo(size_t numPredecessors, size_t maxNewPoints, size_t numPoints)
      : numNotStartedPredecessors(numPredecessors),
        numNotCompletedPredecessors(numPredecessors),
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
        computationStage(ComputationStage::COMPLETED),
        handle(nullptr),
        norm(norm),
        maxNewPoints(maxNewPoints),
        numPoints(numPoints) {}

  /**
   * Updates the priority of this level in the priority queue.
   */
  void setPriority(MultiIndexQueue &queue, double priority) {
    auto entry = *(*handle);
    entry.priority = priority;
    queue.update(*handle, entry);
  }
};

} /* namespace combigrid */
} /* namespace sgpp*/

#endif /* COMBIGRID_SRC_SGPP_COMBIGRID_OPERATION_MULTIDIM_LEVELHELPERS_HPP_ */
