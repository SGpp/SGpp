// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/globaldef.hpp>

#include <sgpp/datadriven/application/learnersgdeonoffparallel/MPITaskScheduler.hpp>

namespace sgpp {
namespace datadriven {
class RoundRobinScheduler : public MPITaskScheduler {
 public:
  /**
   * Create a Round-Robin-Scheduler that will distribute tasks to all workers fairly with
   * specified batch size.
   * @param batchSize The size of one training batch to distribute.
   */
  explicit RoundRobinScheduler(size_t batchSize);

  /**
   * Assign a task of static size to the next worker in the queue.
   * @param taskType Type of task to assign to a worker.
   * @param result The result of determining assignment.
   */
  void assignTaskStaticTaskSize(TaskType taskType, AssignTaskResult &result) override;

  /**
 * Assign a task of variable size equal to the batch size to the next worker in the queue.
 * @param taskType Type of task to assign to a worker.
 * @param result The result of determining assignment.
 */
  void assignTaskVariableTaskSize(TaskType taskType, AssignTaskResult &result) override;

  /**
 * Check whether the master can start to refine. This can only happen if all requests from one
 * cycle ago have successfully completed.
 * @return Whether to start refining.
 */
  bool isReadyForRefinement() override;

  /**
   * Update the number of outstanding requests when a request is completed by a worker.
   * The difference in grid versions is used to determine whether to update the current number
   * of outstanding requests or the previous number of outstanding requests.
   * @param batchOffset Not used.
   * @param batchSize Not used.
   * @param remoteGridVersion The grid version that was used to train the batch.
   * @param localGridVersion The current grid version on the master.
   */
  void onMergeRequestIncoming(size_t batchOffset, size_t batchSize,
                              size_t remoteGridVersion, size_t localGridVersion) override;

  /**
   * Move the number of current outstanding requests into the number of previous outstanding
   * requests.
   */
  void onRefinementStarted() override;

 protected:
  /**
   * Keeps track of the last worker that had work assigned to it. The next assignment will
   * go to a worker with an index of one higher.
   */
  unsigned int lastWorkerID;
  /**
   * The batch size to use for all workers.
   */
  size_t batchSize;

  /**
   * Track how many outstanding requests (batch*class) are pending for the current refinement
   * cycle.
   */
  size_t numOutstandingRequestsCurrentRefinement;
  /**
   * Track how many outstanding requests (batch*class) are pending for the previous refinement
   * cycle.
   */
  size_t numOutstandingRequestsLastRefinement;
};
}  // namespace datadriven
}  // namespace sgpp
