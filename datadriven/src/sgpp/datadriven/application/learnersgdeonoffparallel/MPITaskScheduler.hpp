// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/globaldef.hpp>

#include <cstdlib>

namespace sgpp {
namespace datadriven {

// Forward declare learner, as we use only pointer
class LearnerSGDEOnOffParallel;

/**
 * Used for the MPI Task Scheduler to differentiate between assigning tasks of different types.
 */
enum TaskType {
  TRAIN_FROM_BATCH,
  RECOMPUTE_SYSTEM_MATRIX_DECOMPOSITION
};
/**
 * Used by the MPI Task Scheduler to deliver the result of assigning the requested task.
 */
struct AssignTaskResult {
  /**
   * Holds the MPI rank of the worker to which the task was assigned.
   */
  int workerID;
  /**
   * Holds the requested size of the task (e.g. batch size), if applicable.
   */
  size_t taskSize;
};

class MPITaskScheduler {
 public:
  /**
   * Default destructor to override if necessary
   */
  virtual ~MPITaskScheduler() = default;

  /**
   * Request for a task of specific type to be assigned. The MPI Task Scheduler will decide on the
   * task size.
   * @param taskType The type of the task to be assigned.
   * @param result The result of assigning the task.
   */
  virtual void
  assignTaskVariableTaskSize(TaskType taskType, AssignTaskResult &result) = 0;

  /**
   * Request for a task of specific type to be assigned. This type of task does not include variable
   * task size.
   * @param taskType
   * @param result
   */
  virtual void assignTaskStaticTaskSize(TaskType taskType, AssignTaskResult &result) = 0;

  /**
   * Check whether the master is allowed to transit into the refinement phase.
   * @return Whether to start refinement.
   */
  virtual bool isReadyForRefinement() = 0;

  /**
   * Callback for when refinement is about to begin.
   */
  virtual void onRefinementStarted() = 0;

  /**
   * Callback for when training results are received from workers. Can be used to update tracking
   * of assigned batches.
   * @param batchOffset The offset of the batch that was trained from.
   * @param batchSize The size of the batch used for training.
   * @param remoteGridVersion The grid version of the worker used for training.
   * @param localGridVersion The grid version of the master upond reception.
   */
  virtual void
  onMergeRequestIncoming(size_t batchOffset,
                         size_t batchSize,
                         size_t remoteGridVersion,
                         size_t localGridVersion) = 0;

  /**
   * Set the learner instance for which to task schedule.
   * @param instance The learner instance.
   */
  void setLearnerInstance(LearnerSGDEOnOffParallel *instance);

 protected:
  /**
   * The learner instance used to interrogate system state if necessary.
   */
  LearnerSGDEOnOffParallel *learnerInstance{};
};
}  // namespace datadriven
}  // namespace sgpp
