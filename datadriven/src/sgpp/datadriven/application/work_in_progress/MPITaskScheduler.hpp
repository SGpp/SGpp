// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef SGPP_MPITASKSCHEDULER_HPP
#define SGPP_MPITASKSCHEDULER_HPP

#include <cstdlib>

namespace sgpp {
namespace datadriven {
// Forward declere Learner, as we use only pointer
class LearnerSGDEOnOffParallel;

enum TaskType {
  TRAIN_FROM_BATCH,
  RECOMPUTE_CHOLESKY_DECOMPOSITION
};
struct AssignTaskResult {
  int workerID;
  size_t taskSize;
};

class MPITaskScheduler {
 public:
  virtual void
  assignTaskVariableTaskSize(TaskType taskType, AssignTaskResult &result) = 0;

  virtual void assignTaskStaticTaskSize(TaskType taskType, AssignTaskResult &result) = 0;

  virtual bool isReadyForRefinement() = 0;

  virtual void onRefinementStarted() = 0;

  virtual void
  onMergeRequestIncoming(size_t batchOffset,
                         size_t batchSize,
                         size_t remoteGridVersion,
                         size_t localGridVersion) = 0;

  void setLearnerInstance(LearnerSGDEOnOffParallel *instance);
 protected:
  LearnerSGDEOnOffParallel *learnerInstance{};
};
}  // namespace datadriven
}  // namespace sgpp

#endif  // SGPP_MPITASKSCHEDULER_HPP
