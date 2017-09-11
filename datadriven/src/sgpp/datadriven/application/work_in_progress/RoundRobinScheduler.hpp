// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef SGPP_ROUNDROBINSCHEDULER_H
#define SGPP_ROUNDROBINSCHEDULER_H

#include <sgpp/datadriven/application/work_in_progress/MPITaskScheduler.hpp>

namespace sgpp {
namespace datadriven {
class RoundRobinScheduler : public MPITaskScheduler {
 public:
  explicit RoundRobinScheduler(size_t batchSize);

  void assignTaskStaticTaskSize(TaskType taskType, AssignTaskResult &result) override;

  void assignTaskVariableTaskSize(TaskType taskType, AssignTaskResult &result) override;

 protected:
  unsigned int lastWorkerID;
  size_t batchSize;

  size_t numOutstandingRequestsCurrentRefinement;
  size_t numOutstandingRequestsLastRefinement;
  bool isReadyForRefinement() override;

  void onMergeRequestIncoming(size_t batchOffset, size_t batchSize,
                              size_t remoteGridVersion, size_t localGridVersion) override;

  void onRefinementStarted() override;
};
}  // namespace datadriven
}  // namespace sgpp

#endif  // SGPP_ROUNDROBINSCHEDULER_H
