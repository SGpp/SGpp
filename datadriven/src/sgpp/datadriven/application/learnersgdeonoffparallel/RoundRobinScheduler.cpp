// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/datadriven/application/learnersgdeonoffparallel/MPIMethods.hpp>
#include <sgpp/datadriven/application/learnersgdeonoffparallel/RoundRobinScheduler.hpp>
#include <sgpp/datadriven/application/learnersgdeonoffparallel/LearnerSGDEOnOffParallel.hpp>
#include <sgpp/base/exception/algorithm_exception.hpp>

namespace sgpp {
namespace datadriven {
RoundRobinScheduler::RoundRobinScheduler(size_t batchSize) {
  this->batchSize = batchSize;
  this->lastWorkerID = 0;
  numOutstandingRequestsCurrentRefinement = 0;
  numOutstandingRequestsLastRefinement = 0;
  learnerInstance = nullptr;
}

void RoundRobinScheduler::assignTaskVariableTaskSize(TaskType taskType, AssignTaskResult &result) {
  assignTaskStaticTaskSize(taskType, result);
  result.taskSize = batchSize;
}

void RoundRobinScheduler::assignTaskStaticTaskSize(TaskType taskType, AssignTaskResult &result) {
  if (lastWorkerID + 1 >= MPIMethods::getWorldSize()) {
    lastWorkerID = 0;
  }
  if (taskType == TRAIN_FROM_BATCH) {
    // Assigning Batch
    numOutstandingRequestsCurrentRefinement += learnerInstance->getNumClasses();
  }
  result.workerID = ++lastWorkerID;
}

bool RoundRobinScheduler::isReadyForRefinement() {
  return numOutstandingRequestsLastRefinement == 0;
}

void RoundRobinScheduler::onMergeRequestIncoming(size_t /*batchOffset*/,
                                                 size_t /*batchSize*/,
                                                 size_t remoteGridVersion,
                                                 size_t localGridVersion) {
  if (remoteGridVersion == localGridVersion) {
    numOutstandingRequestsCurrentRefinement--;
  } else if (remoteGridVersion + 1 == localGridVersion) {
    numOutstandingRequestsLastRefinement--;
  } else {
    throw sgpp::base::algorithm_exception("Received merge request that was too old.");
  }
}

void RoundRobinScheduler::onRefinementStarted() {
  if (numOutstandingRequestsLastRefinement != 0) {
    throw sgpp::base::algorithm_exception("Refinement started illegally.");
  }
  numOutstandingRequestsLastRefinement = numOutstandingRequestsCurrentRefinement;
  numOutstandingRequestsCurrentRefinement = 0;
}
}  // namespace datadriven
}  // namespace sgpp

