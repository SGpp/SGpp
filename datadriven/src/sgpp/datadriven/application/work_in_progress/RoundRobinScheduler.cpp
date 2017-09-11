//
// Created by Vincent Bode on 08.08.2017.
//

#include <sgpp/datadriven/application/work_in_progress/MPIMethods.hpp>
#include <sgpp/datadriven/application/work_in_progress/RoundRobinScheduler.hpp>
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

    void RoundRobinScheduler::onMergeRequestIncoming(unsigned long /*batchOffset*/,
                                                     unsigned long /*batchSize*/,
                                                     size_t remoteGridVersion, size_t localGridVersion) {
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
    }
}
