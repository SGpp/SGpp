//
// Created by Vincent Bode on 08.08.2017.
//

#include <sgpp/datadriven/application/work_in_progress/MPIMethods.hpp>
#include <sgpp/datadriven/application/work_in_progress/RoundRobinScheduler.hpp>

namespace sgpp {
    namespace datadriven {
        RoundRobinScheduler::RoundRobinScheduler(size_t batchSize) {
            this->batchSize = batchSize;
            this->lastWorkerID = 0;
        }

        void RoundRobinScheduler::assignTaskVariableTaskSize(TaskType taskType, AssignTaskResult &result) {
            assignTaskStaticTaskSize(taskType, result);
            result.taskSize = batchSize;
            numOutstandingBatchTracker.reserve(1);
        }

        void RoundRobinScheduler::assignTaskStaticTaskSize(TaskType taskType, AssignTaskResult &result) {
            if (lastWorkerID + 1 >= MPIMethods::getWorldSize()) {
                lastWorkerID = 0;
            }
            if (taskType == TRAIN_FROM_BATCH) {
                // Assigning Batch
                numOutstandingBatchTracker.back() += learnerInstance->getNumClasses();
            }
            result.workerID = ++lastWorkerID;
        }

        bool RoundRobinScheduler::isReadyForRefinement() {
            //Limit number of outstanding batch requests
            bool outStandingOK = numOutstandingBatchTracker.back() / learnerInstance->getNumClasses() <
                                 MPIMethods::getWorldSize() * 2;
            if (numOutstandingBatchTracker.size() <= 1) {
                return outStandingOK;
            }
            return numOutstandingBatchTracker[numOutstandingBatchTracker.size() - 2] == 0 && outStandingOK;
        }

        void RoundRobinScheduler::onMergeRequestIncoming(size_t batchOffset, size_t batchSize) {
            numOutstandingBatchTracker.back()--;
        }

        void RoundRobinScheduler::onRefinementStarted() {
            numOutstandingBatchTracker.emplace_back(0);
        }
    }
}
