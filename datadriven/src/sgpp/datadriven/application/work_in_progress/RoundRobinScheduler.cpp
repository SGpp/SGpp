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
        }

        void RoundRobinScheduler::assignTaskStaticTaskSize(TaskType taskType, AssignTaskResult &result) {
            if (lastWorkerID + 1 >= MPIMethods::getWorldSize()) {
                lastWorkerID = 0;
            }
            result.workerID = ++lastWorkerID;
        }
    }
}
