//
// Created by Vincent Bode on 08.08.2017.
//

#include <sgpp/datadriven/application/work_in_progress/MPIMethods.hpp>
#include <sgpp/datadriven/application/work_in_progress/RoundRobinScheduler.hpp>

namespace sgpp {
    namespace datadriven {
        RoundRobinScheduler::RoundRobinScheduler(size_t worldSize, size_t batchSize) {
            this->worldSize = worldSize;
            this->batchSize = batchSize;
            this->lastWorkerID = 0;
        }

        void RoundRobinScheduler::assignTaskVariableTaskSize(TaskType taskType, size_t maximumTaskSize,
                                                             AssignTaskResult &result) {
            assignTaskStaticTaskSize(taskType, result);
            result.taskSize = std::min(maximumTaskSize, batchSize);
        }

        void RoundRobinScheduler::assignTaskStaticTaskSize(TaskType taskType, AssignTaskResult &result) {
            if (lastWorkerID + 1 >= worldSize) {
                lastWorkerID = 0;
            }
            result.workerID = ++lastWorkerID;
        }
    }
}
