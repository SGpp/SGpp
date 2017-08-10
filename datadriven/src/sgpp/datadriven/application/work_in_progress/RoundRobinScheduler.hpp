//
// Created by Vincent Bode on 08.08.2017.
//

#ifndef SGPP_ROUNDROBINSCHEDULER_H
#define SGPP_ROUNDROBINSCHEDULER_H

#include <sgpp/datadriven/application/work_in_progress/MPITaskScheduler.hpp>

class RoundRobinScheduler : public MPITaskScheduler {

public:
    RoundRobinScheduler(size_t worldSize, size_t batchSize);

    void assignTaskVariableTaskSize(TaskType taskType, size_t maximumTaskSize, AssignTaskResult &result);

    void assignTaskStaticTaskSize(TaskType taskType, AssignTaskResult &result);

protected:
    size_t lastWorkerID;
    size_t batchSize;
    size_t worldSize;
};

#endif //SGPP_ROUNDROBINSCHEDULER_H
