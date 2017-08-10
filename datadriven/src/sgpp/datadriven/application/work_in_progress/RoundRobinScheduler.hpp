//
// Created by Vincent Bode on 08.08.2017.
//

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
            size_t lastWorkerID;
            size_t batchSize;

        };

    }
}

#endif //SGPP_ROUNDROBINSCHEDULER_H
