//
// Created by Vincent Bode on 08.08.2017.
//

#ifndef SGPP_MPITASKSCHEDULER_HPP
#define SGPP_MPITASKSCHEDULER_HPP

#include <cstdlib>

namespace sgpp {
    namespace datadriven {

        enum TaskType {
            TRAIN_FROM_BATCH,
            UPDATE_CHOLESKY_DECOMPOSITION
        };
        struct AssignTaskResult {
            unsigned long workerID;
            unsigned long taskSize;
        };

        class MPITaskScheduler {


        public:
            virtual void
            assignTaskVariableTaskSize(TaskType taskType, AssignTaskResult &result);

            virtual void assignTaskStaticTaskSize(TaskType taskType, AssignTaskResult &result);

        protected:

        };

    }
}

#endif //SGPP_MPITASKSCHEDULER_HPP
