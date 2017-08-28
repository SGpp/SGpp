//
// Created by Vincent Bode on 08.08.2017.
//

#ifndef SGPP_MPITASKSCHEDULER_HPP
#define SGPP_MPITASKSCHEDULER_HPP

#include <cstdlib>

namespace sgpp {
    namespace datadriven {

        //Forward declere Learner, as we use only pointer
        class LearnerSGDEOnOffParallel;

        enum TaskType {
            TRAIN_FROM_BATCH,
            RECOMPUTE_CHOLESKY_DECOMPOSITION
        };
        struct AssignTaskResult {
            unsigned long workerID;
            unsigned long taskSize;
        };

        class MPITaskScheduler {


        public:
            virtual void
            assignTaskVariableTaskSize(TaskType taskType, AssignTaskResult &result) = 0;

            virtual void assignTaskStaticTaskSize(TaskType taskType, AssignTaskResult &result) = 0;

            virtual bool isReadyForRefinement() = 0;

            virtual void onRefinementStarted() = 0;

            virtual void
            onMergeRequestIncoming(unsigned long batchOffset, unsigned long batchSize, size_t remoteGridVersion,
                                   size_t localGridVersion) = 0;

            void setLearnerInstance(LearnerSGDEOnOffParallel *instance);
        protected:
            LearnerSGDEOnOffParallel *learnerInstance;
        };

    }
}

#endif //SGPP_MPITASKSCHEDULER_HPP
