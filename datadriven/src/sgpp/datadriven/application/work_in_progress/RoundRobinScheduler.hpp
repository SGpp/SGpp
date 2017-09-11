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
            int lastWorkerID;
            size_t batchSize;

            size_t numOutstandingRequestsCurrentRefinement;
            size_t numOutstandingRequestsLastRefinement;
            bool isReadyForRefinement() override;

            void onMergeRequestIncoming(unsigned long batchOffset, unsigned long batchSize,
                                        size_t remoteGridVersion, size_t localGridVersion) override;

            void onRefinementStarted() override;
        };

    }
}

#endif //SGPP_ROUNDROBINSCHEDULER_H
