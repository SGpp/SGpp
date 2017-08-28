//
// Created by Vincent Bode on 08.08.2017.
//

#include <sgpp/datadriven/application/work_in_progress/MPITaskScheduler.hpp>

namespace sgpp {
    namespace datadriven {

        void MPITaskScheduler::setLearnerInstance(LearnerSGDEOnOffParallel *instance) {
            learnerInstance = instance;
        }
    }
}
