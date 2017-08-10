//
// Created by Vincent Bode on 08.08.2017.
//

#include "MPITaskScheduler.hpp"

namespace sgpp {
    namespace datadriven {
        MPITaskScheduler &instance;


        MPITaskScheduler &MPITaskScheduler::getInstance() {
            return instance;
        }

        void MPITaskScheduler::setScheduler(MPITaskScheduler &newInstance) {
            instance = newInstance;
        }
    }
}
