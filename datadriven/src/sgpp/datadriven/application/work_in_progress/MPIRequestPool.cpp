//
// Created by Vincent_Bode on 11.08.2017.
//

#include "MPIRequestPool.hpp"

namespace sgpp {
    namespace datadriven {
        size_t MPIRequestPool::createMPIRequestHandle() {
            mpiRequestStorage.emplace_back();
            return mpiRequestStorage.size() - 1;
        }

        void MPIRequestPool::deleteMPIRequestHandle(size_t handleIndex) {
            mpiRequestStorage[handleIndex] = MPI_REQUEST_NULL;

            std::cout << "MPI_Request pool size is " << mpiRequestStorage.size();

            if (handleIndex == mpiRequestStorage.size() - 1) {
                auto iterator = mpiRequestStorage.end();
                auto vectorStart = mpiRequestStorage.begin();
                iterator--;
                while (*iterator == MPI_REQUEST_NULL && iterator != vectorStart) {
                    mpiRequestStorage.erase(iterator);
                    iterator--;
                }
                std::cout << "New MPI_Request pool size is " << mpiRequestStorage.size();
            }
        }

        MPI_Request *MPIRequestPool::getMPIRequestHandle(size_t handleIndex) {
            return &mpiRequestStorage[handleIndex];
        }

        MPI_Request *MPIRequestPool::getMPIRequests() {
            return getMPIRequestHandle(0);
        }
    }
}