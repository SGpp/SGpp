//
// Created by Vincent_Bode on 11.08.2017.
//

#include <cstring>
#include "MPIRequestPool.hpp"

namespace sgpp {
    namespace datadriven {
        size_t MPIRequestPool::createMPIRequestHandle() {

            printPoolStatistics();

            if (!freedRequests.empty()) {
                auto iterator = freedRequests.begin();
                size_t index = *iterator;
                freedRequests.erase(iterator);
                std::memset(mpiRequestStorage[index], 0, sizeof(MPI_Request));
//                std::cout << "Reused freed request " << index << std::endl;
                return index;
            }

            mpiRequestStorage.emplace_back();
            return mpiRequestStorage.size() - 1;

        }

        void MPIRequestPool::deleteMPIRequestHandle(size_t handleIndex) {
//            mpiRequestStorage[handleIndex] = MPI_REQUEST_NULL;

//            std::cout << "Received delete request for handle " << handleIndex << std::endl;
//            printPoolStatistics();
            freedRequests.insert(handleIndex);

            if (handleIndex == mpiRequestStorage.size() - 1) {
                auto iterator = mpiRequestStorage.end();
                auto vectorStart = mpiRequestStorage.begin();
                size_t index = mpiRequestStorage.size() - 1;
                iterator--;
                while (freedRequests.count(index) == 1 && iterator != vectorStart) {
                    mpiRequestStorage.erase(iterator);
                    freedRequests.erase(index);

                    //Reset invalidated iterator
                    iterator = mpiRequestStorage.end();
                    iterator--;
                    index--;
//                    std::cout << "Removing unused pool handle" << std::endl;
                }
                printPoolStatistics();
            }
        }

        inline void MPIRequestPool::printPoolStatistics() const {
            std::cout << "MPI_Request pool size is " << mpiRequestStorage.size() << " (vector capacity "
                      << mpiRequestStorage.capacity() << ", secondary " << freedRequests.size() << ")" << std::endl;
        }

        MPI_Request *MPIRequestPool::getMPIRequestHandle(size_t handleIndex) {
            return &mpiRequestStorage[handleIndex];
        }

        MPI_Request *MPIRequestPool::getMPIRequests() {
            return getMPIRequestHandle(0);
        }

        size_t MPIRequestPool::size() {
            return mpiRequestStorage.size();
        }
    }
}