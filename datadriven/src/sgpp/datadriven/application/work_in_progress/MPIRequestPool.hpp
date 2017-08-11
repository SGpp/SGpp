//
// Created by Vincent_Bode on 11.08.2017.
//

#ifndef SGPP_MPIREQUESTPOOL_HPP
#define SGPP_MPIREQUESTPOOL_HPP

#include <sgpp/datadriven/application/work_in_progress/NetworkMessageData.hpp>
#include <vector>
#include "PendingMPIRequest.h"

namespace sgpp {
    namespace datadriven {
        class MPIRequestPool {

        public:
            size_t createMPIRequestHandle();

            MPI_Request *getMPIRequestHandle(size_t handleIndex);

            void deleteMPIRequestHandle(size_t handleIndex);

            MPI_Request *getMPIRequests();

        protected:
            std::vector<MPI_Request> mpiRequestStorage;

        };

    }
}


#endif //SGPP_MPIREQUESTPOOL_HPP
