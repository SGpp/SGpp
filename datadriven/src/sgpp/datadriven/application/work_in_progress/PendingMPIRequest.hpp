//
// Created by Vincent_Bode on 11.08.2017.
//

#ifndef SGPP_PENDINGMPIREQUEST_HPP
#define SGPP_PENDINGMPIREQUEST_HPP


#include <sgpp/datadriven/application/work_in_progress/MPIRequestPool.hpp>

namespace sgpp {
    namespace datadriven {
        class PendingMPIRequest {
        public:
            explicit PendingMPIRequest(sgpp::datadriven::MPIRequestPool *requestPool);

            ~PendingMPIRequest();

            sgpp::datadriven::MPI_Packet *buffer;
            std::function<void(PendingMPIRequest &)> callback;
            bool disposeAfterCallback;

            MPI_Request *getMPIRequestHandle();

            size_t getMPIRequestIndex();

        protected:
            sgpp::datadriven::MPIRequestPool &mpiRequestPool;
            size_t mpiRequestIndex;

        };

    }
}


#endif //SGPP_PENDINGMPIREQUEST_HPP
