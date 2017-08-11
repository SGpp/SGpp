//
// Created by Vincent_Bode on 11.08.2017.
//

#ifndef SGPP_PENDINGMPIREQUEST_H
#define SGPP_PENDINGMPIREQUEST_H


#include "NetworkMessageData.hpp"
#include "MPIRequestPool.hpp"

namespace sgpp {
    namespace datadriven {
        class PendingMPIRequest {
        public:
            explicit PendingMPIRequest(MPIRequestPool &mpiRequestPool);

            ~PendingMPIRequest();

            sgpp::datadriven::MPI_Packet *buffer;
            std::function<void(PendingMPIRequest &)> callback;
            bool disposeAfterCallback;

            MPI_Request *getMPIRequestHandle();

        protected:
            MPIRequestPool &mpiRequestPool;
            size_t mpiRequestIndex;

        };

    }
}


#endif //SGPP_PENDINGMPIREQUEST_H
