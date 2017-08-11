//
// Created by Vincent_Bode on 11.08.2017.
//

#include <sgpp/datadriven/application/work_in_progress/PendingMPIRequest.hpp>

namespace sgpp {
    namespace datadriven {
        PendingMPIRequest::PendingMPIRequest(MPIRequestPool *requestPool)
                : mpiRequestPool(*requestPool) {
            mpiRequestIndex = mpiRequestPool.createMPIRequestHandle();
        }

        PendingMPIRequest::~PendingMPIRequest() {
            mpiRequestPool.deleteMPIRequestHandle(mpiRequestIndex);
        }

        MPI_Request *PendingMPIRequest::getMPIRequestHandle() {
            return mpiRequestPool.getMPIRequestHandle(mpiRequestIndex);
        }

        size_t PendingMPIRequest::getMPIRequestIndex() {
            return mpiRequestIndex;
        }

    }
}