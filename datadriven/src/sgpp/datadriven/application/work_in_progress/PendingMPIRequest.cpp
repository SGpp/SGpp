//
// Created by Vincent_Bode on 11.08.2017.
//

#include "PendingMPIRequest.h"

namespace sgpp {
    namespace datadriven {
        PendingMPIRequest::PendingMPIRequest(MPIRequestPool &mpiRequestPool)
                : mpiRequestPool(mpiRequestPool) {
            mpiRequestIndex = mpiRequestPool.createMPIRequestHandle();
        }

        PendingMPIRequest::~PendingMPIRequest() {
            mpiRequestPool.deleteMPIRequestHandle(0);
        }

        MPI_Request *PendingMPIRequest::getMPIRequestHandle() {
            return mpiRequestPool.getMPIRequestHandle(mpiRequestIndex);
        }

    }
}