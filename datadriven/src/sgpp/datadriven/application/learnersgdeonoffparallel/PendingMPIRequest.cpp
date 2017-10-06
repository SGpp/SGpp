// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/datadriven/application/learnersgdeonoffparallel/PendingMPIRequest.hpp>

namespace sgpp {
namespace datadriven {
PendingMPIRequest::PendingMPIRequest(MPIRequestPool *requestPool)
    : mpiRequestPool(*requestPool) {
  mpiRequestIndex = mpiRequestPool.createMPIRequestHandle();
}

PendingMPIRequest::~PendingMPIRequest() {
  mpiRequestPool.deleteMPIRequestHandle(mpiRequestIndex);
}

MPI_Request *PendingMPIRequest::getMPIRequestFromHandle() {
  return mpiRequestPool.getMPIRequestHandle(mpiRequestIndex);
}

size_t PendingMPIRequest::getMPIRequestIndex() {
  return mpiRequestIndex;
}
}  // namespace datadriven
}  // namespace sgpp
