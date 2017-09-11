// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

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
  bool inbound;

  MPI_Request *getMPIRequestHandle();

  size_t getMPIRequestIndex();

 protected:
  sgpp::datadriven::MPIRequestPool &mpiRequestPool;
  size_t mpiRequestIndex;
};
}  // namespace datadriven
}  // namespace sgpp

#endif  // SGPP_PENDINGMPIREQUEST_HPP
