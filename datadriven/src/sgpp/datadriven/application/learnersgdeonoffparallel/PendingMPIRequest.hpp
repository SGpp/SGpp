// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/globaldef.hpp>

#include <sgpp/datadriven/application/learnersgdeonoffparallel/MPIRequestPool.hpp>

namespace sgpp {
namespace datadriven {
class PendingMPIRequest {
 public:
  /**
   * Create a PendingMPIRequest in the specified MPIRequestPool
   * @param requestPool The MPIRequestPool to create the request in.
   */
  explicit PendingMPIRequest(sgpp::datadriven::MPIRequestPool *requestPool);

  /**
   * Removes the PendingMPIRequest from its MPIRequestPool automatically.
   */
  ~PendingMPIRequest();

  /**
   * The buffer that contains the actual data to be sent or received.
   */
  sgpp::datadriven::MPI_Packet *buffer;
  /**
   * The callback to execute when the incoming/outgoing MPIRequest is completed.
   */
  std::function<void(PendingMPIRequest &)> callback;
  /**
   * Whether to delete the request after calling the callback. This is set to false for requests
   * that reuse their buffers after completing, such as restarting a read request in the same slot.
   */
  bool disposeAfterCallback;
  /**
   * Whether this PendingMPIRequest is a request for incoming messages.
   */
  bool inbound;

  /**
   * Fetch the MPI_Request from its handle that is attached to this PendingMPIRequest.
   * @return The MPI_Request stored in the MPIRequestPool.
   */
  MPI_Request *getMPIRequestFromHandle();

  /**
   * Get the Handle to the MPIRequestPool that can be used to get the associated MPI_Request.
   * @return The handle to the MPI_Request in the MPIRequestPool.
   */
  size_t getMPIRequestIndex();

 protected:
  /**
   * A reference to the MPIRequestPool to which this PendingMPIRequest is associated with.
   */
  sgpp::datadriven::MPIRequestPool &mpiRequestPool;
  /**
   * The handle to the MPI_Request inside the request pool.
   */
  size_t mpiRequestIndex;
};
}  // namespace datadriven
}  // namespace sgpp
