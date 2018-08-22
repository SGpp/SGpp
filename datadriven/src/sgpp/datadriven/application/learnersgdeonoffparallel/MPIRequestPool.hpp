// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/globaldef.hpp>

#include <sgpp/datadriven/application/learnersgdeonoffparallel/NetworkMessageData.hpp>
#include <vector>
#include <set>

namespace sgpp {
namespace datadriven {
class MPIRequestPool {
 public:
  /**
   * Allocate storage for one MPI_Request and return the handle to it.
   * @return The handle to the newly allocated MPI_Request.
   */
  size_t createMPIRequestHandle();

  /**
   * Gets the actual MPI_Request from its handle.
   * @param handleIndex The handle to the MPI_Request.
   * @return The MPI_Request itself.
   */
  MPI_Request *getMPIRequestHandle(size_t handleIndex);

  /**
   * Deallocate one MPI_Request based on its handle.
   * @param handleIndex The handle for which to delete.
   */
  void deleteMPIRequestHandle(size_t handleIndex);

  /**
   * Fetches the start of the MPI_Requests held sequentially in memory for use in MPI operations.
   * @return The start in memory of the MPI_Requests storage.
   */
  MPI_Request *getMPIRequests();

  /**
   * Gets the number of stored MPI_Requests currently held in the pool.
   * @return The number of MPI_Requests
   */
  size_t size();

 protected:
  /**
   * Holds indices to requests that were freed but could not yet be deallocated for later use
   * or deallocation.
   */
  std::set<size_t> freedRequests;
  /**
   * Holds the actual MPI_Requests sequentially in memory.
   */
  std::vector<MPI_Request> mpiRequestStorage;

  /**
   * Prints grid pool usage statistics for performance analysis.
   */
  void printPoolStatistics() const;
};
}  // namespace datadriven
}  // namespace sgpp
