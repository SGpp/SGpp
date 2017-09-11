// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef SGPP_MPIREQUESTPOOL_HPP
#define SGPP_MPIREQUESTPOOL_HPP

#include <sgpp/datadriven/application/work_in_progress/NetworkMessageData.hpp>
#include <vector>
#include <set>

namespace sgpp {
namespace datadriven {
class MPIRequestPool {
 public:
  size_t createMPIRequestHandle();

  MPI_Request *getMPIRequestHandle(size_t handleIndex);

  void deleteMPIRequestHandle(size_t handleIndex);

  MPI_Request *getMPIRequests();

  size_t size();

 protected:
  std::set<size_t> freedRequests;
  std::vector<MPI_Request> mpiRequestStorage;

  void printPoolStatistics() const;
};
}  // namespace datadriven
}  // namespace sgpp

#endif  // SGPP_MPIREQUESTPOOL_HPP
