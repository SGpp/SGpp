// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org
#ifndef OPERATIONGRIDBASEMPI_H
#define OPERATIONGRIDBASEMPI_H

#include <mpi.h>

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/datadriven/operation/hash/OperationMPI/MPIEnviroment.hpp>
#include <sgpp/datadriven/operation/hash/OperationMPI/OperationMPI.hpp>

#include <string>

namespace sgpp {
namespace datadriven {
namespace clusteringmpi {
class GridCommunicator {
 private:
  MPI_Comm &communicator;
  int worker_count;

 public:
  explicit GridCommunicator(MPI_Comm &comm, int worker_count)
      : communicator(comm), worker_count(worker_count) {}

  void send_grid(int *gridpoints, int dimensions, int gridsize) {
    // Send grid to slaves
    for (int dest = 1; dest < worker_count; dest++)
      MPI_Send(gridpoints, static_cast<int>(gridsize * 2 * dimensions), MPI_INT,
               dest, 1, communicator);
    // Send grid dimension to slaves
    for (int dest = 1; dest < worker_count; dest++)
      MPI_Send(&dimensions, 1, MPI_INT, dest, 1, communicator);
  }
  void receive_grid(int *gridpoints, int &dimensions, int &gridsize) {
    MPI_Status stat;
    int complete_gridsize;
    // receive grid
    MPI_Probe(0, 1, communicator, &stat);
    MPI_Get_count(&stat, MPI_INT, &complete_gridsize);
    gridpoints = new int[complete_gridsize];
    MPI_Recv(gridpoints, complete_gridsize, MPI_INT, stat.MPI_SOURCE, stat.MPI_TAG,
             communicator, &stat);

    // Receive grid dimensions
    MPI_Probe(0, 1, communicator, &stat);
    MPI_Recv(&dimensions, 1, MPI_INT, stat.MPI_SOURCE, stat.MPI_TAG,
             communicator, &stat);

    gridsize = complete_gridsize / (2 * dimensions);
  }
};

/// Base class for grid releated mpi operations
class OperationGridMethod : public MPIOperation {
 protected:
  size_t gridsize;
  GridCommunicator *grid_comm;

 public:
  /// Protected constructor - creates slave operations and sends the grid to them
  OperationGridMethod(base::OperationConfiguration conf, base::Grid &grid, std::string slave_id)
      : MPIOperation(conf, slave_id) {
    std::cerr << "ctor grid method" << std::endl;
    grid_comm = new GridCommunicator(node_comm->get_communicator(), node_comm->get_workercount());

    // Store grid in integer array
    sgpp::base::GridStorage& gridStorage = grid.getStorage();
    gridsize = gridStorage.getSize();
    size_t dimensions = gridStorage.getDimension();
    int *gridpoints = new int[gridsize * 2 * dimensions];
    size_t pointscount = 0;
    for (size_t i = 0; i < gridsize; i++) {
      sgpp::base::HashGridIndex *point = gridStorage.get(i);
      pointscount++;
      for (size_t d = 0; d < dimensions; d++) {
        gridpoints[i * 2 * dimensions + 2 * d] = point->getIndex(d);
        gridpoints[i * 2 * dimensions + 2 * d + 1] = point->getLevel(d);
      }
    }
    grid_comm->send_grid(gridpoints, dimensions, gridsize / (2 * dimensions));
    delete [] gridpoints;
  }
};
class MPIGridWorker : public MPIWorkerBase {
 protected:
  int *gridpoints;
  int gridsize;
  int dimensions;
  GridCommunicator *grid_comm;
 public:
  MPIGridWorker(int masternode, std::string operationName,
                base::OperationConfiguration conf)
      : MPIWorkerBase(masternode, operationName, conf), gridpoints(NULL) {
    grid_comm = new GridCommunicator(node_comm->get_communicator(), node_comm->get_workercount());
    grid_comm->receive_grid(gridpoints, dimensions, gridsize);
  }
  virtual ~MPIGridWorker() {
    if (gridpoints != NULL)
      delete [] gridpoints;
    delete grid_comm;
  }
  virtual void slave_code() {}
};

}  // namespace clusteringmpi
}  // namespace datadriven
}  // namespace sgpp
#endif /* OPERATIONGRIDBASEMPI_H */
