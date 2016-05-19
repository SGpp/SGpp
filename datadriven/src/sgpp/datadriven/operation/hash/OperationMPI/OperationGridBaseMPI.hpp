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
/// Base class for grid releated mpi operations
class OperationGridMethod : public MPIOperation {
 private:
  void send_grid(base::Grid &grid) {
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

    // Send grid to slaves
    for (int dest = 1; dest < MPIEnviroment::get_node_count(); dest++)
      MPI_Send(gridpoints, static_cast<int>(gridsize * 2 * dimensions), MPI_INT,
               dest, 1, MPI_COMM_WORLD);
    // Send grid dimension to slaves
    for (int dest = 1; dest < MPIEnviroment::get_node_count(); dest++)
      MPI_Send(&dimensions, 1, MPI_INT, dest, 1, MPI_COMM_WORLD);

    delete [] gridpoints;
  }

 protected:
  size_t gridsize;

  /// Protected constructor - creates slave operations and sends the grid to them
  OperationGridMethod(base::Grid &grid, std::string slave_id)
      : MPIOperation(slave_id) {
    std::cerr << "ctor grid method" << std::endl;
    send_grid(grid);
  }
  /// Protected constructor - expects slaves already exist and just sends the grid to them
  explicit OperationGridMethod(base::Grid &grid) : MPIOperation() {
    send_grid(grid);
  }

 protected:
  /// Abstract base class for grid related slave operations
  class OperationGridMethodSlave : virtual public MPISlaveOperation {
   protected:
    bool verbose;
    int grid_dimensions;
    int complete_gridsize;
    int *gridpoints;
    MPI_Status stat;

   public:
    OperationGridMethodSlave() : verbose(true), complete_gridsize(0) {
      // receive grid
      MPI_Probe(0, 1, MPI_COMM_WORLD, &stat);
      MPI_Get_count(&stat, MPI_INT, &complete_gridsize);
      gridpoints = new int[complete_gridsize];
      MPI_Recv(gridpoints, complete_gridsize, MPI_INT, stat.MPI_SOURCE, stat.MPI_TAG,
               MPI_COMM_WORLD, &stat);

      // Receive grid dimensions
      MPI_Probe(0, 1, MPI_COMM_WORLD, &stat);
      MPI_Recv(&grid_dimensions, 1, MPI_INT, stat.MPI_SOURCE, stat.MPI_TAG,
               MPI_COMM_WORLD, &stat);
      if (verbose) {
        std::cout << "Node " << MPIEnviroment::get_node_rank()
                  << ":  Recevied grid with "<< complete_gridsize / (grid_dimensions *2)
                  << " integers and " << grid_dimensions << " dimensions " << std::endl;
      }
    }
    virtual ~OperationGridMethodSlave() {
      delete [] gridpoints;
    }
    virtual void slave_code() = 0;
  };

 public:
  virtual ~OperationGridMethod() {}
};
}  // namespace clusteringmpi
}  // namespace datadriven
}  // namespace sgpp
#endif /* OPERATIONGRIDBASEMPI_H */
