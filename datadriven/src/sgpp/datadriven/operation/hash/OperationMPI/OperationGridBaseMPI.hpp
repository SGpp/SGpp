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

class MPIWorkerGridBase : virtual public MPIWorkerBase {
 private:
  void receive_grid(void) {
    // receive grid
    MPI_Status stat;
    MPI_Probe(0, 1, MPIEnviroment::get_input_communicator(), &stat);
    MPI_Get_count(&stat, MPI_INT, &complete_gridsize);
    gridpoints = new int[complete_gridsize];
    MPI_Recv(gridpoints, complete_gridsize, MPI_INT, stat.MPI_SOURCE, stat.MPI_TAG,
             MPIEnviroment::get_input_communicator(), &stat);

    // Receive grid dimensions
    MPI_Probe(0, 1, MPIEnviroment::get_input_communicator(), &stat);

    MPI_Recv(&grid_dimensions, 1, MPI_INT, stat.MPI_SOURCE, stat.MPI_TAG,
             MPIEnviroment::get_input_communicator(), &stat);
    gridsize = complete_gridsize / (2 * grid_dimensions);
    if (verbose) {
      std::cout << "Node " << MPIEnviroment::get_node_rank() << ":  Recevied grid with "
                << complete_gridsize / (grid_dimensions * 2) << " integers and " << grid_dimensions
                << " dimensions " << std::endl;
    }
  }
  void send_grid(void) {
    // Send grid to slaves
    for (int i = 1; i < MPIEnviroment::get_sub_worker_count() + 1; i++) {
      MPI_Send(gridpoints, static_cast<int>(complete_gridsize), MPI_INT, i, 1,
               MPIEnviroment::get_communicator());
    }
    // Send grid dimension to slaves
    for (int i = 1; i < MPIEnviroment::get_sub_worker_count() + 1; i++) {
      MPI_Send(&grid_dimensions, 1, MPI_INT, i, 1, MPIEnviroment::get_communicator());
    }
  }

 protected:
  int grid_dimensions;
  int complete_gridsize;
  int gridsize;
  int *gridpoints;
  explicit MPIWorkerGridBase(std::string operationName) : MPIWorkerBase(operationName) {
    receive_grid();
    send_grid();
  }
  MPIWorkerGridBase(std::string operationName, base::Grid &grid) : MPIWorkerBase(operationName) {
    // Store grid in integer array
    std::cout << "IN GridWorker cstr"
              << "\n";
    sgpp::base::GridStorage &gridStorage = grid.getStorage();
    gridsize = static_cast<int>(gridStorage.getSize());
    int dimensions = static_cast<int>(gridStorage.getDimension());
    gridpoints = new int[gridsize * 2 * dimensions];
    size_t pointscount = 0;
    for (int i = 0; i < gridsize; i++) {
      sgpp::base::HashGridPoint &point = gridStorage.getPoint(i);
      pointscount++;
      for (int d = 0; d < dimensions; d++) {
        gridpoints[i * 2 * dimensions + 2 * d] = point.getIndex(d);
        gridpoints[i * 2 * dimensions + 2 * d + 1] = point.getLevel(d);
      }
    }

    complete_gridsize = gridsize * 2 * dimensions;
    grid_dimensions = dimensions;
    send_grid();
  }

 public:
  virtual ~MPIWorkerGridBase() { delete[] gridpoints; }
};

}  // namespace clusteringmpi
}  // namespace datadriven
}  // namespace sgpp
#endif /* OPERATIONGRIDBASEMPI_H */
