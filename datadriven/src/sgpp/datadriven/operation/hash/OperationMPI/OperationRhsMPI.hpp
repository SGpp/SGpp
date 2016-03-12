// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org
#ifndef OPERATIONRHSMPI_H
#define OPERATIONRHSMPI_H

#include <mpi.h>

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/datadriven/operation/hash/OperationMPI/MPIEnviroment.hpp>
#include <sgpp/datadriven/operation/hash/OperationMPI/OperationMPI.hpp>
#include <sgpp/datadriven/operation/hash/OperationDensityOCLMultiPlatform/OpFactory.hpp>

namespace sgpp {
namespace datadriven {
namespace clusteringmpi {

class OperationRhsSlave : public MPISlaveOperation {
 public:
  bool verbose;
  OperationRhsSlave() : MPISlaveOperation() {
    verbose = true;
  }
  virtual ~OperationRhsSlave() {}
  virtual void slave_code(void) {
    MPI_Status stat;
    int dimensions;
    // Receive dataset
    int dataset_size = 0;
    MPI_Probe(0, 1, MPI_COMM_WORLD, &stat);
    MPI_Get_count(&stat, MPI_DOUBLE, &dataset_size);
    double *dataset = new double[dataset_size];
    MPI_Recv(dataset, dataset_size, MPI_DOUBLE, stat.MPI_SOURCE, stat.MPI_TAG,
             MPI_COMM_WORLD, &stat);
    // Receive clustering parameters
    MPI_Probe(0, 1, MPI_COMM_WORLD, &stat);
    MPI_Recv(&dimensions, 1, MPI_INT, stat.MPI_SOURCE, stat.MPI_TAG,
             MPI_COMM_WORLD, &stat);
    base::DataMatrix data_matrix(dataset, dataset_size / dimensions, dimensions);
    if (verbose) {
      std::cout << "Node " << MPIEnviroment::get_node_rank() << ": Received dataset (size "
                << dataset_size / dimensions << " datapoints)" << std::endl;
    }

    // Receive grid
    int complete_gridsize = 0;
    MPI_Probe(0, 1, MPI_COMM_WORLD, &stat);
    MPI_Get_count(&stat, MPI_INT, &complete_gridsize);
    int *grid = new int[complete_gridsize];
    MPI_Recv(grid, complete_gridsize, MPI_INT, stat.MPI_SOURCE, stat.MPI_TAG,
             MPI_COMM_WORLD, &stat);

    DensityOCLMultiPlatform::OperationDensityOCL *op =
        sgpp::datadriven::createDensityOCLMultiPlatformConfigured(grid, complete_gridsize /
                                                                  (2 * dimensions), dimensions,
                                                                  0.0, "MyOCLConf.cfg");
    int datainfo[2];
    do {
      MPI_Probe(0, 1, MPI_COMM_WORLD, &stat);
      MPI_Recv(datainfo, 2, MPI_INT, 0, stat.MPI_TAG, MPI_COMM_WORLD, &stat);
      // Check for exit
      if (datainfo[0] == -2 && datainfo[1] == -2) {
        std::cerr << "Node" << MPIEnviroment::get_node_rank()
                  << " received exit signal" << std::endl;
        break;
      } else {
        if (verbose) {
          std::cout << "Node " << MPIEnviroment::get_node_rank()
                    << ": Received work package" << std::endl;
        }
        // Create partial graph
        base::DataVector partial_rhs(datainfo[1]);
        op->generateb(data_matrix, partial_rhs, datainfo[0], datainfo[1]);
        // Send results back
        MPI_Send(partial_rhs.getPointer(), datainfo[1], MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
      }
    }while(true);
    delete op;
    delete [] grid;
    delete [] dataset;
  }
};

class OperationRhsMPI : public MPIOperation {
 private:
  base::Grid &grid;
  size_t dimensions;
  sgpp::base::DataMatrix &data;
  size_t datasize;
 public:
  OperationRhsMPI(base::Grid& grid, size_t dimensions, base::DataMatrix &data)
      : MPIOperation(typeid(OperationRhsSlave).name()), grid(grid), dimensions(dimensions),
        data(data), datasize(data.getSize()) {}
  base::DataVector create_rhs(void) {
    this->start_slave_code();

    // Sending dataset to slaves
    base::DataVector ret_rhs(grid.getStorage().getSize());
    for (int dest = 1; dest < MPIEnviroment::get_node_count(); dest++)
      MPI_Send(data.getPointer(), static_cast<int>(data.getSize()), MPI_DOUBLE,
               dest, 1, MPI_COMM_WORLD);
    for (int dest = 1; dest < MPIEnviroment::get_node_count(); dest++)
      MPI_Send(&dimensions, 1, MPI_INT, dest, 1, MPI_COMM_WORLD);

    // Store grid in integer array
    sgpp::base::GridStorage& gridStorage = grid.getStorage();
    size_t gridsize = gridStorage.getSize();
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

    // Create packages and let the slaves solve them
    double *partial_result = new double[2000];
    SimpleQueue<double> workitem_queue(gridsize, 2000);
    int chunkid = 0;
    int messagesize = workitem_queue.receive_result(chunkid, partial_result);
    while (messagesize > 0) {
      // Store result
      std::cerr << messagesize << std::endl;
      for (int i = 0; i < messagesize; i++) {
        ret_rhs[chunkid + i] = partial_result[i];
      }
      messagesize = workitem_queue.receive_result(chunkid, partial_result);
    }
    int exitmessage[2] = {-2, -2};
    for (int dest = 1; dest < MPIEnviroment::get_node_count(); dest++)
      MPI_Send(exitmessage, 2, MPI_INT, dest, 1, MPI_COMM_WORLD);
    delete [] partial_result;
    return ret_rhs;
  }
  virtual ~OperationRhsMPI() {}
};


}  // namespace clusteringmpi
}  // namespace datadriven
}  // namespace sgpp
#endif /* OPERATIONRHSMPI_H */ 
