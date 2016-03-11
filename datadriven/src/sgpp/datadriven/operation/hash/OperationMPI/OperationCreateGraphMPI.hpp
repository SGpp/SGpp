// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org
#ifndef OPERATIONCREATEGRAPHMPI_H
#define OPERATIONCREATEGRAPHMPI_H

#include <mpi.h>
#include <sgpp/globaldef.hpp>
#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/datadriven/operation/hash/OperationMPI/MPIEnviroment.hpp>
#include <sgpp/datadriven/operation/hash/OperationMPI/OperationMPI.hpp>
#include <sgpp/datadriven/operation/hash/OperationCreateGraphOCL/OpFactory.hpp>
#include <sgpp/datadriven/operation/hash/OperationCreateGraphOCL/OperationCreateGraphOCL.hpp>

#include <vector>

namespace sgpp {
namespace datadriven {
namespace clusteringmpi {
class OperationCreateGraphSlave : public MPISlaveOperation {
 public:
  bool verbose;
  OperationCreateGraphSlave()
      : MPISlaveOperation() {
    verbose = true;
  }
  virtual ~OperationCreateGraphSlave() {}
  virtual void slave_code(void) {
    MPI_Status stat;
    int k;
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
    MPI_Probe(0, 1, MPI_COMM_WORLD, &stat);
    MPI_Recv(&k, 1, MPI_INT, stat.MPI_SOURCE, stat.MPI_TAG,
             MPI_COMM_WORLD, &stat);
    if (verbose) {
      std::cout << "Node " << MPIEnviroment::get_node_rank() << ": Received dataset (size "
                << dataset_size / dimensions << " datapoints) and graph parameters" << std::endl;
    }
    // Create opencl operation
    sgpp::datadriven::DensityOCLMultiPlatform::OperationCreateGraphOCL* op =
        sgpp::datadriven::createNearestNeighborGraphConfigured(dataset, dataset_size,
                                                               k, dimensions, "MyOCLConf.cfg");

    if (verbose) {
      std::cout << "Node " << MPIEnviroment::get_node_rank()
                << ": Created opencl graph operation"<< std::endl;
    }
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
        std::vector<int> partial_graph(datainfo[1] * k);
        op->create_graph(partial_graph, datainfo[0], datainfo[1]);
        // Send results back
        MPI_Send(partial_graph.data(), datainfo[1] * k, MPI_INT, 0, 1, MPI_COMM_WORLD);
      }
    }while(true);
  }
};

class OperationCreateGraphMPI : public MPIOperation {
 private:
  sgpp::base::DataMatrix &data;
  size_t datasize;
  size_t dimensions;
  size_t k;

 public:
  OperationCreateGraphMPI(sgpp::base::DataMatrix& data, size_t dimensions, size_t k)
      : MPIOperation(typeid(OperationCreateGraphSlave).name()), data(data),
        datasize(data.getSize()), dimensions(dimensions), k(k) {
  }
  std::vector<int> create_graph(void) {
    this->start_slave_code();
    int *graph = new int[datasize / dimensions * k];
    std::vector<int> returngraph(graph, graph + (datasize / dimensions *k));

    // Sending dataset to slaves
    for (int dest = 1; dest < MPIEnviroment::get_node_count(); dest++)
      MPI_Send(data.getPointer(), static_cast<int>(data.getSize()), MPI_DOUBLE,
               dest, 1, MPI_COMM_WORLD);
    for (int dest = 1; dest < MPIEnviroment::get_node_count(); dest++)
      MPI_Send(&dimensions, 1, MPI_INT, dest, 1, MPI_COMM_WORLD);
    for (int dest = 1; dest < MPIEnviroment::get_node_count(); dest++)
      MPI_Send(&k, 1, MPI_INT, dest, 1, MPI_COMM_WORLD);

    // Create packages and let the slaves solve them
    int *partial_result = new int[2000 * k];
    SimpleQueue<int> workitem_queue(datasize / dimensions, 2000);
    int chunkid = 0;
    int messagesize = workitem_queue.receive_result(chunkid, partial_result);
    while (messagesize > 0) {
      // Store result
      for (int i = 0; i < messagesize; i++) {
        returngraph[chunkid*k+i] = partial_result[i];
      }
      messagesize = workitem_queue.receive_result(chunkid, partial_result);
    }
    int exitmessage[2] = {-2, -2};
    for (int dest = 1; dest < MPIEnviroment::get_node_count(); dest++)
      MPI_Send(exitmessage, 2, MPI_INT, dest, 1, MPI_COMM_WORLD);
    delete [] partial_result;
    delete [] graph;
    return returngraph;
  }
  virtual ~OperationCreateGraphMPI() {}
};

}  // namespace clusteringmpi
}  // namespace datadriven
}  // namespace sgpp
#endif /* OPERATIONCREATEGRAPHMPI_H */
