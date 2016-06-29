// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org
#ifndef OPERATIONCREATEPRUNEDGRAPHMPI_H
#define OPERATIONCREATEPRUNEDGRAPHMPI_H

#include <sgpp/datadriven/operation/hash/OperationMPI/OperationGridBaseMPI.hpp>
#include <sgpp/datadriven/operation/hash/OperationMPI/OperationGraphBaseMPI.hpp>
#include <sgpp/datadriven/operation/hash/OperationDensityOCLMultiPlatform/OpFactory.hpp>
#include <sgpp/datadriven/operation/hash/OperationCreateGraphOCL/OpFactory.hpp>
#include <sgpp/datadriven/operation/hash/OperationCreateGraphOCL/OperationCreateGraphOCL.hpp>
#include <sgpp/datadriven/operation/hash/OperationPruneGraphOCL/OpFactory.hpp>
#include <sgpp/datadriven/operation/hash/OperationPruneGraphOCL/OperationPruneGraphOCL.hpp>
#include <string>
#include <vector>

namespace sgpp {
namespace datadriven {
namespace clusteringmpi {
class OperationCreatePrunedGraph : public OperationGridMethod, public OperationGraphMethodMPI {
 public:
  OperationCreatePrunedGraph(base::OperationConfiguration conf, base::Grid& grid,
                             base::DataVector &alpha, base::DataMatrix &data,
                             size_t k, int packagesize)
      : OperationGridMethod(conf, grid, "OperationCreatePrunedGraphSlave"),
        OperationGraphMethodMPI(data, grid.getDimension(), k), packagesize(packagesize) {
    // Send alpha vector
    for (int dest = 1; dest < MPIEnviroment::get_node_count(); dest++)
      MPI_Send(alpha.getPointer(), static_cast<int>(gridsize), MPI_DOUBLE, dest, 1, MPI_COMM_WORLD);
  }
  std::vector<int> createPrunedGraph(double treshold) {
    OperationGridMethod::start_slave_code();
    // Send treshold to slaves
    std::vector<int> ret_graph(datasize / dimensions * k);
    for (int dest = 1; dest < MPIEnviroment::get_node_count(); dest++)
      MPI_Send(&treshold, 1, MPI_DOUBLE, dest, 1, MPI_COMM_WORLD);

    int *partial_result = new int[packagesize * k];
    SimpleQueue<int> workitem_queue(datasize / dimensions, packagesize);
    int chunkid = 0;
    size_t messagesize = workitem_queue.receive_result(chunkid, partial_result);
    while (messagesize > 0) {
      // Store result
      for (size_t i = 0; i < messagesize; i++) {
        ret_graph[chunkid*k+i] = partial_result[i];
      }
      messagesize = workitem_queue.receive_result(chunkid, partial_result);
    }
    int exitmessage[2] = {-2, -2};
    for (int dest = 1; dest < MPIEnviroment::get_node_count(); dest++)
      MPI_Send(exitmessage, 2, MPI_INT, dest, 1, MPI_COMM_WORLD);
    delete [] partial_result;

    return ret_graph;
  }

  static MPISlaveOperation* create_slave(base::OperationConfiguration conf) {
    return new OperationCreatePrunedGraphSlave(conf);
  }
  virtual ~OperationCreatePrunedGraph() {}

 private:
  int packagesize;
  class OperationCreatePrunedGraphSlave  : public OperationGridMethodSlave,
                                           public OperationGraphMethodSlave {
   private:
    double *alpha;
   public:
    bool verbose;
    explicit OperationCreatePrunedGraphSlave(base::OperationConfiguration conf) :
        OperationGridMethodSlave(conf),
        OperationGraphMethodSlave(conf) {
      verbose = true;
      // Receive alpha vector
      int buffer_size = 0;
      MPI_Probe(0, 1, MPI_COMM_WORLD, &stat);
      MPI_Get_count(&stat, MPI_DOUBLE, &buffer_size);
      if (buffer_size != complete_gridsize / (2 * grid_dimensions)) {
        std::stringstream errorString;
        errorString << "Error: Gridsize " << complete_gridsize / (2 * grid_dimensions)
                    << " and the size of the alpha vector "
                    << buffer_size << " should match!" << std::endl;
        throw std::logic_error(errorString.str());
      }
      alpha = new double[complete_gridsize / (2 * grid_dimensions)];
      MPI_Recv(alpha, complete_gridsize / (2 * grid_dimensions), MPI_DOUBLE, stat.MPI_SOURCE,
               stat.MPI_TAG, MPI_COMM_WORLD, &stat);
    }
    virtual ~OperationCreatePrunedGraphSlave() {
      delete [] alpha;
    }

    virtual void slave_code(void) {
      double treshold = 0.0;
      MPI_Probe(0, 1, MPI_COMM_WORLD, &stat);
      MPI_Recv(&treshold, 1, MPI_DOUBLE, stat.MPI_SOURCE,
               stat.MPI_TAG, MPI_COMM_WORLD, &stat);

      base::DataMatrix data_matrix(dataset, dataset_size / dimensions, dimensions);
      // Create opencl operation
      sgpp::datadriven::DensityOCLMultiPlatform::OperationCreateGraphOCL* op_create =
          sgpp::datadriven::createNearestNeighborGraphConfigured(dataset, dataset_size,
                                                                 k, dimensions, "MyOCLConf.cfg",
                                                                 0, 0);
      sgpp::datadriven::DensityOCLMultiPlatform::OperationPruneGraphOCL* op_prune =
          sgpp::datadriven::pruneNearestNeighborGraphConfigured(gridpoints, complete_gridsize /
                                                                (2 * grid_dimensions),
                                                                grid_dimensions, alpha,
                                                                data_matrix, treshold, k,
                                                                "MyOCLConf.cfg", 0, 0);
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
          op_create->create_graph(partial_graph, datainfo[0], datainfo[1]);
          op_prune->prune_graph(partial_graph, datainfo[0], datainfo[1]);
          // Send results back
          MPI_Send(partial_graph.data(), datainfo[1], MPI_INT, 0, 1, MPI_COMM_WORLD);
        }
      }while(true);
      delete op_create;
      delete op_prune;
    }
  };
};

}  // namespace clusteringmpi
}  // namespace datadriven
}  // namespace sgpp
#endif /* OPERATIONCREATEPRUNEDGRAPHMPI_H */
