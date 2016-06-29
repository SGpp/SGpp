// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org
#ifndef OPERATIONCREATEGRAPHMPI_H
#define OPERATIONCREATEGRAPHMPI_H

#include <sgpp/datadriven/operation/hash/OperationMPI/OperationGraphBaseMPI.hpp>
#include <sgpp/datadriven/operation/hash/OperationCreateGraphOCL/OpFactory.hpp>
#include <sgpp/datadriven/operation/hash/OperationCreateGraphOCL/OperationCreateGraphOCL.hpp>

#include <vector>
#include <string>

namespace sgpp {
namespace datadriven {
namespace clusteringmpi {

class OperationCreateGraphMPI : public OperationGraphMethodMPI {
 public:
  static MPISlaveOperation* create_slave(base::OperationConfiguration conf) {
    return new OperationCreateGraphSlave(conf);
  }
  OperationCreateGraphMPI(base::OperationConfiguration conf,
                          sgpp::base::DataMatrix& data,
                          size_t dimensions, size_t k, int packagesize)
      : OperationGraphMethodMPI(conf, data, dimensions, k, "OperationCreateGraphSlave"),
        packagesize(packagesize) {
  }
  std::vector<int> create_graph(void) {
    this->start_slave_code();
    std::vector<int> returngraph(datasize / dimensions *k);

    // Create packages and let the slaves solve them
    int *partial_result = new int[packagesize * k];
    SimpleQueue<int> workitem_queue(datasize / dimensions, packagesize);
    int chunkid = 0;
    size_t messagesize = workitem_queue.receive_result(chunkid, partial_result);
    while (messagesize > 0) {
      // Store result
      for (size_t i = 0; i < messagesize; i++) {
        returngraph[chunkid*k+i] = partial_result[i];
      }
      messagesize = workitem_queue.receive_result(chunkid, partial_result);
    }
    int exitmessage[2] = {-2, -2};
    for (int dest = 1; dest < MPIEnviroment::get_node_count(); dest++)
      MPI_Send(exitmessage, 2, MPI_INT, dest, 1, MPI_COMM_WORLD);
    delete [] partial_result;
    return returngraph;
  }
  virtual ~OperationCreateGraphMPI() {}

 protected:
  int packagesize;
  class OperationCreateGraphSlave : public OperationGraphMethodSlave {
   public:
    explicit OperationCreateGraphSlave(base::OperationConfiguration conf)
        : OperationGraphMethodSlave(conf) {
    }
    virtual ~OperationCreateGraphSlave() {}
    virtual void slave_code(void) {
      MPI_Status stat;
      // Create opencl operation
      sgpp::datadriven::DensityOCLMultiPlatform::OperationCreateGraphOCL* op =
          sgpp::datadriven::createNearestNeighborGraphConfigured(dataset, dataset_size,
                                                                 k, dimensions, "MyOCLConf.cfg",
                                                                 0, 0);
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
      delete op;
    }
  };
};

}  // namespace clusteringmpi
}  // namespace datadriven
}  // namespace sgpp
#endif /* OPERATIONCREATEGRAPHMPI_H */
