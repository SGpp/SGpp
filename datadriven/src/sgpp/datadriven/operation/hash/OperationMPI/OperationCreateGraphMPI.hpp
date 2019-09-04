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

class GraphCreationWorker : public MPIWorkerGraphBase {
 protected:
  bool opencl_node;
  bool overseer_node;
  double lambda;
  DensityOCLMultiPlatform::OperationCreateGraphOCL *op;

  MPI_Comm &master_worker_comm;
  MPI_Comm &sub_worker_comm;
  void divide_workpackages(int *package, std::vector<int> &graph) {
    // Divide into more work packages
    int packagesize = static_cast<int>(MPIEnviroment::get_configuration()
                                       ["PREFERED_PACKAGESIZE"].getInt());
    int *partial_result = new int[package[1] * k];
    SimpleQueue<int> workitem_queue(package[0], package[1], packagesize,
                                    sub_worker_comm,
                                    MPIEnviroment::get_sub_worker_count());
    int chunkid = package[0];
    size_t messagesize = 0;
    while (!workitem_queue.is_finished()) {
      // Store result
      if (verbose) {
        std::cout << "Messagesize: "<< messagesize << std::endl;
        std::cout << partial_result[0] << " at  " << chunkid - package[0] + 0
                  << " with packageid " << chunkid << " on "
                  << MPIEnviroment::get_node_rank() << std::endl;
      }
      messagesize = workitem_queue.receive_result(chunkid, partial_result);
      for (size_t i = 0; i < messagesize; i++) {
        graph[chunkid - package[0] + i] = partial_result[i];
      }
    }
    delete [] partial_result;
  }

 public:
  GraphCreationWorker()
      : MPIWorkerBase("GraphCreationWorker"),
        MPIWorkerGraphBase("GraphCreationWorker"),
        opencl_node(false), overseer_node(false),
        master_worker_comm(MPIEnviroment::get_input_communicator()),
        sub_worker_comm(MPIEnviroment::get_communicator()) {
    if (MPIEnviroment::get_sub_worker_count() > 0) {
      overseer_node = true;
      opencl_node = false;
    } else {
      overseer_node = false;
      opencl_node = true;
    }
    // Create opencl operation
    if (opencl_node) {
      op = createNearestNeighborGraphConfigured(dataset, dataset_size, k, dimensions,
                                                "MyOCLConf.cfg", 0, 0);
    }
  }
  GraphCreationWorker(sgpp::base::DataMatrix &data, int k)
      : MPIWorkerBase("GraphCreationWorker"),
        MPIWorkerGraphBase("GraphCreationWorker", data, k),
        opencl_node(false), overseer_node(false),
        master_worker_comm(MPIEnviroment::get_input_communicator()),
        sub_worker_comm(MPIEnviroment::get_communicator()) {}
  virtual ~GraphCreationWorker(void) {
    if (opencl_node) {
      delete op;
    }
  }
  void start_worker_main(void) {
    base::DataMatrix data_matrix(dataset, dataset_size / dimensions, dimensions);
    // Receive alpha
    MPI_Status stat;

    // Work Loop
    int datainfo[2];
    std::vector<int> partial_graph;
    int old_partial_size = 0;
    do {
      // Receive Workpackage
      MPI_Probe(0, 1, master_worker_comm, &stat);
      MPI_Recv(datainfo, 2, MPI_INT, 0, stat.MPI_TAG, master_worker_comm, &stat);
      if (verbose) {
        std::cout << "Received workpackage [" << datainfo[0] << "," << datainfo[1]
                  << "] on " << MPIEnviroment::get_node_rank() << std::endl;
      }
      // Check for exit
      if (datainfo[0] == -2 && datainfo[1] == -2) {
        std::cerr << "Node" << MPIEnviroment::get_node_rank()
                  << " received exit signal" << std::endl;
        for (int dest = 1; dest < MPIEnviroment::get_sub_worker_count() + 1; dest++)
          MPI_Send(datainfo, 2, MPI_INT, dest, 1, sub_worker_comm);
        break;
      } else {
        if (datainfo[1] * k != old_partial_size) {
          partial_graph.resize(datainfo[1] * k);
          if (verbose)
            std::cout << "New Buffer created!" << std::endl;
          old_partial_size = datainfo[1] * k;
        }
        if (opencl_node) {
          // Create partial rhs
          op->create_graph(partial_graph, datainfo[0], datainfo[1]);
          if (verbose)
            std::cerr << "Workpackage abgeschlossen" << std::endl;
        } else {
          divide_workpackages(datainfo, partial_graph);
        }
        // Send results back
        MPI_Send(partial_graph.data(), datainfo[1] * k, MPI_INT, 0, 1, master_worker_comm);
      }
    }while(true);
  }
};
class OperationGraphCreationMPI : public GraphCreationWorker {
 public:
  OperationGraphCreationMPI(sgpp::base::DataMatrix &data, int k) :
      MPIWorkerBase("GraphCreationWorker"), GraphCreationWorker(data, k) {
  }
  virtual ~OperationGraphCreationMPI() {}
  virtual void create_graph(std::vector<int> &result) {
    start_sub_workers();
    int datainfo[2];
    datainfo[0] = 0;
    datainfo[1] = dataset_size / dimensions;

    result.resize(dataset_size / dimensions * k);
    divide_workpackages(datainfo, result);
    int exitmessage[2] = {-2, -2};
    for (int dest = 1; dest < MPIEnviroment::get_sub_worker_count() + 1; dest++)
      MPI_Send(exitmessage, 2, MPI_INT, dest, 1, sub_worker_comm);
  }
};
}  // namespace clusteringmpi
}  // namespace datadriven
}  // namespace sgpp
#endif /* OPERATIONCREATEGRAPHMPI_H */
