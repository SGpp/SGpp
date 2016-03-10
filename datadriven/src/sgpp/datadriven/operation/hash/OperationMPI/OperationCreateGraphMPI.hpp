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
    SGPP::datadriven::DensityOCLMultiPlatform::OperationCreateGraphOCL* op =
        SGPP::datadriven::createNearestNeighborGraphConfigured(dataset, dataset_size,
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
  SGPP::base::DataMatrix &data;
  size_t datasize;
  size_t dimensions;
  size_t k;

 public:
  OperationCreateGraphMPI(SGPP::base::DataMatrix& data, size_t dimensions, size_t k)
      : MPIOperation(typeid(OperationCreateGraphSlave).name()), data(data),
        datasize(data.getSize()), dimensions(dimensions), k(k) {
  }
  std::vector<int> create_graph(void) {
    this->start_slave_code();
    MPI_Status stat;
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
    unsigned int packagesize = 2000;
    if (packagesize > datasize / dimensions / MPIEnviroment::get_node_count()) {
      packagesize = static_cast<unsigned int>(datasize / dimensions /
                                              MPIEnviroment::get_node_count());
    }
    unsigned int send_packageindex = 0;
    unsigned int received_packageindex = 0;
    unsigned int packagecount = static_cast<unsigned int>(datasize / dimensions / packagesize);
    int packageinfo[2];
    unsigned int *startindices = new unsigned int[MPIEnviroment::get_node_count()-1];
    packageinfo[1] = packagesize;
    // Send first packages
    for (int dest = 1; dest < MPIEnviroment::get_node_count(); dest++) {
      packageinfo[0] = send_packageindex * packagesize;
      MPI_Send(packageinfo, 2, MPI_INT, dest, 1, MPI_COMM_WORLD);
      startindices[dest - 1] = packageinfo[0];
      send_packageindex++;
    }

    // Receive and send packages until all are done
    int messagesize;
    while (received_packageindex != packagecount+1) {
      MPI_Probe(MPI_ANY_SOURCE, 1, MPI_COMM_WORLD, &stat);
      MPI_Get_count(&stat, MPI_INT, &messagesize);  // Count should be packagesize*k
      std::cerr << "Received graph package [" << received_packageindex+1
                << " / " << packagecount+1 << "] from node "<< stat.MPI_SOURCE
                << "! Messagesize: " << messagesize << std::endl;
      int source = stat.MPI_SOURCE;
      int start = startindices[source-1];
      int *partial_result = new int[messagesize];
      MPI_Recv(partial_result, messagesize, MPI_INT, stat.MPI_SOURCE,
               stat.MPI_TAG, MPI_COMM_WORLD, &stat);

      // Send next package
      if (send_packageindex < packagecount) {
        packageinfo[0] = send_packageindex * packagesize;
        MPI_Send(packageinfo, 2, MPI_INT, source, 1, MPI_COMM_WORLD);
        startindices[source - 1] = packageinfo[0];
        send_packageindex++;
      } else if (send_packageindex == packagecount) {
        // Send last package
        packageinfo[0] = send_packageindex * packagesize;
        packageinfo[1] = static_cast<unsigned int>((datasize / dimensions) % packagesize);
        if (packageinfo[1] == 0) {
          send_packageindex++;
          std::cerr << "Received graph package [" << received_packageindex+1
                    << " / " << packagecount+1 << "] (empty package)" << std::endl;
          received_packageindex++;
        } else {
          MPI_Send(packageinfo, 2, MPI_INT, source, 1, MPI_COMM_WORLD);
          startindices[source-1] = packageinfo[0];
          send_packageindex++;
          // std::cerr<<"Last Package"<<send_packageindex-1<<" sent! Size: "<<packageinfo[1];
        }
      } else {
        std::cerr << "No more graph packages to send!" << std::endl;
      }
      // Store result
      for (int i = 0; i < messagesize; i++)
        returngraph[start*k+i] = partial_result[i];
      delete [] partial_result;
      received_packageindex++;
    }
    int exitmessage[2] = {-2, -2};
    for (int dest = 1; dest < MPIEnviroment::get_node_count(); dest++)
      MPI_Send(exitmessage, 2, MPI_INT, dest, 1, MPI_COMM_WORLD);
    delete [] startindices;
    delete [] graph;
    return returngraph;
  }
  virtual ~OperationCreateGraphMPI() {}
};

#endif /* OPERATIONCREATEGRAPHMPI_H */
