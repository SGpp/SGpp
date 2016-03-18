// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org
#ifndef OPERATIONGRAPHBASEMPI_H
#define OPERATIONGRAPHBASEMPI_H

#include <mpi.h>
#include <sgpp/globaldef.hpp>
#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/datadriven/operation/hash/OperationMPI/MPIEnviroment.hpp>
#include <sgpp/datadriven/operation/hash/OperationMPI/OperationMPI.hpp>

namespace sgpp {
namespace datadriven {
namespace clusteringmpi {

/// Base class for kNN graph operations with MPI
class OperationGraphMethodMPI : public MPIOperation {
 public:
  virtual ~OperationGraphMethodMPI(void) {
  }

 protected:
  sgpp::base::DataMatrix &data;
  size_t datasize;
  size_t dimensions;
  size_t k;

  OperationGraphMethodMPI(sgpp::base::DataMatrix& data, size_t dimensions, size_t k,
                          std::string slave_class_id)
      : MPIOperation(slave_class_id), data(data),
        datasize(data.getSize()), dimensions(dimensions), k(k) {
    // Sending dataset to slaves
    for (int dest = 1; dest < MPIEnviroment::get_node_count(); dest++)
      MPI_Send(data.getPointer(), static_cast<int>(data.getSize()), MPI_DOUBLE,
               dest, 1, MPI_COMM_WORLD);
    for (int dest = 1; dest < MPIEnviroment::get_node_count(); dest++)
      MPI_Send(&dimensions, 1, MPI_INT, dest, 1, MPI_COMM_WORLD);
    for (int dest = 1; dest < MPIEnviroment::get_node_count(); dest++)
      MPI_Send(&k, 1, MPI_INT, dest, 1, MPI_COMM_WORLD);
  }

 protected:
  /// Slave base class for kNN graph operations
  class OperationGraphMethodSlave : public MPISlaveOperation {
   protected:
    bool verbose;
    int k;
    int dimensions;
    double *dataset;
    int dataset_size;

   public:
    OperationGraphMethodSlave(void) : verbose(true) {
      MPI_Status stat;
      // Receive dataset
      MPI_Probe(0, 1, MPI_COMM_WORLD, &stat);
      MPI_Get_count(&stat, MPI_DOUBLE, &dataset_size);
      dataset = new double[dataset_size];
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
    }
    virtual ~OperationGraphMethodSlave(void) {
      delete [] dataset;
    }
    virtual void slave_code(void) = 0;
  };
};
}  // namespace clusteringmpi
}  // namespace datadriven
}  // namespace sgpp
#endif /* OPERATIONGRAPHBASEMPI_H */
