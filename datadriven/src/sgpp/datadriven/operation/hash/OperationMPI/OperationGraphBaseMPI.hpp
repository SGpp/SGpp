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
#include <string>

namespace sgpp {
namespace datadriven {
namespace clusteringmpi {
class MPIWorkerGraphBase : virtual public MPIWorkerBase {
 protected:
  double *dataset;
  int dataset_size;
  int k;
  int dimensions;

  MPIWorkerGraphBase(std::string operationName, sgpp::base::DataMatrix &data,
                     int k) : MPIWorkerBase(operationName), dataset(data.getPointer()),
                                 dataset_size(static_cast<int>(data.getSize())), k(k),
                                 dimensions(static_cast<int>(data.getNcols())),
                                 delete_dataset(true) {
      std::cout << "Node " << MPIEnviroment::get_node_rank() << ": Received dataset (size "
                << dataset_size / dimensions << " datapoints) and graph parameters"
                << std::endl;
    send_dataset();
    delete_dataset = false;
  }
  MPIWorkerGraphBase(sgpp::base::DataMatrix &data,
                     int k) : MPIWorkerBase(), dataset(data.getPointer()),
                                 dataset_size(static_cast<int>(data.getSize())), k(k),
                                 dimensions(static_cast<int>(data.getNcols())),
                                 delete_dataset(true) {
    send_dataset();
    delete_dataset = false;
  }
  explicit MPIWorkerGraphBase(std::string operationName) : MPIWorkerBase(operationName) {
    receive_dataset();
    send_dataset();
  }

  virtual ~MPIWorkerGraphBase(void) {
    if (delete_dataset)
      delete [] dataset;
  }

 private:
  bool delete_dataset;
  void send_dataset() {
    // Sending dataset to slaves
    for (int i = 1; i < MPIEnviroment::get_sub_worker_count() + 1; i++) {
      MPI_Send(dataset, dataset_size, MPI_DOUBLE,
               i, 1, MPIEnviroment::get_communicator());
    }
    for (int i = 1; i < MPIEnviroment::get_sub_worker_count() + 1; i++) {
      MPI_Send(&dimensions, 1, MPI_INT, i, 1, MPIEnviroment::get_communicator());
    }
    for (int i = 1; i < MPIEnviroment::get_sub_worker_count() + 1; i++) {
      MPI_Send(&k, 1, MPI_INT, i, 1, MPIEnviroment::get_communicator());
    }
  }
  void receive_dataset(void) {
    MPI_Status stat;
    // Receive dataset
    MPI_Probe(0, 1, MPIEnviroment::get_input_communicator(), &stat);
    MPI_Get_count(&stat, MPI_DOUBLE, &dataset_size);
    dataset = new double[dataset_size];
    MPI_Recv(dataset, dataset_size, MPI_DOUBLE, stat.MPI_SOURCE, stat.MPI_TAG,
             MPIEnviroment::get_input_communicator(), &stat);
    // Receive clustering parameters
    MPI_Probe(0, 1, MPIEnviroment::get_input_communicator(), &stat);
    MPI_Recv(&dimensions, 1, MPI_INT, stat.MPI_SOURCE, stat.MPI_TAG,
             MPIEnviroment::get_input_communicator(), &stat);
    MPI_Probe(0, 1, MPIEnviroment::get_input_communicator(), &stat);
    MPI_Recv(&k, 1, MPI_INT, stat.MPI_SOURCE, stat.MPI_TAG,
             MPIEnviroment::get_input_communicator(), &stat);
    if (verbose) {
      std::cout << "Node " << MPIEnviroment::get_node_rank() << ": Received dataset (size "
                << dataset_size / dimensions << " datapoints) and graph parameters"
                << std::endl;
    }
  }
};

}  // namespace clusteringmpi
}  // namespace datadriven
}  // namespace sgpp
#endif /* OPERATIONGRAPHBASEMPI_H */
