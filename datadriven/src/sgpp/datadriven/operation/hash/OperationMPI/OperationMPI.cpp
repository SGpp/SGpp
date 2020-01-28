// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <mpi.h>
#include <cstring>
#include <string>
#include <typeinfo>
#include <algorithm>
#include <vector>
#include <sgpp/datadriven/operation/hash/OperationMPI/OperationMPI.hpp>
#include <sgpp/datadriven/operation/hash/OperationMPI/MPIEnviroment.hpp>

namespace sgpp {
namespace datadriven {
namespace clusteringmpi {

MPIWorkerBase::MPIWorkerBase(std::string operationName) : object_index(index),
                                                          operationName(operationName),
                                                          verbose(false) {
  std::cout << "In Base cstr" << std::endl;
  std::cout << "Creating operation on " << MPIEnviroment::get_node_rank() << std::endl;
  if (MPIEnviroment::get_configuration().contains("VERBOSE"))
    verbose = MPIEnviroment::get_configuration()["VERBOSE"].getBool();

  int message[1];
  index++;
  // Command for creation and execution of a worker
  message[0] = 1;
  for (int i = 1; i < MPIEnviroment::get_sub_worker_count() + 1; i++) {
    MPI_Send(message, static_cast<int>(1), MPI_INT, i, 1, MPIEnviroment::get_communicator());
  }

  char *class_message = new char[operationName.size() + 1];
  snprintf(class_message, operationName.size() + 1, "%s", operationName.c_str());
  for (int i = 1; i < MPIEnviroment::get_sub_worker_count() + 1; i++) {
    MPI_Send(class_message, static_cast<int>(operationName.size() + 1),
             MPI_CHAR, i, 1, MPIEnviroment::get_communicator());
  }
  delete [] class_message;
  if (verbose) {
    std::cout << "Created Worker on node " << MPIEnviroment::get_node_rank() << " with operation "
              << operationName << " and workerindex" << object_index << std::endl;
  }
}

MPIWorkerBase::MPIWorkerBase(void) {
  std::cout << "In default Base cstr" << std::endl;
}

int MPIWorkerBase::index = 0;
MPIWorkerBase::~MPIWorkerBase(void) {
}

void MPIWorkerBase::start_sub_workers() {
  int message[1];
  // Command for creation and execution of a worker
  message[0] = object_index + 10;
  for (int i = 1; i < MPIEnviroment::get_sub_worker_count() + 1; i++) {
    MPI_Send(message, static_cast<int>(1), MPI_INT, i, 1, MPIEnviroment::get_communicator());
  }
}

void MPIWorkerBase::release_sub_workers(void) {
  int message[1];
  // Kill all slave objects
  message[0] = 2;
  for (int i = 1; i < MPIEnviroment::get_sub_worker_count() + 1; i++) {
    MPI_Send(message, static_cast<int>(1), MPI_INT, i, 1, MPIEnviroment::get_communicator());
  }
  message[0] = object_index + 10;
  for (int i = 1; i < MPIEnviroment::get_sub_worker_count() + 1; i++) {
    MPI_Send(message, static_cast<int>(1), MPI_INT, i, 1, MPIEnviroment::get_communicator());
  }
}


WorkerDummy::WorkerDummy(std::string operationName)
    : MPIWorkerBase(operationName) {
  if (verbose) {
    std::cout << "Created Worker dummy on node " << MPIEnviroment::get_node_rank() << std::endl;
  }
}
void WorkerDummy::start_worker_main() {
  std::cout << "Dummy Main!!" << std::endl;
}

}  // namespace clusteringmpi
}  // namespace datadriven
}  // namespace sgpp
