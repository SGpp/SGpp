// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org
#include <mpi.h>
#include <cstring>
#include <string>
#include <typeinfo>
#include "OperationMPI.hpp"
#include "MPIEnviroment.hpp"

MPISlaveOperation::MPISlaveOperation(void) {
}

MPISlaveOperation::~MPISlaveOperation(void) {
}
MPIOperation::MPIOperation(std::string slave_class_name) {
  int message[1];
  // Command for creation and execution of a slave
  message[0] = 1;
  for (int i = 1; i < MPIEnviroment::get_node_count(); i++) {
    MPI_Send(message, static_cast<int>(1), MPI_INT, i, 1, MPI_COMM_WORLD);
  }

  char *class_message = new char[slave_class_name.size()];
  snprintf(class_message, slave_class_name.size(), "%s", slave_class_name.c_str());
  for (int i = 1; i < MPIEnviroment::get_node_count(); i++) {
    MPI_Send(class_message, static_cast<int>(slave_class_name.size()),
             MPI_CHAR, i, 1, MPI_COMM_WORLD);
  }
  delete [] class_message;
}

void MPIOperation::start_slave_code(void) {
  int message[1];
  // Command for creation and execution of a slave
  message[0] = 2;
  for (int i = 1; i < MPIEnviroment::get_node_count(); i++) {
    MPI_Send(message, static_cast<int>(1), MPI_INT, i, 1, MPI_COMM_WORLD);
  }
}
MPIOperation::~MPIOperation(void) {
}
