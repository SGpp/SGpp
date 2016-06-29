// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org
#include <mpi.h>
#include <cstring>
#include <string>
#include <typeinfo>
#include <algorithm>
#include "OperationMPI.hpp"
#include "MPIEnviroment.hpp"

namespace sgpp {
namespace datadriven {
namespace clusteringmpi {
MPISlaveOperation::MPISlaveOperation(void) {
}
MPISlaveOperation::MPISlaveOperation(base::OperationConfiguration conf) {
}

MPISlaveOperation::~MPISlaveOperation(void) {
}
MPIOperation::MPIOperation(base::OperationConfiguration &conf,
                           std::string slave_class_name) : object_index(index), verbose(false) {
  if (conf.contains("VERBOSE"))
    verbose = conf["VERBOSE"].getBool();

  if (verbose)
    std::cerr << "Created Master Node at process " << MPIEnviroment::get_node_rank() << std::endl;
  int message[1];
  index++;
  // Command for creation and execution of a slave
  message[0] = 1;
  for (int i = 1; i < MPIEnviroment::get_node_count(); i++) {
    MPI_Send(message, static_cast<int>(1), MPI_INT, i, 1, MPI_COMM_WORLD);
  }

  json::Node& slaves = conf["SLAVES"]["SLAVE1"];
  std::ostringstream sstream;
  slaves.serialize(sstream, 0);
  std::string serialized_conf = sstream.str();
  char *conf_message = new char[serialized_conf.length() + 1];
  std::copy(serialized_conf.begin(), serialized_conf.end(), conf_message);
  conf_message[serialized_conf.size()] = '\0';
  for (int i = 1; i < MPIEnviroment::get_node_count(); i++) {
    MPI_Send(conf_message, static_cast<int>(serialized_conf.size() + 1),
             MPI_CHAR, i, 1, MPI_COMM_WORLD);
  }
  delete [] conf_message;

  char *class_message = new char[slave_class_name.size() + 1];
  snprintf(class_message, slave_class_name.size() + 1, "%s", slave_class_name.c_str());
  for (int i = 1; i < MPIEnviroment::get_node_count(); i++) {
    MPI_Send(class_message, static_cast<int>(slave_class_name.size() + 1),
             MPI_CHAR, i, 1, MPI_COMM_WORLD);
  }
  delete [] class_message;
}
MPIOperation::MPIOperation(void) {
}

int MPIOperation::index = 0;
void MPIOperation::start_slave_code(void) {
  int message[1];
  // Command for creation and execution of a slave
  message[0] = object_index + 10;
  for (int i = 1; i < MPIEnviroment::get_node_count(); i++) {
    MPI_Send(message, static_cast<int>(1), MPI_INT, i, 1, MPI_COMM_WORLD);
  }
}
void MPIOperation::release_slave_objects(void) {
  int message[1];
  // Kill all slave objects
  message[0] = 2;
  for (int i = 1; i < MPIEnviroment::get_node_count(); i++) {
    MPI_Send(message, static_cast<int>(1), MPI_INT, i, 1, MPI_COMM_WORLD);
  }
  message[0] = object_index + 10;
  for (int i = 1; i < MPIEnviroment::get_node_count(); i++) {
    MPI_Send(message, static_cast<int>(1), MPI_INT, i, 1, MPI_COMM_WORLD);
  }
}
MPIOperation::~MPIOperation(void) {
}

}  // namespace clusteringmpi
}  // namespace datadriven
}  // namespace sgpp
