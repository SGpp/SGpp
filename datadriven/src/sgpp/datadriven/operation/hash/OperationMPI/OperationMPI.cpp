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
#include "OperationMPI.hpp"
#include "MPIEnviroment.hpp"

namespace sgpp {
namespace datadriven {
namespace clusteringmpi {
void ArchitecturCommunicator::spawn_workers(std::string workerOperation) {
  // Send Configuration to slaves
  size_t slaveid = 0;
  for (std::string &slaveName : configuration["SLAVES"].keys()) {
    json::Node& slave = configuration["SLAVES"][slaveName];
    std::ostringstream sstream;
    slave.serialize(sstream, 0);
    std::string serialized_conf = sstream.str();
    char *conf_message = new char[serialized_conf.length() + 1];
    std::copy(serialized_conf.begin(), serialized_conf.end(), conf_message);
    conf_message[serialized_conf.size()] = '\0';
    MPI_Send(conf_message, static_cast<int>(serialized_conf.size() + 1),
             MPI_CHAR, slaveid, 1, communicator);
    delete [] conf_message;
    slaveid++;
  }

  char *class_message = new char[workerOperation.size() + 1];
  snprintf(class_message, workerOperation.size() + 1, "%s", workerOperation.c_str());
  for (size_t i = 0; i < slaveid; i++) {
    MPI_Send(class_message, static_cast<int>(workerOperation.size() + 1),
             MPI_CHAR, i, 1, communicator);
  }
}


MPISlaveOperation::MPISlaveOperation(void) {
}
MPISlaveOperation::MPISlaveOperation(base::OperationConfiguration conf) {
}

MPISlaveOperation::~MPISlaveOperation(void) {
}
MPIOperation::MPIOperation(base::OperationConfiguration &conf,
                           std::string slave_class_name) : object_index(index), verbose(false),
                                                           conf(conf) {
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

  // Get Slave MPI IDs to construct an MPI_Group
  std::vector<int> slave_ids;
  int slaveid = MPIEnviroment::get_node_rank() + 1;
  for (std::string &slaveName : conf["SLAVES"].keys()) {
    slave_ids.push_back(slaveid);
    slaveid += count_slaves(conf["SLAVES"][slaveName]) + 1;
  }
  MPI_Group world_group;
  MPI_Comm_group(MPI_COMM_WORLD, &world_group);
  MPI_Group_incl(world_group, slave_ids.size(), slave_ids.data(), &slave_group);
  MPI_Comm_create(MPI_COMM_WORLD, slave_group, &slave_comm);

  // Send Configuration to slaves
  slaveid = 0;
  for (std::string &slaveName : conf["SLAVES"].keys()) {
    json::Node& slave = conf["SLAVES"][slaveName];
    std::ostringstream sstream;
    slave.serialize(sstream, 0);
    std::string serialized_conf = sstream.str();
    char *conf_message = new char[serialized_conf.length() + 1];
    std::copy(serialized_conf.begin(), serialized_conf.end(), conf_message);
    conf_message[serialized_conf.size()] = '\0';
    MPI_Send(conf_message, static_cast<int>(serialized_conf.size() + 1),
             MPI_CHAR, slaveid, 1, slave_comm);
    delete [] conf_message;
    slaveid++;
  }
  slave_count = slaveid;

  char *class_message = new char[slave_class_name.size() + 1];
  snprintf(class_message, slave_class_name.size() + 1, "%s", slave_class_name.c_str());
  for (size_t i = 0; i < slave_count; i++) {
    MPI_Send(class_message, static_cast<int>(slave_class_name.size() + 1),
             MPI_CHAR, i, 1, slave_comm);
  }
  delete [] class_message;
}

int MPIOperation::count_slaves(json::Node &currentslave) {
  int slavecount = 0;
  for (std::string &slaveName : currentslave["SLAVES"].keys()) {
    if (currentslave["SLAVES"][slaveName].contains("SLAVES")) {
      slavecount += count_slaves(currentslave["SLAVES"][slaveName]);
    } else {
      slavecount++;
    }
  }
  return slavecount;
}

MPIOperation::MPIOperation(void) {
}

int MPIOperation::index = 0;
void MPIOperation::start_slave_code(void) {
  int message[1];
  // Command for creation and execution of a slave
  message[0] = object_index + 10;
  for (size_t i = 0; i < slave_count; i++) {
    MPI_Send(message, static_cast<int>(1), MPI_INT, i, 1, MPI_COMM_WORLD);
  }
}
void MPIOperation::release_slave_objects(void) {
  int message[1];
  // Kill all slave objects
  message[0] = 2;
  for (size_t i = 0; i < slave_count; i++) {
    MPI_Send(message, static_cast<int>(1), MPI_INT, i, 1, MPI_COMM_WORLD);
  }
  message[0] = object_index + 10;
  for (size_t i = 0; i < slave_count; i++) {
    MPI_Send(message, static_cast<int>(1), MPI_INT, i, 1, MPI_COMM_WORLD);
  }
}
MPIOperation::~MPIOperation(void) {
}

}  // namespace clusteringmpi
}  // namespace datadriven
}  // namespace sgpp
