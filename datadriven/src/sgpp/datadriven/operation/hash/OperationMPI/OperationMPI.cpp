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
NodeCommunicator::NodeCommunicator(int masternode, base::OperationConfiguration &conf) :
    configuration(conf) {
  // Get Slave MPI IDs to construct an MPI_Group
  std::vector<int> slave_ids;
  slave_ids.push_back(masternode);  // Masternode on position 0
  int slaveid = MPIEnviroment::get_node_rank() + 1;
  worker_count = 0;
  for (std::string &slaveName : conf["SLAVES"].keys()) {
    slave_ids.push_back(slaveid);
    slaveid += count_slaves(conf["SLAVES"][slaveName]) + 1;
    worker_count++;
  }
  MPI_Group world_group;
  MPI_Group worker_group;
  MPI_Comm_group(MPI_COMM_WORLD, &world_group);
  MPI_Group_incl(world_group, slave_ids.size(), slave_ids.data(), &worker_group);
  MPI_Comm_create(MPI_COMM_WORLD, worker_group, &communicator);
}

void NodeCommunicator::spawn_workers(std::string &operationName) {
  // Command for creation and execution of a slave
  int message[1];
  message[0] = 1;
  for (int i = 1; i < worker_count; i++) {
    MPI_Send(message, static_cast<int>(1), MPI_INT, i, 1, MPI_COMM_WORLD);
  }

  // Send Configuration to other processes
  size_t slaveid = 1;   // First one is always the master worker
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

  // Send signal to create worker by using the factory
  char *class_message = new char[operationName.size() + 1];
  snprintf(class_message, operationName.size() + 1, "%s", operationName.c_str());
  for (size_t i = 1; i < slaveid; i++) {
    MPI_Send(class_message, static_cast<int>(operationName.size() + 1),
             MPI_CHAR, i, 1, communicator);
  }
}

void NodeCommunicator::start_workers(int object_index) {
  int message[1];
  // Command for creation and execution of a slave
  message[0] = object_index + 10;
  for (int i = 1; i < worker_count; i++) {
    MPI_Send(message, static_cast<int>(1), MPI_INT, i, 1, communicator);
  }
}
void NodeCommunicator::release_workers(int object_index) {
  int message[1];
  // Kill all worker objects
  message[0] = 2;
  for (int i = 1; i < worker_count; i++) {
    MPI_Send(message, static_cast<int>(1), MPI_INT, i, 1, communicator);
  }
  message[0] = object_index + 10;
  for (int i = 1; i < worker_count; i++) {
    MPI_Send(message, static_cast<int>(1), MPI_INT, i, 1, communicator);
  }
}

MPIWorkerBase::MPIWorkerBase(void) {
}
MPIWorkerBase::MPIWorkerBase(int masternode, std::string &operationName,
                             base::OperationConfiguration conf) : configuration(conf) {
  // Create communicator
  node_comm = new NodeCommunicator(masternode, configuration);
  // Create worker processes if necessary
  node_comm->spawn_workers(operationName);
  // Get Configuration
  if (conf.contains("VERBOSE"))
    verbose = configuration["VERBOSE"].getBool();
  if (conf.contains("OPENCL_NODE"))
    opencl_node = configuration["OPENCL_NODE"].getBool();
}

MPIWorkerBase::~MPIWorkerBase(void) {
  delete node_comm;
}
MPIOperation::MPIOperation(base::OperationConfiguration &conf,
                           std::string &operationName) : object_index(index), verbose(false),
                                                           configuration(conf) {
  if (conf.contains("VERBOSE"))
    verbose = conf["VERBOSE"].getBool();

  if (verbose)
    std::cerr << "Created Master Node at process " << MPIEnviroment::get_node_rank() << std::endl;
  index++;

  // Create communicator
  node_comm = new NodeCommunicator(0, configuration);
  // Create worker processes if necessary
  node_comm->spawn_workers(operationName);
}

int NodeCommunicator::count_slaves(json::Node &currentslave) {
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
MPIOperation::~MPIOperation(void) {
  delete node_comm;
}

}  // namespace clusteringmpi
}  // namespace datadriven
}  // namespace sgpp
