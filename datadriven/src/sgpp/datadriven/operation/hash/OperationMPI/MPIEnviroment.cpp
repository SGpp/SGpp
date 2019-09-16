// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <algorithm>
#include <string>
#include <vector>

#include <sgpp/datadriven/operation/hash/OperationMPI/MPIEnviroment.hpp>
#include <sgpp/datadriven/operation/hash/OperationMPI/MPIOperationFactory.hpp>

namespace sgpp {
namespace datadriven {
namespace clusteringmpi {
MPIEnviroment *MPIEnviroment::singleton_instance = nullptr;

MPIEnviroment::MPIEnviroment(int argc, char *argv[], bool verbose)
    : numTasks(0), rank(0), verbose(verbose), initialized(false), initialized_worker_counter(0) {
  MPI_Init(&argc, &argv);
  // Gets number of tasks/processes that this program is running on
  MPI_Comm_size(MPI_COMM_WORLD, &numTasks);
  // Gets the rank (process/task number) that this program is running on
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
}
MPIEnviroment::MPIEnviroment(void) {}
MPIEnviroment::MPIEnviroment(MPIEnviroment &cpy) {}

void MPIEnviroment::slave_mainloop(void) {
  MPI_Status stat;
  if (verbose) {
    std::cout << "Node " << rank << ": Started slave mainloop" << std::endl;
  }
  std::vector<MPIWorkerBase *> slave_ops;
  do {
    int messagesize = 0;
    int message_source = -1;
    if (!initialized) {
      MPI_Probe(MPI_ANY_SOURCE, 1, MPI_COMM_WORLD, &stat);
    } else /*if (!(rank == 0))*/ {
      MPI_Probe(MPI_ANY_SOURCE, 1, input_communicator, &stat);
    }
    MPI_Get_count(&stat, MPI_INT, &messagesize);
    message_source = stat.MPI_SOURCE;
    int *message = new int[messagesize];
    if (!initialized)
      MPI_Recv(message, messagesize, MPI_INT, stat.MPI_SOURCE, stat.MPI_TAG, MPI_COMM_WORLD, &stat);
    else  // if (!(rank == 0))
      MPI_Recv(message, messagesize, MPI_INT, stat.MPI_SOURCE, stat.MPI_TAG, input_communicator,
               &stat);
    if (messagesize != 1) {
      std::cout << "Warning! Received illegal message at node " << rank
                << "! Messagesize: " << messagesize << " - Going to ignore this..." << std::endl;
      delete[] message;
      continue;
    } else if (message[0] == -1) {
      std::cout << "Node " << rank << " received cleanup signal..." << std::endl;
      for (auto p : slave_ops) {
        if (p != nullptr) delete p;
      }
      // Send the order to terminate all slaves processes!
      for (int i = 1; i < worker_count + 1; i++) {
        std::cout << "Sending cleanup signal to" << neighbor_list[i] << std::endl;
        MPI_Send(message, static_cast<int>(1), MPI_INT, i, 1, communicator);
      }
      delete[] message;
      break;
    } else if (message[0] == 1) {
      if (!initialized) {
        std::cout << "Warning! Received illegal message at node " << rank
                  << "! Node not yet initialized with configuration file" << std::endl;
        delete[] message;
        continue;
      }

      // Get Operationname and create operation here
      MPI_Probe(message_source, 1, input_communicator, &stat);
      MPI_Get_count(&stat, MPI_CHAR, &messagesize);
      char *classname = new char[messagesize];
      MPI_Recv(classname, messagesize, MPI_CHAR, stat.MPI_SOURCE, stat.MPI_TAG, input_communicator,
               &stat);
      slave_ops.push_back(create_mpi_operation(message_source, configuration, classname));
      if (verbose) {
        std::cout << "Node " << rank << ": Created slave operation \"" << classname << "\""
                  << std::endl;
      }
      delete[] classname;
    } else if (message[0] == 2) {
      if (!initialized) {
        std::cout << "Warning! Received illegal message at node " << rank
                  << "! Node not yet initialized with configuration file" << std::endl;
        delete[] message;
        continue;
      }
      MPI_Probe(message_source, 1, communicator, &stat);
      MPI_Recv(message, messagesize, MPI_INT, stat.MPI_SOURCE, stat.MPI_TAG, input_communicator,
               &stat);
      delete slave_ops[message[0] - 10];
      slave_ops[message[0] - 10] = nullptr;
      if (verbose) {
        std::cout << "Node " << rank << ": Deleted slave operation with ID: " << message[0] - 10
                  << "" << std::endl;
      }
    } else if (message[0] == 3) {
      // Receive node list
      MPI_Probe(message_source, 1, MPI_COMM_WORLD, &stat);
      MPI_Get_count(&stat, MPI_INT, &messagesize);
      int *nodelist = new int[messagesize];
      MPI_Recv(nodelist, messagesize, MPI_INT, stat.MPI_SOURCE, stat.MPI_TAG, MPI_COMM_WORLD,
               &stat);
      bool is_input_comm = false;
      for (int i = 0; i < messagesize; ++i) {
        if (rank == nodelist[i]) {
          is_input_comm = true;
          break;
        }
      }

      // Create Comm
      MPI_Group world_group, tmp_group;
      MPI_Comm tmp_comm;
      MPI_Comm_group(MPI_COMM_WORLD, &world_group);
      MPI_Group_incl(world_group, messagesize, nodelist, &tmp_group);
      if (is_input_comm) {
        MPI_Comm_create(MPI_COMM_WORLD, tmp_group, &input_communicator);
        std::cout << "Created input_comm on " << rank << std::endl;
      } else {
        MPI_Comm_create(MPI_COMM_WORLD, tmp_group, &tmp_comm);
        if (tmp_comm != MPI_COMM_NULL) {
          std::cout << "ERROR - created unnecessary comm!" << std::endl;
        }
      }
      delete[] nodelist;
    } else if (message[0] == 4) {
      // Get serialized configuration and deserialize
      MPI_Probe(message_source, 1, MPI_COMM_WORLD, &stat);
      MPI_Get_count(&stat, MPI_CHAR, &messagesize);
      char *serialized_conf = new char[messagesize];
      MPI_Recv(serialized_conf, messagesize, MPI_CHAR, stat.MPI_SOURCE, stat.MPI_TAG,
               MPI_COMM_WORLD, &stat);
      base::OperationConfiguration conf;
      conf.deserialize(serialized_conf);

      // Init communicator
      init_communicator(conf);
      // Init subworkers if necessary
      init_worker(0, message_source);
      initial_source = message_source;
      delete[] serialized_conf;
    } else if (message[0] == 5) {
      initialized_worker_counter++;
      if (initialized_worker_counter == worker_count) {
        std::cout << "Node " << rank << " initialized!" << std::endl;
        if (rank == 0) {
          return;
        } else {
          int message[1];
          message[0] = 5;
          MPI_Send(message, static_cast<int>(1), MPI_INT, initial_source, 1, MPI_COMM_WORLD);
        }
      } else {
        // Init next worker
        init_worker(initialized_worker_counter, initial_source);
      }
    } else if (message[0] == 6) {
      std::cout << "Node " << rank << " switched to local communicator!" << std::endl;
      initialized = true;
    } else if (message[0] >= 10) {
      // run operation here
      if (verbose) {
        std::cout << "Node " << rank << ": Starting slave operation "
                  << " with index " << message[0] << std::endl;
      }
      if (slave_ops[message[0] - 10] == nullptr)
        throw std::logic_error("Trying to run an non existing slave operation!");
      slave_ops[message[0] - 10]->start_sub_workers();
      slave_ops[message[0] - 10]->start_worker_main();
    }
    delete[] message;
  } while (true);
  MPI_Finalize();
  if (verbose) {
    std::cout << "Node " << rank << ": Exiting... " << std::endl;
  }
  std::exit(0);
}

int MPIEnviroment::count_slaves(json::Node &currentslave) {
  int slavecount = 0;
  if (currentslave.contains("SLAVES")) {
    for (std::string &slaveName : currentslave["SLAVES"].keys()) {
      if (currentslave["SLAVES"][slaveName].contains("SLAVES")) {
        slavecount += count_slaves(currentslave["SLAVES"][slaveName]);
      } else {
        slavecount++;
      }
    }
  }
  return slavecount;
}

void MPIEnviroment::init_communicator(base::OperationConfiguration conf) {
  configuration = conf;
  // Get Slave MPI IDs to construct an MPI_Group
  neighbor_list.push_back(MPIEnviroment::get_node_rank());  // Self at 0
  int slaveid = MPIEnviroment::get_node_rank() + 1;
  worker_count = 0;
  if (conf.contains("SLAVES")) {
    for (std::string &slaveName : conf["SLAVES"].keys()) {
      neighbor_list.push_back(slaveid);
      slaveid += count_slaves(conf["SLAVES"][slaveName]) + 1;
      worker_count++;
    }
  }

  // Send MPI_Group to Slaves - they need to create the same comm!
  int message[1];
  message[0] = 3;
  for (int i = 0; i < MPIEnviroment::get_node_count(); ++i) {
    if (i != rank) MPI_Send(message, static_cast<int>(1), MPI_INT, i, 1, MPI_COMM_WORLD);
  }
  for (int i = 0; i < MPIEnviroment::get_node_count(); ++i) {
    if (i != rank)
      MPI_Send(neighbor_list.data(), static_cast<int>(neighbor_list.size()), MPI_INT, i, 1,
               MPI_COMM_WORLD);
  }

  // Create Comm
  MPI_Group world_group;
  MPI_Comm_group(MPI_COMM_WORLD, &world_group);
  MPI_Group_incl(world_group, static_cast<int>(neighbor_list.size()), neighbor_list.data(),
                 &node_neighbors);

  std::cout << "Created comm on " << MPIEnviroment::get_node_rank() << "\n";
  MPI_Comm_create(MPI_COMM_WORLD, node_neighbors, &communicator);
  if (communicator == MPI_COMM_NULL) {
    throw std::logic_error("Comm on node could not be created! Check configuration file");
  }
}

void MPIEnviroment::init_worker(int workerid, int source) {
  if (!configuration.contains("SLAVES")) {
    // Send init signal back
    int message[1];
    message[0] = 5;
    MPI_Send(message, static_cast<int>(1), MPI_INT, source, 1, MPI_COMM_WORLD);
    std::cout << "Node " << rank << " initialized! Sending signal to " << neighbor_list[1]
              << std::endl;
    return;
  }
  // Command for creation and execution of a slave
  int message[1];
  message[0] = 4;
  int offset = 1;
  for (int i = offset; i < worker_count + offset; i++) {
  }
  // Send Configuration to other processes
  int workercounter = offset;
  for (std::string &slaveName : configuration["SLAVES"].keys()) {
    if (workerid + 1 == workercounter) {
      // Send init_comm message
      MPI_Send(message, static_cast<int>(1), MPI_INT, neighbor_list[workerid + 1], 1,
               MPI_COMM_WORLD);
      json::Node &slave = configuration["SLAVES"][slaveName];
      std::ostringstream sstream;
      slave.serialize(sstream, 0);
      std::string serialized_conf = sstream.str();
      char *conf_message = new char[serialized_conf.length() + 1];
      std::copy(serialized_conf.begin(), serialized_conf.end(), conf_message);
      conf_message[serialized_conf.size()] = '\0';
      // Send configuration
      MPI_Send(conf_message, static_cast<int>(serialized_conf.size() + 1), MPI_CHAR,
               neighbor_list[workerid + 1], 1, MPI_COMM_WORLD);

      delete[] conf_message;
    }

    workercounter++;
  }
}

void MPIEnviroment::init(int argc, char *argv[], bool verbose) {
  if (singleton_instance == nullptr) {
    singleton_instance = new MPIEnviroment(argc, argv, verbose);
    if (singleton_instance->rank != 0) singleton_instance->slave_mainloop();
  } else {
    throw std::logic_error("Singleton class \"MPIEnviroment\" already initialized!");
  }
}

void MPIEnviroment::connect_nodes(base::OperationConfiguration conf) {
  if (singleton_instance != nullptr) {
    if (singleton_instance->rank == 0) {
      singleton_instance->init_communicator(conf);
      singleton_instance->init_worker(0, 0);
      singleton_instance->slave_mainloop();
      int message[1];
      message[0] = 6;
      for (int i = 1; i < MPIEnviroment::get_node_count(); ++i) {
        MPI_Send(message, 1, MPI_INT, i, 1, MPI_COMM_WORLD);
      }
      std::cout << "Network initialized and ready" << std::endl;
      singleton_instance->initialized = true;
    } else {
      throw std::logic_error("connect_nodes should only be called on the the MPI process 0!");
    }
  } else {
    throw std::logic_error("Singleton class \"MPIEnviroment\" not yet initialized!");
  }
}
void MPIEnviroment::release(void) {
  if (singleton_instance != nullptr) {
    std::cout << "Beginning cleanup..." << std::endl;
    int message[1];
    // Send the order to terminate all slaves processes!
    message[0] = -1;
    for (int i = 1; i < singleton_instance->worker_count + 1; i++) {
      MPI_Send(message, static_cast<int>(1), MPI_INT, i, 1, singleton_instance->communicator);
    }
    MPI_Finalize();
    delete singleton_instance;
    std::cout << "Cleanup done!" << std::endl;
  } else {
    throw std::logic_error("Singleton class \"MPIEnviroment\" not yet initialized!");
  }
}

int MPIEnviroment::get_node_rank(void) {
  if (singleton_instance != nullptr)
    return singleton_instance->rank;
  else
    throw std::logic_error("Singleton class \"MPIEnviroment\" not yet initialized!");
}

int MPIEnviroment::get_node_count(void) {
  if (singleton_instance != nullptr)
    return singleton_instance->numTasks;
  else
    throw std::logic_error("Singleton class \"MPIEnviroment\" not yet initialized!");
}

MPIEnviroment::~MPIEnviroment(void) {}

base::OperationConfiguration MPIEnviroment::createMPIConfiguration(int packagesize_master,
                                                                   int leutnant_nodes,
                                                                   int packagesize_leutnants,
                                                                   int slave_nodes_per_leutnant) {
  base::OperationConfiguration conf;
  conf.addIDAttr("PREFETCHING", true);
  conf.addIDAttr("VERBOSE", false);
  conf.addIDAttr("PREFERED_PACKAGESIZE", static_cast<int64_t>(packagesize_master));
  std::unique_ptr<json::Node> workers(new json::DictNode);
  for (int i = 0; i < leutnant_nodes; ++i) {
    std::unique_ptr<json::Node> node_worker(new json::DictNode);
    node_worker->addIDAttr("VERBOSE", false);
    node_worker->addIDAttr("PREFETCHING", true);
    node_worker->addIDAttr("PREFERED_PACKAGESIZE", static_cast<int64_t>(packagesize_leutnants));
    std::unique_ptr<json::Node> device_workers(new json::DictNode);
    for (int device = 0; device < slave_nodes_per_leutnant; device++) {
      std::unique_ptr<json::Node> device_worker(new json::DictNode);
      device_worker->addIDAttr("VERBOSE", false);
      device_worker->addIDAttr("PREFETCHING", true);
      device_worker->addIDAttr("OPENCL_PLATFORM", UINT64_C(0));
      device_worker->addIDAttr("OPENCL_DEVICE", UINT64_C(0));
      std::string id = std::string("SLAVE_WORKER") + std::to_string(device);
      device_workers->addAttribute(id, std::move(device_worker));
    }
    node_worker->addAttribute("SLAVES", std::move(device_workers));
    std::string id = std::string("LEUTNANT") + std::to_string(i);
    workers->addAttribute(id, std::move(node_worker));
  }
  conf.addAttribute("SLAVES", std::move(workers));
  return conf;
}

base::OperationConfiguration MPIEnviroment::createMPIConfiguration(int packagesize,
                                                                   int compute_nodes) {
  base::OperationConfiguration conf;
  conf.addIDAttr("PREFETCHING", true);
  conf.addIDAttr("VERBOSE", false);
  conf.addIDAttr("PREFERED_PACKAGESIZE", static_cast<int64_t>(packagesize));
  std::unique_ptr<json::Node> workers(new json::DictNode);
  for (int i = 0; i < compute_nodes; ++i) {
    std::unique_ptr<json::Node> node_worker(new json::DictNode);
    node_worker->addIDAttr("VERBOSE", false);
    node_worker->addIDAttr("PREFETCHING", true);
    node_worker->addIDAttr("OPENCL_PLATFORM", UINT64_C(0));
    node_worker->addIDAttr("OPENCL_DEVICE", UINT64_C(0));
    std::string id = std::string("WORKER") + std::to_string(i);
    workers->addAttribute(id, std::move(node_worker));
  }
  conf.addAttribute("SLAVES", std::move(workers));
  return conf;
}
base::OperationConfiguration MPIEnviroment::createMPIConfiguration(
    int compute_nodes, base::OCLOperationConfiguration node_opencl_conf) {
  base::OperationConfiguration conf;
  conf.addIDAttr("PACKAGE_BUFFERING", false);
  conf.addIDAttr("VERBOSE", false);
  conf.addIDAttr("PREFERED_PACKAGESIZE", UINT64_C(10240));
  std::unique_ptr<json::Node> workers(new json::DictNode);
  for (int i = 0; i < compute_nodes; ++i) {
    std::unique_ptr<json::Node> node_worker(new json::DictNode);
    node_worker->addIDAttr("VERBOSE", false);
    node_worker->addIDAttr("PREFERED_PACKAGESIZE", UINT64_C(2560));
    std::unique_ptr<json::Node> device_workers(new json::DictNode);
    int platform = 0;
    for (std::string &platformName : node_opencl_conf["PLATFORMS"].keys()) {
      int device = 0;
      json::Node &platformNode = node_opencl_conf["PLATFORMS"][platformName];
      for (std::string &deviceName : platformNode["DEVICES"].keys()) {
        std::unique_ptr<json::Node> device_worker(new json::DictNode);
        device_worker->addIDAttr("VERBOSE", false);
        device_worker->addIDAttr("OPENCL_PLATFORM", static_cast<int64_t>(platform));
        device_worker->addIDAttr("OPENCL_DEVICE", static_cast<int64_t>(device));
        std::string id =
            std::string("OPENCL_WORKER_") + platformName + std::string("_") + deviceName;
        device_workers->addAttribute(id, std::move(device_worker));
        device++;
      }
      platform++;
    }
    node_worker->addAttribute("SLAVES", std::move(device_workers));
    std::string id = std::string("WORKER") + std::to_string(i);
    workers->addAttribute(id, std::move(node_worker));
  }
  conf.addAttribute("SLAVES", std::move(workers));
  return conf;
}

}  // namespace clusteringmpi
}  // namespace datadriven
}  // namespace sgpp
