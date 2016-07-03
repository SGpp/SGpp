// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org
#ifndef OPERATIONMPI_H
#define OPERATIONMPI_H
#include <mpi.h>


#include <sgpp/base/tools/OperationConfiguration.hpp>

#include <cstring>
#include <string>
#include <typeinfo>
#include <vector>

namespace sgpp {
namespace datadriven {
namespace clusteringmpi {
class NodeCommunicator {
 protected:
  MPI_Comm communicator;
  int worker_count;
  base::OperationConfiguration &configuration;

  int count_slaves(json::Node &currentslave);
 public:
  void spawn_workers(std::string &workerOperation);
  void start_workers(int object_index);
  void release_workers(int object_index);

  MPI_Comm& get_communicator(void) {return communicator;}
  int get_workercount(void) {return worker_count;}
  NodeCommunicator(int masternode, base::OperationConfiguration &configuration);
};

/// Base class for all MPI slave node operations
class MPIWorkerBase {
 protected:
  base::OperationConfiguration configuration;
  NodeCommunicator *node_comm;

  bool verbose;
  bool opencl_node;

 public:
  MPIWorkerBase(void);
  explicit MPIWorkerBase(int masternode, std::string &operationName,
                         base::OperationConfiguration conf);
  virtual ~MPIWorkerBase();

  virtual void slave_code(void) = 0;
};

/// Base class for MPI master node operations
class MPIOperation {
 private:
  static int index;
  int object_index;
  bool verbose;

 protected:
  base::OperationConfiguration configuration;
  NodeCommunicator *node_comm;

 public:
  /// Constructor - creates all slave operations of the given name
  MPIOperation(base::OperationConfiguration &conf, std::string &slave_class_name);
  /// Constructor - does not create any slave operations
  MPIOperation();
  virtual ~MPIOperation(void);
  void start_slave_code(void);
  void release_slave_objects(void);
};

}  // namespace clusteringmpi
}  // namespace datadriven
}  // namespace sgpp
#endif /* OPERATIONMPI_H */
