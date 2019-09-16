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

/// Base class for MPI master node operations
class MPIWorkerBase {
 private:
  static int index;
  int object_index;
  std::string operationName;

 protected:
  bool verbose;

 public:
  explicit MPIWorkerBase(std::string worker_class_name);
  MPIWorkerBase();
  virtual ~MPIWorkerBase(void);
  virtual void start_worker_main(void) = 0;
  void start_sub_workers(void);
  void release_sub_workers(void);
};

class WorkerDummy : public MPIWorkerBase {
 public:
  explicit WorkerDummy(std::string worker_name);
  virtual void start_worker_main(void);
  virtual ~WorkerDummy() {}
};

class OperationDummy : protected WorkerDummy {
 public:
  OperationDummy(void)
      : WorkerDummy("OPDummy") {
  }
  void start_operation(void) {
    std::cout << "Press any key to start dummy operation " << std::endl;
    std::cin.get();
    start_sub_workers();
    std::cout << "Press any key to exit dummy operation " << std::endl;
    std::cin.get();
  }
};





}  // namespace clusteringmpi
}  // namespace datadriven
}  // namespace sgpp
#endif /* OPERATIONMPI_H */
