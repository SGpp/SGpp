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

namespace sgpp {
namespace datadriven {
namespace clusteringmpi {
/// Base class for alls MPI slave node operations
class MPISlaveOperation {
 public:
  MPISlaveOperation(void);
  explicit MPISlaveOperation(base::OperationConfiguration conf);
  virtual ~MPISlaveOperation();

  virtual void slave_code(void) = 0;
};

/// Base class for MPI master node operations
class MPIOperation {
 private:
  static int index;
  int object_index;
  bool verbose;
 public:
  /// Constructor - creates all slave operations of the given name
  explicit MPIOperation(base::OperationConfiguration &conf, std::string slave_class_name);
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
