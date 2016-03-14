// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org
#ifndef OPERATIONMPI_H
#define OPERATIONMPI_H

#include <mpi.h>
#include <cstring>
#include <string>
#include <typeinfo>

namespace sgpp {
namespace datadriven {
namespace clusteringmpi {
class MPISlaveOperation {
 public:
  MPISlaveOperation();
  virtual ~MPISlaveOperation();

  virtual void slave_code(void) = 0;
};

class MPIOperation {
 private:
  static int index;
  int object_index;
 public:
  explicit MPIOperation(std::string slave_class_name);
  virtual ~MPIOperation(void);
  void start_slave_code(void);
  void release_slave_objects(void);
};

}  // namespace clusteringmpi
}  // namespace datadriven
}  // namespace sgpp
#endif /* OPERATIONMPI_H */
