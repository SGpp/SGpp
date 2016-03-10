// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org
#ifndef MPIOPERATIONFACTORY_H
#define MPIOPERATIONFACTORY_H
#include <sgpp/datadriven/operation/hash/OperationMPI/OperationCreateGraphMPI.hpp>

MPISlaveOperation* create_mpi_operation(char *classname) {
  if (std::strcmp(classname, "OperationCreateGraphSlave"))  {
    return new OperationCreateGraphSlave();
  }
  return NULL;
}

#endif /* MPIOPERATIONFACTORY_H */
