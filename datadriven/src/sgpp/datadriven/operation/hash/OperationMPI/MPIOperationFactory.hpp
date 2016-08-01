// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org
#ifndef MPIOPERATIONFACTORY_H
#define MPIOPERATIONFACTORY_H
#include <sgpp/datadriven/operation/hash/OperationMPI/OperationCreateGraphMPI.hpp>
#include <sgpp/datadriven/operation/hash/OperationMPI/OperationRhsMPI.hpp>
#include <sgpp/datadriven/operation/hash/OperationMPI/OperationDensityMPI.hpp>
#include <sgpp/datadriven/operation/hash/OperationMPI/OperationCreatePrunedGraphMPI.hpp>
namespace sgpp {
namespace datadriven {
namespace clusteringmpi {


MPIWorkerBase* create_mpi_operation(int masternode,
                                        base::OperationConfiguration conf, char *classname) {
  std::cout << classname << std::endl;
  /*if (std::strcmp(classname, "OperationCreateGraphSlave")
      == 0)  {
    return OperationCreateGraphMPI::create_slave(conf);
  }
  if (std::strcmp(classname, "OperationRhsSlave") == 0)  {
    return OperationRhsMPI::create_slave(conf);
    }*/
  /*if (std::strcmp(classname, "OperationDensitySlave") == 0)  {
    return OperationDensityMPI::create_slave(masternode, conf);
    }*/
  /*if (std::strcmp(classname, "OperationCreatePrunedGraphSlave") == 0)  {
    return OperationCreatePrunedGraph::create_slave(conf);
    }*/
  if (std::strcmp(classname, "OPDummy") == 0)  {
    return new WorkerDummy(classname);
  }
  if (std::strcmp(classname, "DensityMultiplicationWorker") == 0)  {
    return new DensityWorker();
  }
  if (std::strcmp(classname, "DensityRHSWorker") == 0)  {
    return new DensityRhsWorker();
  }
  return NULL;
}

}  // namespace clusteringmpi
}  // namespace datadriven
}  // namespace sgpp
#endif /* MPIOPERATIONFACTORY_H */
