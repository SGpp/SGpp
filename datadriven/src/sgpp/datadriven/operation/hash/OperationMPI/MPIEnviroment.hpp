// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org
#ifndef MPIENVIROMENT_H
#define MPIENVIROMENT_H

#include <mpi.h>
#include <iostream>
#include <stdexcept>

namespace sgpp {
namespace datadriven {
namespace clusteringmpi {

class MPIEnviroment {
 private:
  static MPIEnviroment *singleton_instance;
  int numTasks;
  int rank;
  bool verbose;

  MPIEnviroment(int argc, char *argv[], bool verbose);
  MPIEnviroment(void);
  MPIEnviroment(MPIEnviroment &cpy);
  void slave_mainloop(void);

 public:
  static void init(int argc, char *argv[], bool verbose = false);
  static void release(void);
  static int get_node_rank(void);
  static int get_node_count(void);
  ~MPIEnviroment();
};

}  // namespace clusteringmpi
}  // namespace datadriven
}  // namespace sgpp
#endif /* MPIENVIROMENT_H */
