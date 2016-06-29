// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org
#ifndef MPIENVIROMENT_H
#define MPIENVIROMENT_H

#include <mpi.h>
#include <iostream>
#include <string>
#include <sstream>
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

template <class T>
class SimpleQueue {
 protected:
  unsigned int send_packageindex;
  unsigned int received_packageindex;
  unsigned int packagecount;
  int packageinfo[2];
  unsigned int *startindices;
  size_t packagesize;
  size_t workitem_count;
  MPI_Datatype mpi_typ;

  bool verbose;

 public:
  SimpleQueue(size_t workitem_count, size_t packagesize)  : packagesize(packagesize),
                                                            workitem_count(workitem_count) {
    if (packagesize > workitem_count / (MPIEnviroment::get_node_count() - 1)) {
      packagesize = static_cast<int>(workitem_count /
                                              MPIEnviroment::get_node_count());
    }
    send_packageindex = 0;
    received_packageindex = 0;
    packagecount = static_cast<unsigned int>(workitem_count / packagesize);
    startindices = new unsigned int[MPIEnviroment::get_node_count()-1];
    packageinfo[0] = 0;
    packageinfo[1] = static_cast<int>(packagesize);

    if (std::is_same<T, int>::value) {
      mpi_typ = MPI_INT;
    } else if (std::is_same<T, float>::value) {
      mpi_typ = MPI_FLOAT;
    } else if (std::is_same<T, double>::value) {
      mpi_typ = MPI_DOUBLE;
    } else {
      std::stringstream errorString;
      errorString << "Unsupported datatyp in class SimpleQueue." << std::endl
                  << "Template class needs to be int, float or double." << std::endl;
      throw std::logic_error(errorString.str());
    }
    // Send first packages
    std::cout << "Sending size: " << packageinfo[1] << std::endl;
    for (int dest = 1; dest < MPIEnviroment::get_node_count(); dest++) {
      packageinfo[0] = static_cast<int>(send_packageindex * packagesize);
      MPI_Send(packageinfo, 2, MPI_INT, dest, 1, MPI_COMM_WORLD);
      startindices[dest - 1] = packageinfo[0];
      send_packageindex++;
    }
    verbose = true;
  }
  size_t receive_result(int &startid, T *partial_result) {
    MPI_Status stat;
    int messagesize = 0;
    if (received_packageindex != packagecount+1) {
      MPI_Probe(MPI_ANY_SOURCE, 1, MPI_COMM_WORLD, &stat);
      MPI_Get_count(&stat, mpi_typ, &messagesize);  // Count should be packagesize*k
      if (verbose) {
        std::cerr << "Received work package [" << received_packageindex+1
                  << " / " << packagecount+1 << "] from node "<< stat.MPI_SOURCE
                  << "! Messagesize: " << messagesize << std::endl;
      }
      int source = stat.MPI_SOURCE;
      startid = startindices[source-1];
      MPI_Recv(partial_result, messagesize, mpi_typ, stat.MPI_SOURCE,
               stat.MPI_TAG, MPI_COMM_WORLD, &stat);

      // Send next package
      if (send_packageindex < packagecount) {
        packageinfo[0] = static_cast<int>(send_packageindex * packagesize);
        MPI_Send(packageinfo, 2, MPI_INT, source, 1, MPI_COMM_WORLD);
        startindices[source - 1] = packageinfo[0];
        send_packageindex++;
      } else if (send_packageindex == packagecount) {
        // Send last package
        packageinfo[0] = static_cast<int>(send_packageindex * packagesize);
        packageinfo[1] = static_cast<int>((workitem_count) % packagesize);
        if (packageinfo[1] == 0) {
          send_packageindex++;
          received_packageindex++;
          if (verbose)
            std::cerr << "Received work package [" << received_packageindex+1
                      << " / " << packagecount+1 << "] (empty package)" << std::endl;
          return 0;
        } else {
          MPI_Send(packageinfo, 2, MPI_INT, source, 1, MPI_COMM_WORLD);
          startindices[source-1] = packageinfo[0];
          send_packageindex++;
        }
      } else {
      }
    }
    received_packageindex++;
    return messagesize;
  }
  virtual ~SimpleQueue() {
    delete [] startindices;
  }
};


}  // namespace clusteringmpi
}  // namespace datadriven
}  // namespace sgpp
#endif /* MPIENVIROMENT_H */
