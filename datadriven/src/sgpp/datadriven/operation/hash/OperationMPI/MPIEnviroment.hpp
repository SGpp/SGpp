// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org
#ifndef MPIENVIROMENT_H
#define MPIENVIROMENT_H

#include <mpi.h>
#include <sgpp/base/tools/OperationConfiguration.hpp>
#include <iostream>
#include <string>
#include <sstream>
#include <stdexcept>
#include <vector>

namespace sgpp {
namespace datadriven {
namespace clusteringmpi {

class MPIEnviroment {
 private:
  static MPIEnviroment *singleton_instance;
  base::OperationConfiguration configuration;
  int numTasks;
  int rank;
  bool verbose;
  bool initialized;
  int initialized_worker_counter;
  int initial_source;

  MPI_Comm communicator;
  MPI_Comm input_communicator;
  MPI_Group node_neighbors;
  std::vector<int> neighbor_list;
  int worker_count;

  MPIEnviroment(int argc, char *argv[], bool verbose);
  MPIEnviroment(void);
  MPIEnviroment(MPIEnviroment &cpy);

  void slave_mainloop(void);
  int count_slaves(json::Node &currentslave);
  void init_communicator(base::OperationConfiguration conf);
  void init_worker(int workerid, int source);

 public:
  static void init(int argc, char *argv[], bool verbose = false);
  static void connect_nodes(base::OperationConfiguration conf);
  static void release(void);
  static int get_node_rank(void);
  static int get_node_count(void);

  static MPI_Comm& get_communicator(void) {return singleton_instance->communicator;}
  static base::OperationConfiguration& get_configuration(void)
  {return singleton_instance->configuration;}
  static int get_sub_worker_count(void) {return singleton_instance->worker_count;}
  static MPI_Comm& get_input_communicator(void) {return singleton_instance->input_communicator;}
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
  size_t startindex;
  size_t packagesize;
  size_t workitem_count;
  MPI_Datatype mpi_typ;

  bool verbose;
  MPI_Comm &comm;
  int commsize;

 public:
  SimpleQueue(size_t startindex, size_t workitem_count, size_t node_packagesize, MPI_Comm &comm,
              int commsize) :
      startindex(startindex), packagesize(node_packagesize),
      workitem_count(workitem_count), comm(comm),
      commsize(commsize) {
    // Adapt packagesize
    if (packagesize > workitem_count) {
      commsize = 1;
      packagesize = workitem_count;
    } else if (packagesize > workitem_count / (commsize)) {
      packagesize = static_cast<int>(workitem_count / commsize);
    }
    send_packageindex = 0;
    received_packageindex = 0;
    packagecount = static_cast<unsigned int>(workitem_count / packagesize) + 1;
    startindices = new unsigned int[commsize];
    packageinfo[0] = startindex;
    packageinfo[1] = static_cast<int>(packagesize);

    // Returnvector type
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
    for (int dest = 1; dest < commsize + 1 ; dest++) {
      packageinfo[0] = static_cast<int>(startindex + send_packageindex * packagesize);
      MPI_Send(packageinfo, 2, MPI_INT, dest, 1, comm);
      startindices[dest - 1] = packageinfo[0];
      send_packageindex++;
    }
    verbose = true;
  }
  size_t receive_result(int &startid, T *partial_result) {
    MPI_Status stat;
    int messagesize = 0;
    if (received_packageindex < packagecount+1) {
      MPI_Probe(MPI_ANY_SOURCE, 1, comm, &stat);
      MPI_Get_count(&stat, mpi_typ, &messagesize);
      if (verbose) {
        std::cout << "Received work package [" << received_packageindex+1
                  << " / " << packagecount << "] from node "<< stat.MPI_SOURCE
                  << "! Messagesize: " << messagesize << std::endl;
      }
      int source = stat.MPI_SOURCE;
      startid = startindices[source-1];
      MPI_Recv(partial_result, messagesize, mpi_typ, stat.MPI_SOURCE,
               stat.MPI_TAG, comm, &stat);
      received_packageindex++;

      // Send next package
      if (send_packageindex < packagecount - 1) {
        packageinfo[0] = static_cast<int>(startindex + send_packageindex * packagesize);
        MPI_Send(packageinfo, 2, MPI_INT, source, 1, comm);
        startindices[source - 1] = packageinfo[0];
        send_packageindex++;
      } else if (send_packageindex == packagecount - 1) {
        // Send last package
        packageinfo[0] = static_cast<int>(startindex + send_packageindex * packagesize);
        packageinfo[1] = static_cast<int>((workitem_count) % packagesize);
        if (packageinfo[1] == 0) {
          send_packageindex++;
          received_packageindex++;
          if (verbose)
            std::cout << "Received work package [" << received_packageindex
                      << " / " << packagecount << "] (empty package)" << std::endl;
        } else {
          MPI_Send(packageinfo, 2, MPI_INT, source, 1, comm);
          startindices[source-1] = packageinfo[0];
          send_packageindex++;
        }
      } else {
      }
    } else {
      std::cout << "Error - packagecount: " << packagecount << " Received:"
                << received_packageindex << "\n";
      throw "bla";
    }
    return messagesize;
  }
  bool is_finished(void) {
    if (received_packageindex == packagecount)
      return true;
    else
      return false;
  }
  virtual ~SimpleQueue() {
    delete [] startindices;
  }
};

}  // namespace clusteringmpi
}  // namespace datadriven
}  // namespace sgpp
#endif /* MPIENVIROMENT_H */
