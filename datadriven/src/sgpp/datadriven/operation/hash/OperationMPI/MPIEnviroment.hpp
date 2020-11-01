// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef MPIENVIROMENT_H
#define MPIENVIROMENT_H

#include <mpi.h>
#include <sgpp/base/tools/OperationConfiguration.hpp>
#include <sgpp/base/opencl/OCLOperationConfiguration.hpp>
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
  /**
   * Slave main - every MPI node except the master should run in one of those
   *
   * In this main loop, every slave receives the messages to create or run operations.
   */
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
  static base::OperationConfiguration createMPIConfiguration(int packagesize_master,
                                                             int leutnant_nodes,
                                                             int packagesize_leutnants,
                                                             int slave_nodes_per_leutnant);
  static base::OperationConfiguration createMPIConfiguration(int compute_nodes,
                                                             int opencl_devices_per_compute_node);
  /**
   * Create MPI configuration file based on a number of compute nodes their OpenCL configuration file
   *
   * Creates a MPI configuration with 3 Levels. One master node to control all other nodes, and
   * one leutenant node on every compute node which controls all MPI nodes for this compute
   * nodes. For each OpenCL device specified in the opencl configuration it will also create
   * one worker. Duplicate devices will not be recognized hoewever.
   *
   * @param compute_nodes Number of compute nodes
   * @param node_opencl_configuration OpenCL configuration file
   * @return Finished MPI configuration
   */
  static base::OperationConfiguration createMPIConfiguration(int compute_nodes,
                                                             base::OCLOperationConfiguration
                                                             node_opencl_configuration);
  ~MPIEnviroment();
};

template <class T>
class SimpleQueue {
 protected:
  /// Number of sent workapckages
  unsigned int send_packageindex;
  /// Number of received workapckages
  unsigned int received_packageindex;
  /// Number of all workpackages
  unsigned int packagecount;
  /// Array with for sending indices and packagesizes
  int packageinfo[2];
  /// Stores the current starting indices for each work package
  unsigned int *startindices;
  /// Stores the current starting indices for each prefetched work package
  unsigned int *secondary_indices;
  /// Start index for the entire problem chunk of this queue
  size_t startindex;
  size_t packagesize;
  size_t workitem_count;
  /// Datatyp for MPI_SEND and MPI_RECV, only supports float, double and int right now
  MPI_Datatype mpi_typ;

  bool verbose;
  MPI_Comm &comm;
  int commsize;
  /// Prefetching activated?
  bool prefetching;

 public:
  /**
   * Constructor for a simple queue
   *
   * This constructor check whether there are enough packages for all subworkers (workers which
   * receive workpackages from this queue). If there are not enough packages, the packagesize
   * will be modified such that every worker receives at least 1 package (2 if prefetching is
   * activated). This constructor will also already send the first workpackages with a blocking
   * send so use with care!
   *
   * @param startindex Index where the chunk of the problem begins
   * @param workitem_count Number of subworkers
   * @param node_packagesize Prefered packagesize - will be modified if to small
   * @param comm Communicator to the subworkers
   * @param commsize Size of the communictor
   * @param verbose Verbosity on or off
   * @param prefetching Activates package prefetching if this parameter is true
   */
  SimpleQueue(size_t startindex, size_t workitem_count, size_t node_packagesize, MPI_Comm &comm,
              int commsize, bool verbose = false, bool prefetching = false) :
      startindex(startindex), packagesize(node_packagesize),
      workitem_count(workitem_count), verbose(verbose), comm(comm),
      commsize(commsize), prefetching(prefetching) {
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
    send_packageindex = 0;
    received_packageindex = 0;

    // Adapt packagesize
    if (packagesize > workitem_count) {
      commsize = 1;
      packagesize = workitem_count;
    } else if (packagesize > workitem_count / (commsize * 2)) {
      packagesize = static_cast<int>(workitem_count / (commsize * 2));
    }
    if (packagesize % 128 != 0)
       packagesize -= packagesize % 128;

    packagecount = static_cast<unsigned int>(workitem_count / packagesize) + 1;
    startindices = new unsigned int[commsize];
    secondary_indices = new unsigned int[commsize];
    packageinfo[0] = static_cast<int>(startindex);
    packageinfo[1] = static_cast<int>(packagesize);

    // Send first packages
    if (verbose)
      std::cout << "Sending size: " << packageinfo[1] << std::endl;
    for (int dest = 1; dest < commsize + 1 ; dest++) {
      packageinfo[0] = static_cast<int>(startindex + send_packageindex * packagesize);
      MPI_Send(packageinfo, 2, MPI_INT, dest, 1, comm);
      startindices[dest - 1] = packageinfo[0];
      send_packageindex++;
    }
    if (prefetching) {
      if (packagesize == workitem_count) {
        packageinfo[0] = -1;
        packageinfo[1] = -1;
        MPI_Send(packageinfo, 2, MPI_INT, 1, 1, comm);
      } else {
        // Send secondary packages
        for (int dest = 1; dest < commsize + 1 ; dest++) {
          packageinfo[0] = static_cast<int>(startindex + send_packageindex * packagesize);
          MPI_Send(packageinfo, 2, MPI_INT, dest, 1, comm);
          secondary_indices[dest - 1] = packageinfo[0];
          send_packageindex++;
        }
      }
    }
  }
  /**
   * Wait for a result from any of the other workers, send a new workpackage and return said result
   *
   * @param startid Returns the starting index of the package result
   * @param partial_result Buffer containing the package result
   * @return Size of the result array
   */
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
        if (prefetching) {
          startindices[source - 1] = secondary_indices[source - 1];
          secondary_indices[source - 1] = packageinfo[0];
        } else {
          startindices[source - 1] = packageinfo[0];
        }
        send_packageindex++;
      } else if (send_packageindex == packagecount - 1) {
        // Send last package
        packageinfo[0] = static_cast<int>(startindex + send_packageindex * packagesize);
        packageinfo[1] = static_cast<int>((workitem_count) % packagesize);
        if (packageinfo[1] == 0) {
          if (prefetching) {
            startindices[source - 1] = secondary_indices[source - 1];
            packageinfo[0] = -1;
            packageinfo[1] = -1;
            MPI_Send(packageinfo, 2, MPI_INT, source, 1, comm);
          }
          send_packageindex++;
          received_packageindex++;
          if (verbose)
            std::cout << "Received work package [" << received_packageindex
                      << " / " << packagecount << "] (empty package)" << std::endl;
        } else {
          MPI_Send(packageinfo, 2, MPI_INT, source, 1, comm);
          if (prefetching) {
            startindices[source - 1] = secondary_indices[source - 1];
            secondary_indices[source - 1] = packageinfo[0];
          } else {
            startindices[source - 1] = packageinfo[0];
          }
          send_packageindex++;
        }
      } else {
        if (prefetching) {
          startindices[source - 1] = secondary_indices[source - 1];
          // Send empty package to avoid deadlock on clients
          packageinfo[0] = -1;
          packageinfo[1] = -1;
          MPI_Send(packageinfo, 2, MPI_INT, source, 1, comm);
        }
      }
    } else {
      std::cout << "Error - packagecount: " << packagecount << " Received:"
                << received_packageindex << "\n";
      throw std::logic_error("Queue error! Received too many packages!");
    }
    return messagesize;
  }
  /**
   * Checks whether the queue has already finished
   *
   * @return Returns true if the queue has received results for all packages
   */
  bool is_finished() {
    if (received_packageindex == packagecount)
      return true;
    else
      return false;
  }
  virtual ~SimpleQueue() {
    delete [] startindices;
    delete [] secondary_indices;
    std::cerr << "Queue deleted" << std::endl;
  }
};

}  // namespace clusteringmpi
}  // namespace datadriven
}  // namespace sgpp
#endif /* MPIENVIROMENT_H */
