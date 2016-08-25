
// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org
#ifndef OPERATIONPACKAGEBASEMPI_H
#define OPERATIONPACKAGEBASEMPI_H


#include <mpi.h>

#include <sgpp/datadriven/operation/hash/OperationMPI/MPIEnviroment.hpp>
#include <sgpp/datadriven/operation/hash/OperationMPI/OperationMPI.hpp>
#include <sgpp/datadriven/operation/hash/OperationDensityOCLMultiPlatform/OpFactory.hpp>

#include <string>
#include <sstream>
#include <exception>

namespace sgpp {
namespace datadriven {
namespace clusteringmpi {

template <class T>
class MPIWorkerPackageBase : virtual public MPIWorkerBase {
 protected:
  int packagesize_multiplier;
  MPI_Datatype mpi_typ;
  bool opencl_node;
  bool overseer_node;

  MPI_Comm &master_worker_comm;
  MPI_Comm &sub_worker_comm;
  bool prefetching;
  int secondary_workpackage[2];

  int opencl_platform = 0;
  int opencl_device = 0;

  void divide_workpackages(int *package, T *erg) {
    // Divide into more work packages
    int packagesize = static_cast<int>(MPIEnviroment::get_configuration()
                                       ["PREFERED_PACKAGESIZE"].getInt());
    T *package_result = new T[packagesize * packagesize_multiplier];
    SimpleQueue<T> workitem_queue(package[0], package[1], packagesize,
                                  sub_worker_comm,
                                  MPIEnviroment::get_sub_worker_count(),
                                  verbose,
                                  prefetching);
    int chunkid = package[0];
    size_t messagesize = 0;
    while (!workitem_queue.is_finished()) {
      // Store result
      messagesize = workitem_queue.receive_result(chunkid, package_result);
      if (verbose) {
        std::cout << "Messagesize: "<< messagesize << std::endl;
        std::cout << package_result[0] << " at  " << chunkid - package[0] + 0
                  << " with packageid " << chunkid << " on "
                  << MPIEnviroment::get_node_rank() << std::endl;
      }
      for (size_t i = 0; i < messagesize; i++) {
        erg[(chunkid - package[0]) * packagesize_multiplier + i] = package_result[i];
      }
    }
    delete [] package_result;
  }

  virtual void receive_and_send_initial_data(void) = 0;
  virtual void begin_opencl_operation(int *workpackage) = 0;
  virtual void finalize_opencl_operation(T *result_buffer, int *workpackage) = 0;

 public:
  MPIWorkerPackageBase(std::string operationName, int multiplier)
      : MPIWorkerBase(operationName),
        opencl_node(false), packagesize_multiplier(multiplier),
        overseer_node(false), master_worker_comm(MPIEnviroment::get_input_communicator()),
        sub_worker_comm(MPIEnviroment::get_communicator()), prefetching(false) {
    if (std::is_same<T, int>::value) {
      mpi_typ = MPI_INT;
    } else if (std::is_same<T, float>::value) {
      mpi_typ = MPI_FLOAT;
    } else if (std::is_same<T, double>::value) {
      mpi_typ = MPI_DOUBLE;
    } else {
      std::stringstream errorString;
      errorString << "Unsupported datatyp in class MPIWorkerPackageBase." << std::endl
                  << "Template class needs to be int, float or double." << std::endl;
      throw std::logic_error(errorString.str());
    }
    if (MPIEnviroment::get_sub_worker_count() > 0) {
      overseer_node = true;
      opencl_node = false;
    } else {
      overseer_node = false;
      opencl_node = true;
    }
    if (MPIEnviroment::get_configuration().contains("PREFETCHING"))
      prefetching = MPIEnviroment::get_configuration()["PREFETCHING"].getBool();
    if (MPIEnviroment::get_configuration().contains("OPENCL_PLATFORM"))
      opencl_platform = MPIEnviroment::get_configuration()["OPENCL_PLATFORM"].getUInt();
    if (MPIEnviroment::get_configuration().contains("OPENCL_DEVICE"))
      opencl_device = MPIEnviroment::get_configuration()["OPENCL_DEVICE"].getUInt();
  }
  virtual ~MPIWorkerPackageBase() {}

  void start_worker_main(void) {
    MPI_Status stat;
    receive_and_send_initial_data();
    // Work Loop
    int datainfo[2];
    T *partial_result = NULL;
    int old_partial_size = 0;
    bool first_package = true;
    do {
      if (!prefetching || first_package) {
        // Receive Workpackage
        MPI_Probe(0, 1, master_worker_comm, &stat);
        MPI_Recv(datainfo, 2, MPI_INT, 0, stat.MPI_TAG, master_worker_comm, &stat);
        if (verbose) {
          std::cout << "Received workpackage [" << datainfo[0] << "," << datainfo[1]
                    << "] on " << MPIEnviroment::get_node_rank() << std::endl;
        }
        first_package = false;
      } else {
        datainfo[0] = secondary_workpackage[0];
        datainfo[1] = secondary_workpackage[1];
      }
      // Check for exit
      if (datainfo[0] == -2 && datainfo[1] == -2) {
        std::cerr << "Node" << MPIEnviroment::get_node_rank()
                  << " received exit signal" << std::endl;
        for (int dest = 1; dest < MPIEnviroment::get_sub_worker_count() + 1; dest++)
          MPI_Send(datainfo, 2, MPI_INT, dest, 1, sub_worker_comm);
        break;
      } else if (datainfo[0] == -1 && datainfo[1] == -1) {
        // Receive exitpackage
        MPI_Probe(0, 1, master_worker_comm, &stat);
        MPI_Recv(secondary_workpackage, 2, MPI_INT, 0, stat.MPI_TAG, master_worker_comm, &stat);
        continue;
      } else {
        if (datainfo[1] * packagesize_multiplier != old_partial_size || partial_result == NULL) {
          if (partial_result != NULL)
            delete [] partial_result;
          partial_result = new T[datainfo[1] * packagesize_multiplier];
          if (verbose)
            std::cout << "New Buffer created!" << std::endl;
          old_partial_size = datainfo[1] * packagesize_multiplier;
        }
        if (opencl_node) {
          // Run partial opencl operation
          begin_opencl_operation(datainfo);
          if (prefetching) {
            // Prefetch secondary workpackage
            MPI_Probe(0, 1, master_worker_comm, &stat);
            MPI_Recv(secondary_workpackage, 2, MPI_INT, 0, stat.MPI_TAG,
                     master_worker_comm, &stat);
            if (verbose) {
              std::cout << "Received secondary workpackage [" << secondary_workpackage[0] << ","
                        << secondary_workpackage[1] << "] on " << MPIEnviroment::get_node_rank()
                        << std::endl;
            }
          }
          // Finish opencl operation
          finalize_opencl_operation(partial_result, datainfo);
        } else {
          if (prefetching) {
            // Prefetch secondary workpackage
            MPI_Probe(0, 1, master_worker_comm, &stat);
            MPI_Recv(secondary_workpackage, 2, MPI_INT, 0, stat.MPI_TAG,
                     master_worker_comm, &stat);
            if (verbose) {
              std::cout << "Received workpackage [" << secondary_workpackage[0] << ","
                        << secondary_workpackage[1] << "] on " << MPIEnviroment::get_node_rank()
                        << std::endl;
            }
          }
          divide_workpackages(datainfo, partial_result);
        }
        std::cerr << "Trying to send back results on "
                  << MPIEnviroment::get_node_rank() << std::endl;
        // Send results back
        MPI_Send(partial_result, datainfo[1] * packagesize_multiplier, mpi_typ,
                 0, 1, master_worker_comm);
      }
    }while(true);
    delete [] partial_result;
  }
};


}  // namespace clusteringmpi
}  // namespace datadriven
}  // namespace sgpp
#endif /* OPERATIONPACKAGEBASEMPI_H */
