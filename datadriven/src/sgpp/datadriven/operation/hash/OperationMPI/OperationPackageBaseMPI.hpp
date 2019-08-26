
// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org
#ifndef OPERATIONPACKAGEBASEMPI_H
#define OPERATIONPACKAGEBASEMPI_H

#include <mpi.h>

#include <sgpp/datadriven/operation/hash/OperationDensityOCLMultiPlatform/OpFactory.hpp>
#include <sgpp/datadriven/operation/hash/OperationMPI/MPIEnviroment.hpp>
#include <sgpp/datadriven/operation/hash/OperationMPI/OperationMPI.hpp>

#include <algorithm>
#include <exception>
#include <sstream>
#include <string>

namespace sgpp {
namespace datadriven {
namespace clusteringmpi {

template <class T>
class MPIWorkerPackageBase : virtual public MPIWorkerBase {
 protected:
  bool opencl_node;
  int packagesize_multiplier;
  bool overseer_node;
  MPI_Comm &master_worker_comm;
  MPI_Comm &sub_worker_comm;
  bool prefetching;
  bool redistribute;
  base::OCLOperationConfiguration *parameters;

  MPI_Datatype mpi_typ;
  int secondary_workpackage[2];
  size_t size;
  size_t opencl_platform = 0;
  size_t opencl_device = 0;

  void divide_workpackages(int *package, T *erg) {
    // Divide into more work packages
    size_t packagesize = size;
    if (redistribute) {
      size_t logical_package_count =
          (package[1] / (packagesize * MPIEnviroment::get_sub_worker_count()));
      if (logical_package_count < 1) logical_package_count = 1;
      packagesize += (package[1] % (packagesize * MPIEnviroment::get_sub_worker_count())) /
                     (MPIEnviroment::get_sub_worker_count() * (logical_package_count));
    }
    T *package_result = new T[packagesize * packagesize_multiplier];
    SimpleQueue<T> workitem_queue(package[0], package[1], packagesize, sub_worker_comm,
                                  MPIEnviroment::get_sub_worker_count(), verbose, prefetching);
    int chunkid = package[0];
    size_t messagesize = 0;
    while (!workitem_queue.is_finished()) {
      // Store result
      messagesize = workitem_queue.receive_result(chunkid, package_result);
      if (verbose) {
        std::cout << "Messagesize: " << messagesize << std::endl;
        std::cout << package_result[0] << " at  " << chunkid - package[0] + 0 << " with packageid "
                  << chunkid << " on " << MPIEnviroment::get_node_rank() << std::endl;
      }
      // for (size_t i = 0; i < messagesize; i++) {
      //   erg[(chunkid - package[0]) * packagesize_multiplier + i] = package_result[i];
      // }
      memcpy(erg + (chunkid - package[0]) * packagesize_multiplier, package_result,
             sizeof(T) * messagesize);
    }
    delete[] package_result;
  }

  virtual void receive_and_send_initial_data(void) = 0;
  virtual void begin_opencl_operation(int *workpackage) = 0;
  virtual void finalize_opencl_operation(T *result_buffer, int *workpackage) = 0;

 public:
  MPIWorkerPackageBase(std::string operationName, int multiplier)
      : MPIWorkerBase(operationName),
        opencl_node(false),
        packagesize_multiplier(multiplier),
        overseer_node(false),
        master_worker_comm(MPIEnviroment::get_input_communicator()),
        sub_worker_comm(MPIEnviroment::get_communicator()),
        prefetching(false),
        redistribute(true),
        parameters(nullptr) {
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
    if (MPIEnviroment::get_configuration().contains("PREFERED_PACKAGESIZE"))
      size = MPIEnviroment::get_configuration()["PREFERED_PACKAGESIZE"].getUInt();
    else
      size = 2048;
    if (MPIEnviroment::get_configuration().contains("REDISTRIBUTE"))
      redistribute = MPIEnviroment::get_configuration()["REDISTRIBUTE"].getBool();
    if (MPIEnviroment::get_configuration().contains("PREFETCHING"))
      prefetching = MPIEnviroment::get_configuration()["PREFETCHING"].getBool();
    if (MPIEnviroment::get_configuration().contains("OPENCL_PLATFORM"))
      opencl_platform = MPIEnviroment::get_configuration()["OPENCL_PLATFORM"].getUInt();
    if (MPIEnviroment::get_configuration().contains("OPENCL_DEVICE"))
      opencl_device = MPIEnviroment::get_configuration()["OPENCL_DEVICE"].getUInt();

    // receive opencl configuration
    std::cerr << "waiting for opencl configuration on" << MPIEnviroment::get_node_rank()
              << std::endl;
    MPI_Status stat;
    int messagesize = 0;
    MPI_Probe(0, 1, master_worker_comm, &stat);
    MPI_Get_count(&stat, MPI_CHAR, &messagesize);
    char *serialized_conf = new char[messagesize];
    MPI_Recv(serialized_conf, messagesize, MPI_CHAR, stat.MPI_SOURCE, stat.MPI_TAG,
             master_worker_comm, &stat);
    // Send OCL configuration
    for (int dest = 1; dest < MPIEnviroment::get_sub_worker_count() + 1; dest++)
      MPI_Send(serialized_conf, messagesize, MPI_CHAR, dest, 1, sub_worker_comm);
    parameters = new base::OCLOperationConfiguration();
    parameters->deserialize(serialized_conf);
    delete[] serialized_conf;
    std::cerr << "Received ocl config on" << MPIEnviroment::get_node_rank() << std::endl;
  }
  MPIWorkerPackageBase(std::string operationName, int multiplier, std::string ocl_conf_filename)
      : MPIWorkerBase(operationName),
        opencl_node(false),
        packagesize_multiplier(multiplier),
        overseer_node(false),
        master_worker_comm(MPIEnviroment::get_input_communicator()),
        sub_worker_comm(MPIEnviroment::get_communicator()),
        prefetching(false),
        redistribute(true),
        parameters(nullptr) {
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
    if (MPIEnviroment::get_configuration().contains("PREFERED_PACKAGESIZE"))
      size = MPIEnviroment::get_configuration()["PREFERED_PACKAGESIZE"].getUInt();
    else
      size = 2048;
    if (MPIEnviroment::get_configuration().contains("REDISTRIBUTE"))
      redistribute = MPIEnviroment::get_configuration()["REDISTRIBUTE"].getBool();
    if (MPIEnviroment::get_configuration().contains("PREFETCHING"))
      prefetching = MPIEnviroment::get_configuration()["PREFETCHING"].getBool();
    if (MPIEnviroment::get_configuration().contains("OPENCL_PLATFORM"))
      opencl_platform = MPIEnviroment::get_configuration()["OPENCL_PLATFORM"].getUInt();
    if (MPIEnviroment::get_configuration().contains("OPENCL_DEVICE"))
      opencl_device = MPIEnviroment::get_configuration()["OPENCL_DEVICE"].getUInt();

    std::cerr << "Start sending ocl conf" << MPIEnviroment::get_node_rank() << std::endl;
    base::OCLOperationConfiguration *parameters =
        new base::OCLOperationConfiguration(ocl_conf_filename);
    std::ostringstream sstream;
    parameters->serialize(sstream, 0);
    std::string serialized_conf = sstream.str();
    char *conf_message = new char[serialized_conf.length() + 1];
    std::copy(serialized_conf.begin(), serialized_conf.end(), conf_message);
    conf_message[serialized_conf.size()] = '\0';
    // Send OCL configuration
    for (int dest = 1; dest < MPIEnviroment::get_sub_worker_count() + 1; dest++)
      MPI_Send(conf_message, static_cast<int>(serialized_conf.size() + 1), MPI_CHAR, dest, 1,
               sub_worker_comm);
    delete[] conf_message;
  }
  virtual ~MPIWorkerPackageBase() {}

  void start_worker_main(void) {
    MPI_Status stat;
    MPI_Request request[2];
    receive_and_send_initial_data();
    // Work Loop
    int datainfo[2];
    T *partial_results[2];
    partial_results[0] = nullptr;
    partial_results[1] = nullptr;
    int buffersizes[2];
    buffersizes[0] = 0;
    buffersizes[1] = 0;
    bool first_package = true;
    unsigned int currentbuffer = 0;
    do {
      if (!prefetching || first_package) {
        // Receive Workpackage
        MPI_Probe(0, 1, master_worker_comm, &stat);
        MPI_Recv(datainfo, 2, MPI_INT, 0, stat.MPI_TAG, master_worker_comm, &stat);
        if (verbose) {
          std::cout << "Received workpackage [" << datainfo[0] << "," << datainfo[1] << "] on "
                    << MPIEnviroment::get_node_rank() << std::endl;
        }
        first_package = false;
      } else {
        datainfo[0] = secondary_workpackage[0];
        datainfo[1] = secondary_workpackage[1];
      }
      // Check for exit
      if (datainfo[0] == -2 && datainfo[1] == -2) {
        std::cerr << "Node" << MPIEnviroment::get_node_rank() << " received exit signal"
                  << std::endl;
        for (int dest = 1; dest < MPIEnviroment::get_sub_worker_count() + 1; dest++)
          MPI_Send(datainfo, 2, MPI_INT, dest, 1, sub_worker_comm);
        break;
      } else if (datainfo[0] == -1 && datainfo[1] == -1) {
        // Receive exitpackage
        MPI_Probe(0, 1, master_worker_comm, &stat);
        MPI_Recv(secondary_workpackage, 2, MPI_INT, 0, stat.MPI_TAG, master_worker_comm, &stat);
        continue;
      } else {
        if (datainfo[1] * packagesize_multiplier != buffersizes[currentbuffer] ||
            partial_results[currentbuffer] == nullptr) {
          if (partial_results[currentbuffer] != nullptr) delete[] partial_results[currentbuffer];
          partial_results[currentbuffer] = new T[datainfo[1] * packagesize_multiplier];
          if (verbose) std::cout << "New Buffer created!" << std::endl;
          buffersizes[currentbuffer] = datainfo[1] * packagesize_multiplier;
        }
        if (opencl_node) {
          // Run partial opencl operation
          begin_opencl_operation(datainfo);
          if (prefetching) {
            // Prefetch secondary workpackage
            MPI_Probe(0, 1, master_worker_comm, &stat);
            MPI_Recv(secondary_workpackage, 2, MPI_INT, 0, stat.MPI_TAG, master_worker_comm, &stat);
            if (verbose) {
              std::cout << "Received secondary workpackage [" << secondary_workpackage[0] << ","
                        << secondary_workpackage[1] << "] on " << MPIEnviroment::get_node_rank()
                        << std::endl;
            }
          }
          // Finish opencl operation
          finalize_opencl_operation(partial_results[currentbuffer], datainfo);
        } else {
          if (prefetching) {
            // Prefetch secondary workpackage
            MPI_Probe(0, 1, master_worker_comm, &stat);
            MPI_Recv(secondary_workpackage, 2, MPI_INT, 0, stat.MPI_TAG, master_worker_comm, &stat);
            if (verbose) {
              std::cout << "Received workpackage [" << secondary_workpackage[0] << ","
                        << secondary_workpackage[1] << "] on " << MPIEnviroment::get_node_rank()
                        << std::endl;
            }
          }
          divide_workpackages(datainfo, partial_results[currentbuffer]);
        }
        std::cout << "Trying to send back results on " << MPIEnviroment::get_node_rank()
                  << std::endl;
        // Send results back
        MPI_Isend(partial_results[currentbuffer], datainfo[1] * packagesize_multiplier, mpi_typ, 0,
                  1, master_worker_comm, request + currentbuffer);
        currentbuffer = (currentbuffer + 1) % 2;
        // Wait for the old message beofre reusing the buffer
        if (buffersizes[currentbuffer] != 0) MPI_Wait(request + currentbuffer, &stat);
      }
    } while (true);
    if (partial_results[0] != nullptr) delete[] partial_results[0];
    if (partial_results[1] != nullptr) delete[] partial_results[1];
  }
};

}  // namespace clusteringmpi
}  // namespace datadriven
}  // namespace sgpp
#endif /* OPERATIONPACKAGEBASEMPI_H */
