// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org
#ifndef OPERATIONDENSITYMPI_H
#define OPERATIONDENSITYMPI_H

#include <mpi.h>

#include <sgpp/datadriven/operation/hash/OperationMPI/OperationGridBaseMPI.hpp>
#include <sgpp/datadriven/operation/hash/OperationDensityOCLMultiPlatform/OpFactory.hpp>

#include <sstream>
#include <exception>

namespace sgpp {
namespace datadriven {
namespace clusteringmpi {
class DensityWorker : public MPIWorkerGridBase {
 protected:
  bool opencl_node;
  bool overseer_node;

  double lambda;
  DensityOCLMultiPlatform::OperationDensityOCL *op;

  MPI_Comm &master_worker_comm;
  MPI_Comm &sub_worker_comm;

  bool prefetching;
  int secondary_workpackage[2];
  bool secondary_package_prefetched;
 public:
  DensityWorker()
      : MPIWorkerBase("DensityMultiplicationWorker"),
        MPIWorkerGridBase("DensityMultiplicationWorker"), opencl_node(false),
        overseer_node(false), master_worker_comm(MPIEnviroment::get_input_communicator()),
        sub_worker_comm(MPIEnviroment::get_communicator()), prefetching(false) {
    if (MPIEnviroment::get_sub_worker_count() > 0) {
      overseer_node = true;
      opencl_node = false;
    } else {
      overseer_node = false;
      opencl_node = true;
    }
    if (MPIEnviroment::get_configuration().contains("PREFETCHING"))
      prefetching = MPIEnviroment::get_configuration()["PREFETCHING"].getBool();

    // Receive lambda
    MPI_Status stat;
    MPI_Probe(0, 1, master_worker_comm, &stat);
    MPI_Recv(&lambda, 1, MPI_DOUBLE, stat.MPI_SOURCE, stat.MPI_TAG,
             master_worker_comm, &stat);

    // Send lambda
    for (int dest = 1; dest < MPIEnviroment::get_sub_worker_count() + 1; dest++)
      MPI_Send(&lambda, 1, MPI_DOUBLE, dest, 1, sub_worker_comm);

    // Create opencl operation
    if (opencl_node) {
      op = createDensityOCLMultiPlatformConfigured(gridpoints, complete_gridsize /
                                                   (2 * grid_dimensions), grid_dimensions,
                                                   lambda, "MyOCLConf.cfg", 0, 0);
    }
    if (verbose) {
      std::cout << "Created mpi opencl density operation on "
                << MPIEnviroment::get_node_rank() << std::endl;
    }
  }
  DensityWorker(base::Grid &grid, double lambda)
      : MPIWorkerBase("DensityMultiplicationWorker"),
        MPIWorkerGridBase("DensityMultiplicationWorker", grid), opencl_node(false),
        overseer_node(false), lambda(lambda),
        master_worker_comm(MPIEnviroment::get_input_communicator()),
        sub_worker_comm(MPIEnviroment::get_communicator()), prefetching(false) {
    if (MPIEnviroment::get_configuration().contains("PREFETCHING"))
      prefetching = MPIEnviroment::get_configuration()["PREFETCHING"].getBool();
    // Send lambda to slaves
    for (int dest = 1; dest < MPIEnviroment::get_sub_worker_count() + 1; dest++)
      MPI_Send(&lambda, 1, MPI_DOUBLE, dest, 1, sub_worker_comm);
    if (verbose) {
      std::cout << "Density master node "
                << MPIEnviroment::get_node_rank() << std::endl;
    }
  }
  void start_worker_main(void) {
    // Receive alpha
    MPI_Status stat;
    double *alpha = NULL;
    receive_alpha(&alpha);
    send_alpha(&alpha);
    if (verbose) {
      std::cout << "Received alpha on " << MPIEnviroment::get_node_rank() << std::endl;
    }

    // Work Loop
    int datainfo[2];
    double *partial_result = NULL;
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
      } else if (datainfo[0] == -1 && datainfo[1] == -1){
        // Receive exitpackage
        MPI_Probe(0, 1, master_worker_comm, &stat);
        MPI_Recv(secondary_workpackage, 2, MPI_INT, 0, stat.MPI_TAG, master_worker_comm, &stat);
        continue;
      } else {
        if (datainfo[1] != old_partial_size || partial_result == NULL) {
          if (partial_result != NULL)
            delete [] partial_result;
          partial_result = new double[datainfo[1]];
          if (verbose)
            std::cout << "New Buffer created!" << std::endl;
          old_partial_size = datainfo[1];
        }
        if (opencl_node) {
          // Run partial multiplication
          op->start_partial_mult(alpha, datainfo[0], datainfo[1]);
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
          // Finish multiplication
          op->finish_partial_mult(partial_result, datainfo[0], datainfo[1]);
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
        // Send results back
        MPI_Send(partial_result, datainfo[1], MPI_DOUBLE, 0, 1, master_worker_comm);
      }
    }while(true);
    delete [] alpha;
    delete [] partial_result;
  }
  virtual ~DensityWorker() {
    if (opencl_node) {
      delete op;
    }
  }

 protected:
  void send_alpha(double **alpha) {
    // Send alpha vector
    for (int dest = 1; dest < MPIEnviroment::get_sub_worker_count() + 1; dest++)
      MPI_Send(*alpha, static_cast<int>(complete_gridsize / (2 * grid_dimensions)),
               MPI_DOUBLE, dest, 1, sub_worker_comm);
  }
  void receive_alpha(double **alpha) {
    // Receive alpha vector
    int gridsize = complete_gridsize / (2 * grid_dimensions);
    int buffer_size = 0;
    MPI_Status stat;
    MPI_Probe(0, 1, master_worker_comm, &stat);
    MPI_Get_count(&stat, MPI_DOUBLE, &buffer_size);
    if (buffer_size != gridsize) {
      std::stringstream errorString;
      errorString << "Error: Gridsize " << gridsize << " and the size of the alpha vector "
                  << buffer_size << " should match!" << std::endl;
      throw std::logic_error(errorString.str());
    }
    *alpha = new double[gridsize];
    MPI_Recv(*alpha, gridsize, MPI_DOUBLE, stat.MPI_SOURCE, stat.MPI_TAG,
             master_worker_comm, &stat);
  }
  void divide_workpackages(int *package, double *erg) {
    // Divide into more work packages
    int packagesize = static_cast<int>(MPIEnviroment::get_configuration()
                                       ["PREFERED_PACKAGESIZE"].getInt());
    double *package_result = new double[packagesize];
    SimpleQueue<double> workitem_queue(package[0], package[1], packagesize,
                                       sub_worker_comm,
                                       MPIEnviroment::get_sub_worker_count(),
                                       verbose,
                                       prefetching);
    int chunkid = package[0];
    size_t messagesize = 0;
    while (!workitem_queue.is_finished()) {
      // Store result
      if (verbose) {
        std::cout << "Messagesize: "<< messagesize << std::endl;
        std::cout << package_result[0] << " at  " << chunkid - package[0] + 0
                  << " with packageid " << chunkid << " on "
                  << MPIEnviroment::get_node_rank() << std::endl;
      }
      messagesize = workitem_queue.receive_result(chunkid, package_result);
      for (size_t i = 0; i < messagesize; i++) {
        erg[chunkid - package[0] + i] = package_result[i];
      }
    }
    delete [] package_result;
  }
};

class OperationDensityMultMPI : public base::OperationMatrix, public DensityWorker {
 public:
  OperationDensityMultMPI(base::Grid &grid, double lambda) :
      MPIWorkerBase("DensityMultiplicationWorker"), DensityWorker(grid, lambda) {
  }
  virtual ~OperationDensityMultMPI() {}
  virtual void mult(base::DataVector& alpha, base::DataVector& result) {
    start_sub_workers();
    double *alpha_ptr = alpha.getPointer();
    send_alpha(&alpha_ptr);
    int datainfo[2];
    datainfo[0] = 0;
    datainfo[1] = static_cast<int>(result.getSize());
    divide_workpackages(datainfo, result.getPointer());
    int exitmessage[2] = {-2, -2};
    for (int dest = 1; dest < MPIEnviroment::get_sub_worker_count() + 1; dest++)
      MPI_Send(exitmessage, 2, MPI_INT, dest, 1, sub_worker_comm);
  }
};

}  // namespace clusteringmpi
}  // namespace datadriven
}  // namespace sgpp

#endif /* OPERATIONDENSITYMPI_H */
