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
 public:
  DensityWorker()
      : MPIWorkerGridBase("DensityMultiplicationWorker"), opencl_node(false),
        overseer_node(false), master_worker_comm(MPIEnviroment::get_input_communicator()),
        sub_worker_comm(MPIEnviroment::get_communicator()) {
    if (MPIEnviroment::get_sub_worker_count() > 0) {
      overseer_node = true;
      opencl_node = false;
    } else {
      overseer_node = false;
      opencl_node = true;
    }

    // Receive lambda
    MPI_Status stat;
    MPI_Probe(0, 1, MPI_COMM_WORLD, &stat);
    MPI_Recv(&lambda, 1, MPI_DOUBLE, stat.MPI_SOURCE, stat.MPI_TAG,
             master_worker_comm, &stat);

    // Create opencl operation
    op = createDensityOCLMultiPlatformConfigured(gridpoints, complete_gridsize /
                                                 (2 * grid_dimensions), grid_dimensions,
                                                 lambda, "MyOCLConf.cfg", 0, 0);
    if (verbose) {
      std::cout << "Created opencl density operation on "
                << MPIEnviroment::get_node_rank() << std::endl;
    }
  }
  DensityWorker(base::Grid &grid, double lambda)
      : MPIWorkerGridBase("DensityMultiplicationWorker", grid), opencl_node(false),
        overseer_node(false), master_worker_comm(MPIEnviroment::get_input_communicator()),
        sub_worker_comm(MPIEnviroment::get_communicator()) {
    // Send lambda to slaves
    for (int dest = 1; dest < MPIEnviroment::get_sub_worker_count() + 1; dest++)
      MPI_Send(&lambda, 1, MPI_DOUBLE, dest, 1, sub_worker_comm);
  }
  void start_worker_main(void) {
    // Receive alpha
    MPI_Status stat;
    double *alpha = NULL;
    receive_alpha(alpha);

    // Work Loop
    int datainfo[2];
    double *partial_result = NULL;
    int old_partial_size = 0;
    do {
      // Receive Workpackage
      MPI_Probe(0, 1, master_worker_comm, &stat);
      MPI_Recv(datainfo, 2, MPI_INT, 0, stat.MPI_TAG, master_worker_comm, &stat);
      // Check for exit
      if (datainfo[0] == -2 && datainfo[1] == -2) {
        std::cerr << "Node" << MPIEnviroment::get_node_rank()
                  << " received exit signal" << std::endl;
        for (int dest = 1; dest < MPIEnviroment::get_sub_worker_count() + 1; dest++)
          MPI_Send(datainfo, 2, MPI_INT, dest, 1, sub_worker_comm);
        break;
      } else {
        if (verbose) {
          std::cout << "Node " << MPIEnviroment::get_node_rank()
                    << ": Received work package" << std::endl;
        }
        if (datainfo[1] != old_partial_size || partial_result == NULL) {
          if (partial_result != NULL)
            delete [] partial_result;
          partial_result = new double[datainfo[1]];
          old_partial_size = datainfo[1];
        }
        if (opencl_node) {
          // Run partial multiplication
          op->partial_mult(alpha, partial_result, datainfo[0], datainfo[1]);
        } else {
          // Divide into more work packages
          int packagesize = 1280;
          double *package_result = new double[packagesize];
          SimpleQueue<double> workitem_queue(datainfo[0], datainfo[1], packagesize, sub_worker_comm,
                                             MPIEnviroment::get_sub_worker_count());
          int chunkid = datainfo[0];
          size_t messagesize = workitem_queue.receive_result(chunkid, partial_result);
          while (messagesize > 0) {
            // Store result
            std::cerr << messagesize << std::endl;
            for (size_t i = 0; i < messagesize; i++) {
              partial_result[chunkid - datainfo[0] + i] = package_result[i];
            }
            messagesize = workitem_queue.receive_result(chunkid, partial_result);
          }
          // Send results back
          MPI_Send(partial_result, datainfo[1], MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
        }
      }
    }while(true);
    delete [] alpha;
  }
  virtual ~DensityWorker() {}

 private:
  void send_alpha(double *alpha) {
    // Send alpha vector
    for (int dest = 1; dest < MPIEnviroment::get_sub_worker_count() + 1; dest++)
      MPI_Send(alpha, static_cast<int>(complete_gridsize / (2 * grid_dimensions)),
               MPI_DOUBLE, dest, 1, sub_worker_comm);
  }
  void receive_alpha(double *alpha) {
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
    alpha = new double[gridsize];
    MPI_Recv(alpha, gridsize, MPI_DOUBLE, stat.MPI_SOURCE, stat.MPI_TAG,
             master_worker_comm, &stat);
  }
};

class OperationDensityMultMPI : public base::OperationMatrix {
 public:
  OperationDensityMultMPI(base::Grid &grid, double lambda) {
  }
  virtual ~OperationDensityMultMPI();
  virtual void mult(base::DataVector& alpha, base::DataVector& result) {}
};

}  // namespace clusteringmpi
}  // namespace datadriven
}  // namespace sgpp

#endif /* OPERATIONDENSITYMPI_H */
