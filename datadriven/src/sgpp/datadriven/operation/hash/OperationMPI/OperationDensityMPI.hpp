// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef OPERATIONDENSITYMPI_H
#define OPERATIONDENSITYMPI_H

#include <mpi.h>

#include <sgpp/datadriven/operation/hash/OperationDensityOCLMultiPlatform/OpFactory.hpp>
#include <sgpp/datadriven/operation/hash/OperationMPI/OperationGridBaseMPI.hpp>
#include <sgpp/datadriven/operation/hash/OperationMPI/OperationPackageBaseMPI.hpp>

#include <exception>
#include <sstream>
#include <string>

namespace sgpp {
namespace datadriven {
namespace clusteringmpi {
class DensityWorker : public MPIWorkerGridBase, public MPIWorkerPackageBase<double> {
 protected:
  double lambda;
  DensityOCLMultiPlatform::OperationDensity *op;
  double *alpha;
  size_t oldgridsize;

  void receive_and_send_initial_data(void) {
    receive_alpha(&alpha);
    if (opencl_node) op->initialize_alpha(alpha);
    send_alpha(&alpha);
    if (verbose) {
      std::cout << "Received alpha on " << MPIEnviroment::get_node_rank() << std::endl;
    }
  }
  void begin_opencl_operation(int *workpackage) {
    op->start_partial_mult(workpackage[0], workpackage[1]);
  }
  void finalize_opencl_operation(double *result_buffer, int *workpackage) {
    op->finish_partial_mult(result_buffer, workpackage[0], workpackage[1]);
  }

 public:
  DensityWorker()
      : MPIWorkerBase("DensityMultiplicationWorker"),
        MPIWorkerGridBase("DensityMultiplicationWorker"),
        MPIWorkerPackageBase("DensityMultiplicationWorker", 1) {
    alpha = nullptr;
    oldgridsize = 0;
    // Receive lambda
    MPI_Status stat;
    MPI_Probe(0, 1, master_worker_comm, &stat);
    MPI_Recv(&lambda, 1, MPI_DOUBLE, stat.MPI_SOURCE, stat.MPI_TAG, master_worker_comm, &stat);

    // Send lambda
    for (int dest = 1; dest < MPIEnviroment::get_sub_worker_count() + 1; dest++)
      MPI_Send(&lambda, 1, MPI_DOUBLE, dest, 1, sub_worker_comm);

    // Create opencl operation
    if (opencl_node) {
      op = createDensityOCLMultiPlatformConfigured(
          gridpoints, complete_gridsize / (2 * grid_dimensions), grid_dimensions, lambda,
          parameters, opencl_platform, opencl_device);
    }
    if (verbose) {
      std::cout << "Created mpi opencl density operation on " << MPIEnviroment::get_node_rank()
                << std::endl;
    }
  }
  DensityWorker(base::Grid &grid, double lambda)
      : MPIWorkerBase("DensityMultiplicationWorker"),
        MPIWorkerGridBase("DensityMultiplicationWorker", grid),
        MPIWorkerPackageBase("DensityMultiplicationWorker", 1) {
    alpha = nullptr;
    // Send lambda to slaves
    for (int dest = 1; dest < MPIEnviroment::get_sub_worker_count() + 1; dest++)
      MPI_Send(&lambda, 1, MPI_DOUBLE, dest, 1, sub_worker_comm);
    if (verbose) {
      std::cout << "Density master node " << MPIEnviroment::get_node_rank() << std::endl;
    }
  }
  DensityWorker(base::Grid &grid, double lambda, std::string ocl_conf_filename)
      : MPIWorkerBase("DensityMultiplicationWorker"),
        MPIWorkerGridBase("DensityMultiplicationWorker", grid),
        MPIWorkerPackageBase("DensityMultiplicationWorker", 1, ocl_conf_filename) {
    alpha = nullptr;
    // Send lambda to slaves
    for (int dest = 1; dest < MPIEnviroment::get_sub_worker_count() + 1; dest++)
      MPI_Send(&lambda, 1, MPI_DOUBLE, dest, 1, sub_worker_comm);
    if (verbose) {
      std::cout << "Density master node " << MPIEnviroment::get_node_rank() << std::endl;
    }
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
      MPI_Send(*alpha, static_cast<int>(complete_gridsize / (2 * grid_dimensions)), MPI_DOUBLE,
               dest, 1, sub_worker_comm);
  }
  void receive_alpha(double **alpha) {
    // Receive alpha vector
    size_t gridsize = complete_gridsize / (2 * grid_dimensions);
    int buffer_size = 0;
    MPI_Status stat;
    MPI_Probe(0, 1, master_worker_comm, &stat);

    MPI_Get_count(&stat, MPI_DOUBLE, &buffer_size);
    if (static_cast<size_t>(buffer_size) != gridsize) {
      std::stringstream errorString;
      errorString << "Error: Gridsize " << gridsize << " and the size of the alpha vector "
                  << buffer_size << " should match!" << std::endl;
      throw std::logic_error(errorString.str());
    }
    if (gridsize != oldgridsize) {
      if (*alpha != nullptr) delete[](*alpha);
      *alpha = new double[gridsize];
      oldgridsize = gridsize;
    }
    MPI_Recv(*alpha, static_cast<int>(gridsize), MPI_DOUBLE, stat.MPI_SOURCE, stat.MPI_TAG,
             master_worker_comm, &stat);
  }
};

class OperationDensityMultMPI : public base::OperationMatrix, public DensityWorker {
 public:
  OperationDensityMultMPI(base::Grid &grid, double lambda, std::string ocl_conf_filename)
      : MPIWorkerBase("DensityMultiplicationWorker"),
        DensityWorker(grid, lambda, ocl_conf_filename) {}
  virtual ~OperationDensityMultMPI() {}
  virtual void mult(base::DataVector &alpha, base::DataVector &result) {
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
