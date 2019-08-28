// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org
#ifndef OPERATIONRHSMPI_H
#define OPERATIONRHSMPI_H

#include <sgpp/datadriven/operation/hash/OperationDensityOCLMultiPlatform/OpFactory.hpp>
#include <sgpp/datadriven/operation/hash/OperationMPI/OperationGraphBaseMPI.hpp>
#include <sgpp/datadriven/operation/hash/OperationMPI/OperationGridBaseMPI.hpp>
#include <sgpp/datadriven/operation/hash/OperationMPI/OperationPackageBaseMPI.hpp>
#include <string>

namespace sgpp {
namespace datadriven {
namespace clusteringmpi {
class DensityRhsWorker : public MPIWorkerGridBase,
                         public MPIWorkerGraphBase,
                         public MPIWorkerPackageBase<double> {
 protected:
  double lambda;
  DensityOCLMultiPlatform::OperationDensity *op;
  base::DataMatrix *data_matrix;

  void receive_and_send_initial_data(void) {
    if (data_matrix != nullptr) delete data_matrix;
    data_matrix = new base::DataMatrix(dataset, dataset_size / dimensions, dimensions);
    if (opencl_node) op->initialize_dataset(*data_matrix);
  }
  void begin_opencl_operation(int *workpackage) {
    op->start_rhs_generation(workpackage[0], workpackage[1]);
  }
  void finalize_opencl_operation(double *result_buffer, int *workpackage) {
    base::DataVector partial_rhs(workpackage[1]);
    op->finalize_rhs_generation(partial_rhs, workpackage[0], workpackage[1]);
    for (int i = 0; i < workpackage[1]; ++i) {
      result_buffer[i] = partial_rhs[i];
    }
  }

 public:
  DensityRhsWorker()
      : MPIWorkerBase("DensityRHSWorker"),
        MPIWorkerGridBase("DensityRHSWorker"),
        MPIWorkerGraphBase("DensityRHSWorker"),
        MPIWorkerPackageBase("DensityMultiplicationWorker", 1) {
    // Create opencl operation
    if (opencl_node) {
      op = createDensityOCLMultiPlatformConfigured(
          gridpoints, complete_gridsize / (2 * grid_dimensions), grid_dimensions, 0.0, parameters,
          opencl_platform, opencl_device);
    }
    data_matrix = nullptr;
  }
  DensityRhsWorker(base::Grid &grid, sgpp::base::DataMatrix &data, std::string ocl_config_file)
      : MPIWorkerBase("DensityRHSWorker"),
        MPIWorkerGridBase("DensityRHSWorker", grid),
        MPIWorkerGraphBase("DensityRHSWorker", data, 0),
        MPIWorkerPackageBase("DensityMultiplicationWorker", 1, ocl_config_file) {
    data_matrix = nullptr;
  }
  virtual ~DensityRhsWorker(void) {
    if (opencl_node) delete op;
  }
};
class OperationDensityRhsMPI : public DensityRhsWorker {
 public:
  OperationDensityRhsMPI(base::Grid &grid, sgpp::base::DataMatrix &data,
                         std::string ocl_config_file)
      : MPIWorkerBase("DensityRHSWorker"), DensityRhsWorker(grid, data, ocl_config_file) {}
  virtual ~OperationDensityRhsMPI() {}
  virtual void generate_b(base::DataVector &b) {
    DensityRhsWorker::MPIWorkerGridBase::start_sub_workers();
    int datainfo[2];
    datainfo[0] = 0;
    datainfo[1] = static_cast<int>(b.getSize());
    divide_workpackages(datainfo, b.getPointer());
    int exitmessage[2] = {-2, -2};
    for (int dest = 1; dest < MPIEnviroment::get_sub_worker_count() + 1; dest++)
      MPI_Send(exitmessage, 2, MPI_INT, dest, 1, sub_worker_comm);
  }
};
}  // namespace clusteringmpi
}  // namespace datadriven
}  // namespace sgpp
#endif /* OPERATIONRHSMPI_H */
