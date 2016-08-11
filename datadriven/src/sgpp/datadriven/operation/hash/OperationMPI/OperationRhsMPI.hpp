// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org
#ifndef OPERATIONRHSMPI_H
#define OPERATIONRHSMPI_H

#include <sgpp/datadriven/operation/hash/OperationMPI/OperationGridBaseMPI.hpp>
#include <sgpp/datadriven/operation/hash/OperationMPI/OperationGraphBaseMPI.hpp>
#include <sgpp/datadriven/operation/hash/OperationDensityOCLMultiPlatform/OpFactory.hpp>
#include <string>

namespace sgpp {
namespace datadriven {
namespace clusteringmpi {
class DensityRhsWorker : public MPIWorkerGridBase, public MPIWorkerGraphBase {
 protected:
  bool opencl_node;
  bool overseer_node;
  double lambda;
  DensityOCLMultiPlatform::OperationDensityOCL *op;

  MPI_Comm &master_worker_comm;
  MPI_Comm &sub_worker_comm;
  void divide_workpackages(int *package, double *erg) {
    // Divide into more work packages
    int packagesize = static_cast<int>(MPIEnviroment::get_configuration()
                                       ["PREFERED_PACKAGESIZE"].getInt());
    double *package_result = new double[packagesize];
    SimpleQueue<double> workitem_queue(package[0], package[1], packagesize,
                                       sub_worker_comm,
                                       MPIEnviroment::get_sub_worker_count());
    int chunkid = package[0];
    size_t messagesize = 0;
    while (!workitem_queue.is_finished()) {
      // Store result
      if (MPIWorkerGridBase::verbose) {
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

 public:
  DensityRhsWorker()
      : MPIWorkerBase("DensityRHSWorker"),
        MPIWorkerGridBase("DensityRHSWorker"), MPIWorkerGraphBase("DensityRHSWorker"),
        opencl_node(false), overseer_node(false),
        master_worker_comm(MPIEnviroment::get_input_communicator()),
        sub_worker_comm(MPIEnviroment::get_communicator()) {
    if (MPIEnviroment::get_sub_worker_count() > 0) {
      overseer_node = true;
      opencl_node = false;
    } else {
      overseer_node = false;
      opencl_node = true;
    }
    // Create opencl operation
    op = createDensityOCLMultiPlatformConfigured(gridpoints, complete_gridsize /
                                                 (2 * grid_dimensions), grid_dimensions,
                                                 0.0, "MyOCLConf.cfg", 0, 0);
  }
  DensityRhsWorker(base::Grid &grid, sgpp::base::DataMatrix &data)
      : MPIWorkerBase("DensityRHSWorker"),
        MPIWorkerGridBase("DensityRHSWorker", grid),
        MPIWorkerGraphBase("DensityRHSWorker", data, 0),
        opencl_node(false), overseer_node(false),
        master_worker_comm(MPIEnviroment::get_input_communicator()),
        sub_worker_comm(MPIEnviroment::get_communicator()) {}
  void start_worker_main(void) {
    base::DataMatrix data_matrix(dataset, dataset_size / dimensions, dimensions);
    // Receive alpha
    MPI_Status stat;

    // Work Loop
    int datainfo[2];
    double *partial_result = NULL;
    int old_partial_size = 0;
    do {
      // Receive Workpackage
      MPI_Probe(0, 1, master_worker_comm, &stat);
      MPI_Recv(datainfo, 2, MPI_INT, 0, stat.MPI_TAG, master_worker_comm, &stat);
      if (MPIWorkerGridBase::verbose) {
        std::cout << "Received workpackage [" << datainfo[0] << "," << datainfo[1]
                  << "] on " << MPIEnviroment::get_node_rank() << std::endl;
      }
      // Check for exit
      if (datainfo[0] == -2 && datainfo[1] == -2) {
        std::cerr << "Node" << MPIEnviroment::get_node_rank()
                  << " received exit signal" << std::endl;
        for (int dest = 1; dest < MPIEnviroment::get_sub_worker_count() + 1; dest++)
          MPI_Send(datainfo, 2, MPI_INT, dest, 1, sub_worker_comm);
        break;
      } else {
        if (datainfo[1] != old_partial_size || partial_result == NULL) {
          if (partial_result != NULL)
            delete [] partial_result;
          partial_result = new double[datainfo[1]];
          if (MPIWorkerGridBase::verbose)
            std::cout << "New Buffer created!" << std::endl;
          old_partial_size = datainfo[1];
        }
        if (opencl_node) {
          // Create partial rhs
          base::DataVector partial_rhs(datainfo[1]);
          op->generateb(data_matrix, partial_rhs, datainfo[0], datainfo[1]);
          for (auto i = 0; i < datainfo[1]; ++i) {
            partial_result[i] = partial_rhs[i];
          }

          if (MPIWorkerGridBase::verbose)
            std::cerr << "Workpackage abgeschlossen" << std::endl;
        } else {
          divide_workpackages(datainfo, partial_result);
        }
        // Send results back
        MPI_Send(partial_result, datainfo[1], MPI_DOUBLE, 0, 1, master_worker_comm);
      }
    }while(true);
  }
};
class OperationDensityRhsMPI : public DensityRhsWorker {
 public:
  OperationDensityRhsMPI(base::Grid &grid, sgpp::base::DataMatrix &data) :
      MPIWorkerBase("DensityRHSWorker"), DensityRhsWorker(grid, data) {
  }
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
