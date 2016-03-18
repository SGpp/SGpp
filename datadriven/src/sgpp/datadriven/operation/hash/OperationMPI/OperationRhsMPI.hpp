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




/// MPI operation for creating right hand side clustering vector
class OperationRhsMPI : public OperationGridMethod, public OperationGraphMethodMPI {
 public:
  OperationRhsMPI(base::Grid& grid, size_t dimensions, base::DataMatrix &data)
      : OperationGridMethod(grid, "OperationRhsSlave"),
        OperationGraphMethodMPI(data, dimensions, 12) {
  }
  base::DataVector create_rhs(void) {
    // Start slave operations
    this->OperationGridMethod::start_slave_code();

    // Return vector
    base::DataVector ret_rhs(gridsize);

    // Create packages and let the slaves solve them
    double *partial_result = new double[2000];
    SimpleQueue<double> workitem_queue(gridsize, 2000);
    int chunkid = 0;
    int messagesize = workitem_queue.receive_result(chunkid, partial_result);
    while (messagesize > 0) {
      // Store result
      std::cerr << messagesize << std::endl;
      for (int i = 0; i < messagesize; i++) {
        ret_rhs[chunkid + i] = partial_result[i];
      }
      messagesize = workitem_queue.receive_result(chunkid, partial_result);
    }
    int exitmessage[2] = {-2, -2};
    for (int dest = 1; dest < MPIEnviroment::get_node_count(); dest++)
      MPI_Send(exitmessage, 2, MPI_INT, dest, 1, MPI_COMM_WORLD);
    delete [] partial_result;
    return ret_rhs;
  }
  virtual ~OperationRhsMPI() {}
  static MPISlaveOperation* create_slave(void) {
    return new OperationRhsSlave();
  }

 private:
  class OperationRhsSlave : public OperationGridMethodSlave,
                            public OperationGraphMethodSlave {
   public:
    bool verbose;
    OperationRhsSlave() : OperationGridMethodSlave(),
                          OperationGraphMethodSlave() {
      verbose = true;
    }
    virtual ~OperationRhsSlave() {}
    virtual void slave_code(void) {
      DensityOCLMultiPlatform::OperationDensityOCL *op =
          sgpp::datadriven::createDensityOCLMultiPlatformConfigured(gridpoints, complete_gridsize /
                                                                    (2 * grid_dimensions), dimensions,
                                                                    0.0, "MyOCLConf.cfg");

      base::DataMatrix data_matrix(dataset, dataset_size / dimensions, dimensions);
      int datainfo[2];
      do {
        MPI_Probe(0, 1, MPI_COMM_WORLD, &stat);
        MPI_Recv(datainfo, 2, MPI_INT, 0, stat.MPI_TAG, MPI_COMM_WORLD, &stat);
        // Check for exit
        if (datainfo[0] == -2 && datainfo[1] == -2) {
          std::cerr << "Node" << MPIEnviroment::get_node_rank()
                    << " received exit signal" << std::endl;
          break;
        } else {
          if (verbose) {
            std::cout << "Node " << MPIEnviroment::get_node_rank()
                      << ": Received work package" << std::endl;
          }
          // Create partial graph
          base::DataVector partial_rhs(datainfo[1]);
          op->generateb(data_matrix, partial_rhs, datainfo[0], datainfo[1]);
          // Send results back
          MPI_Send(partial_rhs.getPointer(), datainfo[1], MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
        }
      }while(true);
      delete op;
    }
  };
};


}  // namespace clusteringmpi
}  // namespace datadriven
}  // namespace sgpp
#endif /* OPERATIONRHSMPI_H */
