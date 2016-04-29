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

class OperationDensityMPI : public OperationGridMethod, public base::OperationMatrix {
 public:
  OperationDensityMPI(base::Grid& grid, double lambda) :
      OperationGridMethod(grid, "OperationDensitySlave") {
    // Send lambda to slaves
    for (int dest = 1; dest < MPIEnviroment::get_node_count(); dest++)
      MPI_Send(&lambda, 1, MPI_DOUBLE, dest, 1, MPI_COMM_WORLD);
  }
  virtual ~OperationDensityMPI() {}
  virtual void mult(base::DataVector& alpha, base::DataVector& result) {
    this->start_slave_code();
    // Send alpha vector
    for (int dest = 1; dest < MPIEnviroment::get_node_count(); dest++)
      MPI_Send(alpha.getPointer(), gridsize, MPI_DOUBLE, dest, 1, MPI_COMM_WORLD);
    // Create packages and let the slaves solve them
    double *partial_result = new double[2000];
    SimpleQueue<double> workitem_queue(gridsize, 2000);
    int chunkid = 0;
    int messagesize = workitem_queue.receive_result(chunkid, partial_result);
    while (messagesize > 0) {
      // Store result
      std::cerr << messagesize << std::endl;
      for (int i = 0; i < messagesize; i++) {
        result[chunkid + i] = partial_result[i];
      }
      messagesize = workitem_queue.receive_result(chunkid, partial_result);
    }
    int exitmessage[2] = {-2, -2};
    for (int dest = 1; dest < MPIEnviroment::get_node_count(); dest++)
      MPI_Send(exitmessage, 2, MPI_INT, dest, 1, MPI_COMM_WORLD);
    delete [] partial_result;
  }

  static MPISlaveOperation* create_slave(void) {
    return new OperationDensitySlave();
  }

 private:
  class OperationDensitySlave : public OperationGridMethodSlave {
   private:
    double lambda;
    DensityOCLMultiPlatform::OperationDensityOCL *op;

   public:
    bool verbose;
    OperationDensitySlave()
        : OperationGridMethodSlave() {
      // Receive lambda
      MPI_Probe(0, 1, MPI_COMM_WORLD, &stat);
      MPI_Recv(&lambda, 1, MPI_DOUBLE, stat.MPI_SOURCE, stat.MPI_TAG,
               MPI_COMM_WORLD, &stat);

      // Create opencl operation
      op = createDensityOCLMultiPlatformConfigured(gridpoints, complete_gridsize /
                                                  (2 * grid_dimensions), grid_dimensions,
                                                   lambda, "MyOCLConf.cfg", 0, 0);
    }
    virtual ~OperationDensitySlave() {
      delete op;
    }
    virtual void slave_code(void) {
      MPI_Status stat;

      // Receive alpha vector
      int gridsize = complete_gridsize / (2 * grid_dimensions);
      int buffer_size = 0;
      MPI_Probe(0, 1, MPI_COMM_WORLD, &stat);
      MPI_Get_count(&stat, MPI_DOUBLE, &buffer_size);
      if (buffer_size != gridsize) {
        std::stringstream errorString;
        errorString << "Error: Gridsize " << gridsize << " and the size of the alpha vector "
                    << buffer_size << " should match!" << std::endl;
        throw std::logic_error(errorString.str());
      }
      double *alpha = new double[gridsize];
      MPI_Recv(alpha, gridsize, MPI_DOUBLE, stat.MPI_SOURCE, stat.MPI_TAG,
               MPI_COMM_WORLD, &stat);


      int datainfo[2];
      double *partial_result = NULL;
      int old_partial_size = 0;
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
          if (datainfo[1] != old_partial_size || partial_result == NULL) {
            if (partial_result != NULL)
              delete [] partial_result;
            partial_result = new double[datainfo[1]];
            old_partial_size = datainfo[1];
          }
          // Run partial multiplication
          op->partial_mult(alpha, partial_result, datainfo[0], datainfo[1]);
          // Send results back
          MPI_Send(partial_result, datainfo[1], MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
        }
      }while(true);
      delete [] alpha;
    }
  };
};


}  // namespace clusteringmpi
}  // namespace datadriven
}  // namespace sgpp

#endif /* OPERATIONDENSITYMPI_H */
