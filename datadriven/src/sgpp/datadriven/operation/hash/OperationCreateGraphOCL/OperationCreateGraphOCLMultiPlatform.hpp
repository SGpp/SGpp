// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/tools/SGppStopwatch.hpp>
#include <sgpp/base/exception/operation_exception.hpp>
#include <sgpp/globaldef.hpp>
#include <sgpp/base/operation/hash/OperationMatrix.hpp>
#include <sgpp/base/opencl/OCLOperationConfiguration.hpp>
#include <sgpp/base/opencl/OCLManager.hpp>
#include <sgpp/base/opencl/QueueLoadBalancer.hpp>

#include <chrono>
#include <vector>

#include "OperationCreateGraphOCL.hpp"
#include "KernelCreateGraph.hpp"

namespace SGPP {
namespace datadriven {
namespace DensityOCLMultiPlatform {

template<typename T>
class OperationCreateGraphOCLMultiPlatform : public OperationCreateGraphOCL {
 private:
  size_t dims;
  json::Node &configuration;
  KernelCreateGraph<T> *graph_kernel;
  std::vector<std::shared_ptr<base::OCLDevice>> devices;
  bool verbose;
  std::shared_ptr<base::OCLManagerMultiPlatform> manager;
  std::vector<T> dataVector;

 public:
  OperationCreateGraphOCLMultiPlatform(base::DataMatrix& data, size_t dimensions,
                                       std::shared_ptr<base::OCLManagerMultiPlatform> manager,
                                       json::Node &configuration, size_t k) :
      OperationCreateGraphOCL(), dims(dimensions), configuration(configuration),
      devices(manager->getDevices()), manager(manager), dataVector(data.getSize()) {
    verbose = true;
    double *data_raw = data.getPointer();
    for (size_t i = 0; i < data.getSize(); i++)
      dataVector[i] = static_cast<T>(data_raw[i]);
    graph_kernel = new KernelCreateGraph<T>(devices[0], dims, k, dataVector,
                                            manager, configuration);
  }
  OperationCreateGraphOCLMultiPlatform(double *dataset, size_t datasize, size_t dimensions,
                                       std::shared_ptr<base::OCLManagerMultiPlatform> manager,
                                       json::Node &configuration, size_t k) :
      OperationCreateGraphOCL(), dims(dimensions), configuration(configuration),
      devices(manager->getDevices()), manager(manager), dataVector(datasize) {
    verbose = true;
    for (size_t i = 0; i < datasize; i++)
      dataVector[i] = static_cast<T>(dataset[i]);
    graph_kernel = new KernelCreateGraph<T>(devices[0], dims, k, dataVector,
                                            manager, configuration);
  }

  ~OperationCreateGraphOCLMultiPlatform(void) {
    delete graph_kernel;
  }
  void set_problemchunk(size_t start_id, size_t chunksize) {
    this->start_id = start_id;
    this->chunksize = chunksize;
  }

  void create_graph(std::vector<int> &resultVector, int startid = 0, int chunksize = -1) {
    if (verbose)
      std::cout << "Creating graph for " << dataVector.size() << " datapoints" << std::endl;
    std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();
    try {
      this->graph_kernel->create_graph(resultVector, startid, chunksize);
    }
    catch(base::operation_exception &e) {
      std::cerr << "Error! Could not create graph." << std::endl
                <<"Error Message: " << e.what() << std::endl;
      return;
    }
    end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;

    if (verbose)
      std::cout << "duration create graph" << elapsed_seconds.count() << std::endl;
  }
};

}  // namespace DensityOCLMultiPlatform
}  // namespace datadriven
}  // namespace SGPP
