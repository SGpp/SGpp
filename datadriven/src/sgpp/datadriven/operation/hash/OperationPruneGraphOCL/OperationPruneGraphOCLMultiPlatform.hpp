// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <omp.h>
#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/operation/hash/OperationMultipleEval.hpp>
#include <sgpp/base/tools/SGppStopwatch.hpp>
#include <sgpp/base/exception/operation_exception.hpp>
#include <sgpp/globaldef.hpp>
#include <sgpp/base/operation/hash/OperationMatrix.hpp>
#include <sgpp/base/opencl/OCLOperationConfiguration.hpp>
#include <sgpp/base/opencl/OCLManager.hpp>
#include <sgpp/base/opencl/QueueLoadBalancer.hpp>

#include <chrono>
#include <vector>

#include "OperationPruneGraphOCL.hpp"
#include "KernelPruneGraph.hpp"

namespace sgpp {
namespace datadriven {
namespace DensityOCLMultiPlatform {

template<typename T>
class OperationPruneGraphOCLMultiPlatform : public OperationPruneGraphOCL {
 private:
  size_t dims;
  size_t gridSize;
  size_t dataSize;
  json::Node &configuration;
  sgpp::datadriven::DensityOCLMultiPlatform::KernelPruneGraph<T> *graph_kernel;
  std::vector<std::shared_ptr<base::OCLDevice>> devices;
  bool verbose;
  std::shared_ptr<base::OCLManagerMultiPlatform> manager;

  std::vector<int> pointsVector;
  std::vector<T> alphaVector;
  std::vector<T> dataVector;

  size_t start_id;
  size_t chunksize;

 public:
  OperationPruneGraphOCLMultiPlatform(base::Grid& grid, base::DataVector& alpha,
                                      base::DataMatrix &data, size_t dims,
                                      std::shared_ptr<base::OCLManagerMultiPlatform> manager,
                                      json::Node &configuration, T treshold, size_t k) :
      OperationPruneGraphOCL(), dims(dims), gridSize(grid.getStorage().getSize()),
      dataSize(data.getSize()), configuration(configuration), devices(manager->getDevices()),
      manager(manager), alphaVector(gridSize), dataVector(data.getSize()) {
    verbose = true;
    // Store Grid in a opencl compatible buffer
    sgpp::base::GridStorage& gridStorage = grid.getStorage();
    size_t pointscount = 0;
    for (int i = 0; i < gridSize; i++) {
      sgpp::base::HashGridIndex *point = gridStorage.get(i);
      pointscount++;
      for (int d = 0; d < dims; d++) {
        pointsVector.push_back(point->getIndex(d));
        pointsVector.push_back(point->getLevel(d));
      }
    }
    if (verbose)
      std::cout << "Grid stored into integer array! Number of gridpoints: "
                << pointscount << std::endl;
    for (size_t i = 0; i < gridSize; i++)
      alphaVector[i] = static_cast<T>(alpha[i]);
    double *data_raw = data.getPointer();
    for (size_t i = 0; i < data.getSize(); i++)
      dataVector[i] = static_cast<T>(data_raw[i]);
    graph_kernel = new KernelPruneGraph<T>(devices[0], dims, treshold, k, manager,
                                           configuration, pointsVector, alphaVector,
                                           dataVector);
  }

  ~OperationPruneGraphOCLMultiPlatform() {
    delete graph_kernel;
  }

  virtual void prune_graph(std::vector<int> &graph) {
    if (verbose)
      std::cout << "Pruning graph for" << graph.size() << " nodes" << std::endl;
    std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();
    this->graph_kernel->prune_graph(graph, 0, 0);
    end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;

    if (verbose)
      std::cout << "duration prune graph: " << elapsed_seconds.count() << std::endl;
  }
};

}  // namespace DensityOCLMultiPlatform
}  // namespace datadriven
}  // namespace sgpp
