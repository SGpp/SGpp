// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org
#pragma once

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/tools/SGppStopwatch.hpp>
#include <sgpp/base/exception/operation_exception.hpp>
#include <sgpp/globaldef.hpp>
#include <sgpp/base/operation/hash/OperationMatrix.hpp>
#include <sgpp/base/opencl/OCLOperationConfiguration.hpp>
#include <sgpp/base/opencl/OCLManager.hpp>
#include <sgpp/base/opencl/QueueLoadBalancer.hpp>

#include <omp.h>
#include <chrono>
#include <vector>

#include "OperationDensityOCL.hpp"
#include "KernelMult.hpp"
#include "KernelB.hpp"

namespace sgpp {
namespace datadriven {
namespace DensityOCLMultiPlatform {

template<typename T>
class OperationDensityOCLMultiPlatform: public OperationDensityOCL {
 private:
  size_t dims;
  size_t gridSize;
  json::Node &firstKernelConfig;
  json::Node &secondKernelConfig;
  sgpp::datadriven::DensityOCLMultiPlatform::KernelDensityMult<T> *multKernel;
  sgpp::datadriven::DensityOCLMultiPlatform::KernelDensityB<T> *bKernel;
  std::vector<std::shared_ptr<base::OCLDevice>> devices;
  bool verbose;
  std::vector<int> points;
  std::shared_ptr<base::OCLManagerMultiPlatform> manager;
  T lambda;
  size_t start_id;
  size_t chunksize;

 public:
  OperationDensityOCLMultiPlatform(base::Grid& grid, size_t dimensions,
                                   std::shared_ptr<base::OCLManagerMultiPlatform> manager,
                                   json::Node &firstKernelConfig, json::Node &secondKernelConfig,
                                   T lambda) : OperationDensityOCL(), dims(dimensions),
                                               gridSize(grid.getStorage().getSize()),
                                               firstKernelConfig(firstKernelConfig),
                                               secondKernelConfig(secondKernelConfig),
                                               devices(manager->getDevices()),
                                               manager(manager), lambda(lambda), start_id(0),
                                               chunksize(-1) {
    verbose = true;
    // Store Grid in a opencl compatible buffer
    sgpp::base::GridStorage& gridStorage = grid.getStorage();
    size_t pointscount = 0;
    for (int i = 0; i < gridSize; i++) {
      sgpp::base::HashGridIndex *point = gridStorage.get(i);
      pointscount++;
      for (int d = 0; d < dims; d++) {
        points.push_back(point->getIndex(d));
        points.push_back(point->getLevel(d));
      }
    }
    if (verbose)
      std::cout << "Grid stored into integer array! Number of gridpoints: "
                << pointscount << std::endl;
    multKernel = new KernelDensityMult<T>(devices[0], dims, manager, firstKernelConfig,
                                          points, lambda);
    bKernel = new KernelDensityB<T>(devices[0], dims, manager, secondKernelConfig,
                                    points);
  }
  OperationDensityOCLMultiPlatform(int *gridpoints, size_t gridsize, size_t dimensions,
                                   std::shared_ptr<base::OCLManagerMultiPlatform> manager,
                                   json::Node &firstKernelConfig, json::Node &secondKernelConfig,
                                   T lambda) : OperationDensityOCL(), dims(dimensions),
                                               gridSize(gridsize),
                                               firstKernelConfig(firstKernelConfig),
                                               secondKernelConfig(secondKernelConfig),
                                               devices(manager->getDevices()),
                                               manager(manager), lambda(lambda), start_id(0),
                                               chunksize(-1) {
    verbose = true;
    for (int i = 0; i < gridSize; i++) {
      for (int d = 0; d < dims; d++) {
        points.push_back(gridpoints[2 * dimensions * i + 2 * d]);
        points.push_back(gridpoints[2 * dimensions * i + 2 * d + 1]);
      }
    }
    multKernel = new KernelDensityMult<T>(devices[0], dims, manager, firstKernelConfig,
                                          points, lambda);
    bKernel = new KernelDensityB<T>(devices[0], dims, manager, secondKernelConfig,
                                    points);
  }

  ~OperationDensityOCLMultiPlatform() {
    delete multKernel;
    delete bKernel;
  }

  void set_problemchunk(size_t start_id, size_t chunksize) {
    this->start_id = start_id;
    this->chunksize = chunksize;
  }
  void partial_mult(double *alpha, double *result, size_t start_id, size_t chunksize) {
    std::vector<T> alphaVector(gridSize);
    std::vector<T> resultVector(chunksize);
    for (size_t i = 0; i < chunksize; i++) {
      resultVector[i] =(result[i]);
    }
    for (auto i = 0; i < gridSize; ++i) {
      alphaVector[i] = alpha[i];
    }

    if (verbose)
      std::cout << "starting multiplication with " << gridSize << " entries" << std::endl;
    std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();
    this->multKernel->mult(alphaVector, resultVector, start_id, chunksize);
    end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;

    if (verbose) {
      std::cout << "duration mult ocl: " << elapsed_seconds.count() << std::endl;
    }
    for (size_t i = 0; i < chunksize; i++)
      result[i] = resultVector[i];
  }

  void mult(base::DataVector& alpha, base::DataVector& result) override {
    std::vector<T> alphaVector(gridSize);
    std::vector<T> resultVector(gridSize);
    for (size_t i = 0; i < gridSize; i++) {
      alphaVector[i] = static_cast<T>(alpha[i]);
      resultVector[i] = static_cast<T>(result[i]);
    }
    if (verbose)
      std::cout << "starting multiplication with " << gridSize << " entries" << std::endl;
    std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();
    this->multKernel->mult(alphaVector, resultVector, 0, -1);
    end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;

    if (verbose) {
      std::cout << "duration mult ocl: " << elapsed_seconds.count() << std::endl;
    }
    for (size_t i = 0; i < gridSize; i++)
      result[i] = resultVector[i];
  }

  void generateb(base::DataMatrix &dataset, sgpp::base::DataVector &b,
                 size_t start_id = 0, size_t chunksize = -1) {
    std::chrono::time_point<std::chrono::system_clock> start, end;
    if (verbose) {
      if (chunksize == -1)
        std::cout << "starting rhs kernel methode! datasize: " << b.getSize() << std::endl;
      else
        std::cout << "starting rhs kernel methode! datasize: " << chunksize << std::endl;
    }
    std::vector<T> bVector(b.getSize());
    std::vector<T> datasetVector(dataset.getSize());
    double *data_raw = dataset.getPointer();
    for (size_t i = 0; i < dataset.getSize(); i++)
      datasetVector[i] = static_cast<T>(data_raw[i]);
    start = std::chrono::system_clock::now();
    bKernel->rhs(datasetVector, bVector, start_id, chunksize);
    end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;

    for (size_t i = 0; i < b.getSize(); i++)
      b[i] = bVector[i];
    if (verbose) {
      std::cout << "duration rhs ocl: " << elapsed_seconds.count() << std::endl;
    }
  }
};

}  // namespace DensityOCLMultiPlatform
}  // namespace datadriven
}  // namespace sgpp
