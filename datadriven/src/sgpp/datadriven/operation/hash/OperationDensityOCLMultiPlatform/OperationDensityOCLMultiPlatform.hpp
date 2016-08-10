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

#include <chrono>
#include <vector>
#include <string>

#include "OperationDensityOCL.hpp"
#include "KernelMult.hpp"
#include "KernelB.hpp"

namespace sgpp {
namespace datadriven {
namespace DensityOCLMultiPlatform {

/// Class for opencl density multiplication and density right hand side vector
template<typename T>
class OperationDensityOCLMultiPlatform: public OperationDensityOCL {
 private:
  size_t dims;
  size_t gridSize;
  /// OpenCL kernel which executes the matrix-vector density multiplications
  sgpp::datadriven::DensityOCLMultiPlatform::KernelDensityMult<T> *multKernel;
  /// OpenCL kernel which generates the right hand side vector of the density equation
  sgpp::datadriven::DensityOCLMultiPlatform::KernelDensityB<T> *bKernel;
  /// Vector with all OpenCL devices
  std::vector<std::shared_ptr<base::OCLDevice>> devices;
  /// Verbosity
  bool verbose;
  /// Contains all levels and indices of the sparse grid
  std::vector<int> points;
  /// OpenCL Manager
  std::shared_ptr<base::OCLManagerMultiPlatform> manager;
  /// Lambda for the density multiplication
  T lambda;

 public:
  /// Normal constructor
  OperationDensityOCLMultiPlatform(base::Grid& grid, size_t dimensions,
                                   std::shared_ptr<base::OCLManagerMultiPlatform> manager,
                                   sgpp::base::OCLOperationConfiguration *parameters,
                                   T lambda, size_t platform_id, size_t device_id)
      : OperationDensityOCL(), dims(dimensions),
        gridSize(grid.getStorage().getSize()),
        devices(manager->getDevices()),
        verbose(false),
        manager(manager), lambda(lambda) {
    // Store Grid in a opencl compatible buffer
    sgpp::base::GridStorage& gridStorage = grid.getStorage();
    size_t pointscount = 0;
    for (int i = 0; i < gridSize; i++) {
      sgpp::base::HashGridPoint &point = gridStorage.getPoint(i);
      pointscount++;
      for (int d = 0; d < dims; d++) {
        points.push_back(point.getIndex(d));
        points.push_back(point.getLevel(d));
      }
    }

    // look for the chosen platform and device and create kernel with it
    size_t platformcounter = 0;
    size_t devicecounter = 0;
    std::string currentplatformName = devices[0]->platformName;
    std::string currentdeviceName = devices[0]->deviceName;
    size_t counter = 0;
    bool success = false;
    for (auto device : devices) {
      if (devices[counter]->platformName != currentplatformName) {
        platformcounter++;
        currentplatformName = devices[counter]->platformName;
        devicecounter = 0;
        currentdeviceName = devices[counter]->deviceName;
      } else if (devices[counter]->deviceName != currentdeviceName) {
        devicecounter++;
        currentdeviceName = devices[counter]->deviceName;
      }
      if (platformcounter == platform_id &&
          devicecounter == device_id) {
        json::Node &deviceNode =
            (*parameters)["PLATFORMS"][currentplatformName]["DEVICES"][currentdeviceName];
        json::Node &firstKernelConfig = deviceNode["KERNELS"]["multdensity"];
        json::Node &secondKernelConfig = deviceNode["KERNELS"]["cscheme"];
        bKernel = new KernelDensityB<T>(devices[counter], dims, manager, secondKernelConfig,
                                        points);
        multKernel = new KernelDensityMult<T>(devices[counter], dims, manager, firstKernelConfig,
                                              points, lambda);
        if (firstKernelConfig["VERBOSE"].getBool())
          verbose = true;
        success = true;
        break;
      }
      counter++;
    }
    // Check whether a kernel was created or not
    if (!success) {
      std::stringstream errorString;
      errorString << "OCL Error: Platform with index " << platform_id
                  << " and the device with index " << device_id << std::endl
                  << " not found! Please check your OpenCL installation!" << std::endl;
      throw base::operation_exception(errorString.str());
    }
  }
  /// Constructor for mpi nodes - accepts grid als integer array
  OperationDensityOCLMultiPlatform(int *gridpoints, size_t gridsize, size_t dimensions,
                                   std::shared_ptr<base::OCLManagerMultiPlatform> manager,
                                   sgpp::base::OCLOperationConfiguration *parameters,
                                   T lambda, size_t platform_id, size_t device_id) :
      OperationDensityOCL(), dims(dimensions),
      gridSize(gridsize),
      devices(manager->getDevices()),
      verbose(false),
      manager(manager), lambda(lambda) {
    for (int i = 0; i < gridSize; i++) {
      for (int d = 0; d < dims; d++) {
        points.push_back(gridpoints[2 * dimensions * i + 2 * d]);
        points.push_back(gridpoints[2 * dimensions * i + 2 * d + 1]);
      }
    }
    // look for the chosen platform and device and create kernel with it
    size_t platformcounter = 0;
    size_t devicecounter = 0;
    std::string currentplatformName = devices[0]->platformName;
    std::string currentdeviceName = devices[0]->deviceName;
    size_t counter = 0;
    bool success = false;
    for (auto device : devices) {
      if (devices[counter]->platformName != currentplatformName) {
        platformcounter++;
        currentplatformName = devices[counter]->platformName;
        devicecounter = 0;
        currentdeviceName = devices[counter]->deviceName;
      } else if (devices[counter]->deviceName != currentdeviceName) {
        devicecounter++;
        currentdeviceName = devices[counter]->deviceName;
      }
      if (platformcounter == platform_id &&
          devicecounter == device_id) {
        json::Node &deviceNode =
            (*parameters)["PLATFORMS"][currentplatformName]["DEVICES"][currentdeviceName];
        json::Node &firstKernelConfig = deviceNode["KERNELS"]["multdensity"];
        json::Node &secondKernelConfig = deviceNode["KERNELS"]["cscheme"];
        bKernel = new KernelDensityB<T>(devices[counter], dims, manager, secondKernelConfig,
                                        points);
        multKernel = new KernelDensityMult<T>(devices[counter], dims, manager, firstKernelConfig,
                                              points, lambda);
        if (firstKernelConfig["VERBOSE"].getBool())
          verbose = true;
        success = true;
        break;
      }
      counter++;
    }
    // Check whether a kernel was created or not
    if (!success) {
      std::stringstream errorString;
      errorString << "OCL Error: Platform with index " << platform_id
                  << " and the device with index " << device_id << std::endl
                  << " not found! Please check your OpenCL installation!" << std::endl;
      throw base::operation_exception(errorString.str());
    }
  }

  ~OperationDensityOCLMultiPlatform() {
    delete multKernel;
    delete bKernel;
  }

  /// Execute a partial (startindex to startindex+chunksize) multiplication with the density matrix
  void partial_mult(double *alpha, double *result, size_t start_id, size_t chunksize) override {
    std::vector<T> alphaVector(gridSize);
    std::vector<T> resultVector(chunksize);
    for (size_t i = 0; i < chunksize; i++) {
      resultVector[i] = static_cast<T>(result[i]);
    }
    for (auto i = 0; i < gridSize; ++i) {
      alphaVector[i] = static_cast<T>(alpha[i]);
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

  /// Execute one matrix-vector multiplication with the density matrix
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
    this->multKernel->mult(alphaVector, resultVector, 0, 0);
    end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;

    if (verbose) {
      std::cout << "duration mult ocl: " << elapsed_seconds.count() << std::endl;
    }
    for (size_t i = 0; i < gridSize; i++)
      result[i] = resultVector[i];
  }

  /// Generates the right hand side vector for the density equation
  void generateb(base::DataMatrix &dataset, sgpp::base::DataVector &b,
                 size_t start_id = 0, size_t chunksize = 0) override {
    std::chrono::time_point<std::chrono::system_clock> start, end;
    if (verbose) {
      if (chunksize == 0)
        std::cout << "starting rhs kernel methode! datasize: " << b.getSize() << std::endl;
      else
        std::cout << "starting rhs kernel methode! chunksize: " << chunksize << std::endl;
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
