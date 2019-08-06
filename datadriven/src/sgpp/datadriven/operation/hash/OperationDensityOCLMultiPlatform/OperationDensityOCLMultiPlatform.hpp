// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org
#pragma once

#include <sgpp/base/exception/operation_exception.hpp>
#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/grid/GridStorage.hpp>
#include <sgpp/base/opencl/OCLManager.hpp>
#include <sgpp/base/opencl/OCLOperationConfiguration.hpp>
#include <sgpp/base/operation/hash/OperationMatrix.hpp>
#include <sgpp/base/tools/SGppStopwatch.hpp>
#include <sgpp/globaldef.hpp>

#include <chrono>
#include <string>
#include <vector>
#include <algorithm>

#include <sgpp/datadriven/operation/hash/OperationDensityOCLMultiPlatform/KernelB.hpp>
#include <sgpp/datadriven/operation/hash/OperationDensityOCLMultiPlatform/KernelMult.hpp>
#include <sgpp/datadriven/operation/hash/OperationDensityOCLMultiPlatform/OperationDensityOCL.hpp>

namespace sgpp {
namespace datadriven {
namespace DensityOCLMultiPlatform {

/// Class for opencl density multiplication and density right hand side vector
template <typename T>
class OperationDensityOCLMultiPlatform : public OperationDensity {
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
  OperationDensityOCLMultiPlatform(base::Grid &grid, size_t dimensions,
                                   std::shared_ptr<base::OCLManagerMultiPlatform> manager,
                                   sgpp::base::OCLOperationConfiguration *parameters, T lambda,
                                   size_t platform_id, size_t device_id)
      : OperationDensity(),
        dims(dimensions),
        gridSize(grid.getStorage().getSize()),
        devices(manager->getDevices()),
        verbose(false),
        manager(manager),
        lambda(lambda) {
    // Store Grid in a opencl compatible buffer
    sgpp::base::GridStorage &gridStorage = grid.getStorage();
    size_t pointscount = 0;
    for (size_t i = 0; i < gridSize; i++) {
      sgpp::base::HashGridPoint &point = gridStorage.getPoint(i);
      pointscount++;
      for (size_t d = 0; d < dims; d++) {
        points.push_back(point.getIndex(d));
        points.push_back(point.getLevel(d));
      }
    }

    // look for the chosen platform and device and create kernel with it
    size_t platformcounter = 0;
    size_t devicecounter = 0;
    cl_device_id old_device_id = devices[0]->deviceId;
    cl_platform_id old_platform_id = devices[0]->platformId;
    size_t counter = 0;
    bool success = false;
    for (auto device : devices) {
      if (device->platformId != old_platform_id) {
        platformcounter++;
        old_platform_id = device->platformId;
        devicecounter = 0;
      }
      if (device->deviceId != old_device_id) {
        devicecounter++;
        old_device_id = device->deviceId;
      }
      if (platformcounter == platform_id && devicecounter == device_id) {
        json::Node &deviceNode =
            (*parameters)["PLATFORMS"][device->platformName]["DEVICES"][device->deviceName];
        json::Node &firstKernelConfig = deviceNode["KERNELS"]["multdensity"];
        json::Node &secondKernelConfig = deviceNode["KERNELS"]["cscheme"];
        bKernel =
            new KernelDensityB<T>(devices[counter], dims, manager, secondKernelConfig, points);
        multKernel = new KernelDensityMult<T>(devices[counter], dims, manager, firstKernelConfig,
                                              points, lambda);
        if (firstKernelConfig["VERBOSE"].getBool()) verbose = true;
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
                                   sgpp::base::OCLOperationConfiguration *parameters, T lambda,
                                   size_t platform_id, size_t device_id)
      : OperationDensity(),
        dims(dimensions),
        gridSize(gridsize),
        devices(manager->getDevices()),
        verbose(false),
        manager(manager),
        lambda(lambda) {
    for (size_t i = 0; i < gridSize; i++) {
      for (size_t d = 0; d < dims; d++) {
        points.push_back(gridpoints[2 * dimensions * i + 2 * d]);
        points.push_back(gridpoints[2 * dimensions * i + 2 * d + 1]);
      }
    }
    // look for the chosen platform and device and create kernel with it
    size_t platformcounter = 0;
    size_t devicecounter = 0;
    cl_device_id old_device_id = devices[0]->deviceId;
    cl_platform_id old_platform_id = devices[0]->platformId;
    size_t counter = 0;
    bool success = false;
    for (auto device : devices) {
      if (device->platformId != old_platform_id) {
        platformcounter++;
        old_platform_id = device->platformId;
        devicecounter = 0;
      }
      if (device->deviceId != old_device_id) {
        devicecounter++;
        old_device_id = device->deviceId;
      }
      if (platformcounter == platform_id && devicecounter == device_id) {
        json::Node &deviceNode =
            (*parameters)["PLATFORMS"][device->platformName]["DEVICES"][device->deviceName];
        json::Node &firstKernelConfig = deviceNode["KERNELS"]["multdensity"];
        json::Node &secondKernelConfig = deviceNode["KERNELS"]["cscheme"];
        bKernel =
            new KernelDensityB<T>(devices[counter], dims, manager, secondKernelConfig, points);
        multKernel = new KernelDensityMult<T>(devices[counter], dims, manager, firstKernelConfig,
                                              points, lambda);
        if (firstKernelConfig["VERBOSE"].getBool()) verbose = true;
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

  /// Use before calling partial_mult directly
  void initialize_alpha(double *alpha) override {
    if (std::is_same<T, double>::value) {
      std::vector<T> alphaVector(alpha, alpha + gridSize);
      this->multKernel->initialize_alpha_buffer(alphaVector);
    } else {
      std::vector<T> alphaVector(gridSize);
      for (size_t i = 0; i < gridSize; ++i) {
        alphaVector[i] = static_cast<T>(alpha[i]);
      }
      this->multKernel->initialize_alpha_buffer(alphaVector);
    }
  }

  /// Execute a partial (startindex to startindex+chunksize) multiplication with the density matrix
  void start_partial_mult(int start_id, int chunksize) override {
    this->multKernel->start_mult(start_id, chunksize);
  }

  void finish_partial_mult(double *result, int start_id, int chunksize) override {
    if (std::is_same<T, double>::value) {
      std::vector<T> resultVector(result, result + chunksize);
      this->multKernel->finish_mult(resultVector, start_id, chunksize);
      std::copy(resultVector.begin(), resultVector.end(), result);
    } else {
      std::vector<T> resultVector(chunksize);
      for (int i = 0; i < chunksize; i++) {
        resultVector[i] = static_cast<T>(result[i]);
      }
      this->multKernel->finish_mult(resultVector, start_id, chunksize);
      std::copy(resultVector.begin(), resultVector.end(), result);
    }
  }

  /// Execute one matrix-vector multiplication with the density matrix
  void mult(base::DataVector &alpha, base::DataVector &result) override {
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
    this->multKernel->initialize_alpha_buffer(alphaVector);
    this->multKernel->start_mult(0, 0);
    this->multKernel->finish_mult(resultVector, 0, 0);
    end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;

    if (verbose) {
      std::cout << "duration mult ocl: " << elapsed_seconds.count() << std::endl;
    }
    for (size_t i = 0; i < gridSize; i++) result[i] = resultVector[i];
  }
  void initialize_dataset(base::DataMatrix &dataset) override {
    if (std::is_same<T, double>::value) {
      double *data_raw = dataset.getPointer();
      std::vector<T> datasetVector(data_raw, data_raw + dataset.getSize());
      bKernel->initialize_dataset(datasetVector);
    } else {
      std::vector<T> datasetVector(dataset.getSize());
      double *data_raw = dataset.getPointer();
      for (size_t i = 0; i < dataset.getSize(); i++) datasetVector[i] = static_cast<T>(data_raw[i]);
      bKernel->initialize_dataset(datasetVector);
    }
  }
  void start_rhs_generation(size_t start_id, size_t chunksize) override {
    bKernel->start_rhs_generation(start_id, chunksize);
  }

  void finalize_rhs_generation(sgpp::base::DataVector &b, size_t start_id,
                               size_t chunksize) override {
    std::vector<T> bVector(b.getSize());
    bKernel->finalize_rhs_generation(bVector, start_id, chunksize);
    for (size_t i = 0; i < b.getSize(); i++) b[i] = bVector[i];
  }

  /// Generates the right hand side vector for the density equation
  void generateb(base::DataMatrix &dataset, sgpp::base::DataVector &b, size_t start_id = 0,
                 size_t chunksize = 0) override {
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
    for (size_t i = 0; i < dataset.getSize(); i++) datasetVector[i] = static_cast<T>(data_raw[i]);
    start = std::chrono::system_clock::now();

    bKernel->initialize_dataset(datasetVector);
    bKernel->start_rhs_generation(start_id, chunksize);
    bKernel->finalize_rhs_generation(bVector, start_id, chunksize);
    end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;

    for (size_t i = 0; i < b.getSize(); i++) b[i] = bVector[i];
    if (verbose) {
      std::cout << "duration rhs ocl: " << elapsed_seconds.count() << std::endl;
    }
  }
};

}  // namespace DensityOCLMultiPlatform
}  // namespace datadriven
}  // namespace sgpp
