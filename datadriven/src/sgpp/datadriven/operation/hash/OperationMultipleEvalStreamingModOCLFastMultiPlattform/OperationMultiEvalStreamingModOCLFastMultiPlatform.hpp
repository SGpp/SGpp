// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <omp.h>

#include <chrono>
#include <vector>
#include <algorithm>

#include "sgpp/base/operation/hash/OperationMultipleEval.hpp"
#include "sgpp/base/tools/SGppStopwatch.hpp"
#include "sgpp/base/exception/operation_exception.hpp"
#include "sgpp/globaldef.hpp"
#include "sgpp/base/opencl/OCLOperationConfiguration.hpp"
#include "sgpp/base/opencl/OCLManager.hpp"
#include "sgpp/base/opencl/QueueLoadBalancer.hpp"
#include "KernelMult.hpp"
#include "KernelMultTranspose.hpp"
#include "Configuration.hpp"

namespace SGPP {
namespace datadriven {

template <typename T>
class OperationMultiEvalStreamingModOCLFastMultiPlatform : public base::OperationMultipleEval {
 protected:
  size_t dims;

  SGPP::base::DataMatrix preparedDataset;

  std::shared_ptr<base::OCLOperationConfiguration> parameters;

  std::vector<T> kernelDataset;
  size_t datasetSize = 0;
  /// Member to store the sparse grid's levels for better vectorization
  std::vector<T> level;
  /// Member to store the sparse grid's indices for better vectorization
  std::vector<T> index;

  size_t gridSize = 0;

  /// Timer object to handle time measurements
  SGPP::base::SGppStopwatch myTimer;

  float_t duration;

  std::shared_ptr<base::QueueLoadBalancer> queueLoadBalancerMult;
  std::shared_ptr<base::QueueLoadBalancer> queueLoadBalancerMultTranspose;

  std::shared_ptr<base::OCLManagerMultiPlatform> manager;
  std::vector<std::shared_ptr<base::OCLDevice>> devices;

  std::vector<StreamingModOCLFastMultiPlatform::KernelMult<T>> multKernels;
  std::vector<StreamingModOCLFastMultiPlatform::KernelMultTranspose<T>> multTransposeKernels;

  json::Node& configuration;

  bool verbose;

  size_t commonDatasetPadding;
  size_t commonGridPadding;

 public:
  OperationMultiEvalStreamingModOCLFastMultiPlatform(
      base::Grid& grid, base::DataMatrix& dataset,
      std::shared_ptr<base::OCLManagerMultiPlatform> manager,
      std::shared_ptr<base::OCLOperationConfiguration> parameters)
      : OperationMultipleEval(grid, dataset),
        preparedDataset(dataset),
        parameters(parameters),
        myTimer(base::SGppStopwatch()),
        duration(-1.0),
        manager(manager),
        devices(manager->getDevices()),
        configuration(configuration) {
    this->dims = dataset.getNcols();  // be aware of transpose!
    this->verbose = configuration["VERBOSE"].getBool();

    this->commonDatasetPadding = calculateCommonDatasetPadding();
    this->commonGridPadding = calculateCommonGridPadding();

    queueLoadBalancerMult = std::make_shared<base::QueueLoadBalancer>();
    queueLoadBalancerMultTranspose = std::make_shared<base::QueueLoadBalancer>();

    this->padDataset(this->preparedDataset, datasetSize);
    this->preparedDataset.transpose();

    //    std::cout << "dims: " << this->dims << std::endl;
    //    std::cout << "padded instances: " << this->datasetSize << std::endl;

    this->kernelDataset =
        std::vector<T>(this->preparedDataset.getNrows() * this->preparedDataset.getNcols());

    for (size_t i = 0; i < this->preparedDataset.getSize(); i++) {
      this->kernelDataset[i] = (T) this->preparedDataset[i];
    }

    for (size_t deviceIndex = 0; deviceIndex < devices.size(); deviceIndex++) {
      json::Node& platformConfiguration =
          (*parameters)["PLATFORMS"][devices[deviceIndex]->platformName];
      json::Node& deviceConfiguration =
          platformConfiguration["DEVICES"][devices[deviceIndex]->deviceName];
      json::Node& kernelConfiguration = deviceConfiguration
          ["KERNELS"][StreamingModOCLFastMultiPlatform::Configuration::getKernelName()];

      multKernels.emplace_back(devices[deviceIndex], dims, this->manager, kernelConfiguration,
                               queueLoadBalancerMult);
      multTransposeKernels.emplace_back(devices[deviceIndex], dims, this->manager,
                                        kernelConfiguration, queueLoadBalancerMultTranspose);
    }

    // create the kernel specific data structures and initialize gridSize and gridSizeExtra
    this->prepare();
  }

  ~OperationMultiEvalStreamingModOCLFastMultiPlatform() {}

  void mult(base::DataVector& alpha, base::DataVector& result) override {
    this->myTimer.start();

    size_t gridFrom = 0;
    size_t gridTo = this->gridSize;
    size_t datasetFrom = 0;
    size_t datasetTo = this->datasetSize;

    queueLoadBalancerMult->initialize(datasetFrom, datasetTo);

    std::vector<T> alphaArray(this->gridSize);

    for (size_t i = 0; i < alpha.getSize(); i++) {
      alphaArray[i] = (T)alpha[i];
    }

    for (size_t i = alpha.getSize(); i < this->gridSize; i++) {
      alphaArray[i] = 0.0;
    }

    // additional padding to allow for devices with different block sizes
    std::vector<T> resultArray(this->datasetSize);

    for (size_t i = 0; i < this->datasetSize; i++) {
      resultArray[i] = 0.0;
    }

    std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();

    omp_set_num_threads(static_cast<int>(devices.size()));

    for (size_t i = 0; i < devices.size(); i++) {
      std::cout << devices[i]->deviceName << std::endl;
    }

#pragma omp parallel
    {
      size_t threadId = omp_get_thread_num();
      this->multKernels[threadId].mult(this->level, this->index, this->gridSize,
                                       this->kernelDataset, this->datasetSize, alphaArray,
                                       resultArray, gridFrom, gridTo, datasetFrom, datasetTo);
    }
    end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;

    if (verbose) {
      std::cout << "duration mult ocl: " << elapsed_seconds.count() << std::endl;
    }

    for (size_t i = 0; i < result.getSize(); i++) {
      result[i] = resultArray[i];
    }

    this->duration = this->myTimer.stop();
  }

  void multTranspose(base::DataVector& source, base::DataVector& result) override {
    this->prepare();

    this->myTimer.start();

    size_t gridFrom = 0;
    size_t gridTo = this->gridSize;
    size_t datasetFrom = 0;
    size_t datasetTo = this->datasetSize;

    queueLoadBalancerMultTranspose->initialize(gridFrom, gridTo);

    std::vector<T> sourceArray(this->datasetSize);

    for (size_t i = 0; i < source.getSize(); i++) {
      sourceArray[i] = (T)source[i];
    }

    for (size_t i = source.getSize(); i < this->datasetSize; i++) {
      sourceArray[i] = 0.0;
    }

    std::vector<T> resultArray(this->gridSize);

    std::fill(resultArray.begin(), resultArray.end(), 0.0);

    std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();

    omp_set_num_threads(static_cast<int>(devices.size()));

#pragma omp parallel
    {
      size_t threadId = omp_get_thread_num();

      this->multTransposeKernels[threadId].multTranspose(
          this->level, this->index, this->kernelDataset, sourceArray, resultArray, gridFrom, gridTo,
          datasetFrom, datasetTo);
    }
    end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;
    if (verbose) {
      std::cout << "duration multTranspose ocl: " << elapsed_seconds.count() << std::endl;
    }

    for (size_t i = 0; i < result.getSize(); i++) {
      result[i] = resultArray[i];
    }

    this->duration = this->myTimer.stop();
  }

  float_t getDuration() override { return this->duration; }

  void prepare() override {
    this->recalculateLevelAndIndex(gridSize);

    for (size_t deviceIndex = 0; deviceIndex < devices.size(); deviceIndex++) {
      this->multKernels[deviceIndex].resetKernel();
      this->multTransposeKernels[deviceIndex].resetKernel();
    }
  }

 private:
  void padDataset(base::DataMatrix& dataset, size_t& datasetSize) {
    // Assure that data has a even number of instances -> padding might be
    // needed
    size_t remainder = dataset.getNrows() % commonDatasetPadding;
    size_t padding = commonDatasetPadding - remainder;

    datasetSize = dataset.getNrows() + padding;

    if (padding != commonDatasetPadding) {
      base::DataVector lastRow(dataset.getNcols());
      size_t oldSize = dataset.getNrows();
      dataset.getRow(dataset.getNrows() - 1, lastRow);
      dataset.resize(dataset.getNrows() + padding);

      for (size_t i = 0; i < padding; i++) {
        dataset.setRow(oldSize + i, lastRow);
      }
    }
  }

  void recalculateLevelAndIndex(size_t& gridSize) {
    base::GridStorage* storage = grid.getStorage();

    size_t remainder = storage->size() % commonGridPadding;
    size_t padding = 0;

    if (remainder != 0) {
      padding = commonGridPadding - remainder;
    }

    gridSize = storage->size() + padding;

    level = std::vector<T>(gridSize * dims);
    index = std::vector<T>(gridSize * dims);

    base::HashGridIndex::level_type curLevel;
    base::HashGridIndex::index_type curIndex;

    /// pointer to index_type
    base::HashGridStorage::index_pointer gridPoint;

    for (size_t i = 0; i < storage->size(); i++) {
      gridPoint = storage->get(i);
      for (size_t dim = 0; dim < dims; dim++) {
        gridPoint->get(dim, curLevel, curIndex);
        level[i * dims + dim] = static_cast<T>(1 << curLevel);
        index[i * dims + dim] = static_cast<T>(curIndex);
      }
    }

    for (size_t i = storage->size(); i < this->gridSize; i++) {
      for (size_t dim = 0; dim < storage->dim(); dim++) {
        level[i * dims + dim] = 1.0;
        index[i * dims + dim] = 1.0;
      }
    }
  }

  size_t calculateCommonDatasetPadding() {
    size_t commonPaddingRequiredment = 1;
    for (size_t deviceIndex = 0; deviceIndex < devices.size(); deviceIndex++) {
      json::Node& platformConfiguration =
          (*parameters)["PLATFORMS"][devices[deviceIndex]->platformName];
      json::Node& deviceConfiguration =
          platformConfiguration["DEVICES"][devices[deviceIndex]->deviceName];
      json::Node& kernelConfiguration = deviceConfiguration
          ["KERNELS"][StreamingModOCLFastMultiPlatform::Configuration::getKernelName()];

      commonPaddingRequiredment = std::max(commonPaddingRequiredment,
                                           kernelConfiguration["KERNEL_DATA_BLOCK_SIZE"].getUInt() *
                                               kernelConfiguration["LOCAL_SIZE"].getUInt());
    }
    return commonPaddingRequiredment;
  }

  size_t calculateCommonGridPadding() {
    size_t commonPaddingRequiredment = 1;
    for (size_t deviceIndex = 0; deviceIndex < devices.size(); deviceIndex++) {
      json::Node& platformConfiguration =
          (*parameters)["PLATFORMS"][devices[deviceIndex]->platformName];
      json::Node& deviceConfiguration =
          platformConfiguration["DEVICES"][devices[deviceIndex]->deviceName];
      json::Node& kernelConfiguration = deviceConfiguration
          ["KERNELS"][StreamingModOCLFastMultiPlatform::Configuration::getKernelName()];

      commonPaddingRequiredment = std::max(
          commonPaddingRequiredment, kernelConfiguration["KERNEL_TRANS_GRID_BLOCK_SIZE"].getUInt() *
                                         kernelConfiguration["LOCAL_SIZE"].getUInt());
    }
    return commonPaddingRequiredment;
  }
};
}  // namespace datadriven
}  // namespace SGPP
