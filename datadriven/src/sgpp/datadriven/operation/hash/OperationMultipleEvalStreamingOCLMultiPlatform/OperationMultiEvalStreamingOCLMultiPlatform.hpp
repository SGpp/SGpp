// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <omp.h>

#include <algorithm>
#include <chrono>
#include <mutex>
#include <vector>

#include "Configuration.hpp"
#include "KernelMult.hpp"
#include "KernelMultTranspose.hpp"
#include "sgpp/base/exception/operation_exception.hpp"
#include "sgpp/base/opencl/OCLManager.hpp"
#include "sgpp/base/opencl/OCLOperationConfiguration.hpp"
#include "sgpp/base/opencl/QueueLoadBalancer.hpp"
#include "sgpp/base/operation/hash/OperationMultipleEval.hpp"
#include "sgpp/base/tools/SGppStopwatch.hpp"
#include "sgpp/globaldef.hpp"

namespace sgpp {
namespace datadriven {
namespace StreamingOCLMultiPlatform {

template <typename T>
class OperationMultiEvalStreamingOCLMultiPlatform : public base::OperationMultipleEval {
 protected:
  size_t dims;

  base::DataMatrix preparedDataset;

  std::shared_ptr<base::OCLOperationConfiguration> parameters;
  std::vector<T> kernelDataset;

  size_t datasetSizeUnpadded;
  size_t datasetSizePadded;
  size_t datasetSizeBuffers;

  //  // includes padding
  //  size_t datasetSize;

  // Member to store the sparse grid's levels for better vectorization
  std::vector<T> level;
  // Member to store the sparse grid's indices for better vectorization
  std::vector<T> index;

  size_t gridSizeUnpadded;
  size_t gridSizePadded;
  size_t gridSizeBuffers;

  //  // includes padding
  //  size_t gridSize;

  // Timer object to handle time measurements
  base::SGppStopwatch myTimer;

  double duration;

  std::shared_ptr<base::QueueLoadBalancer> queueLoadBalancerMult;
  std::shared_ptr<base::QueueLoadBalancer> queueLoadBalancerMultTranspose;

  std::shared_ptr<base::OCLManagerMultiPlatform> manager;
  std::vector<std::shared_ptr<base::OCLDevice>> devices;

  std::vector<StreamingOCLMultiPlatform::KernelMult<T>> multKernels;
  std::vector<StreamingOCLMultiPlatform::KernelMultTranspose<T>> multTransposeKernels;

  json::Node &configuration;

  bool verbose;

  size_t commonDatasetPadding;
  size_t commonGridPadding;

 public:
  OperationMultiEvalStreamingOCLMultiPlatform(
      base::Grid &grid, base::DataMatrix &dataset,
      std::shared_ptr<base::OCLManagerMultiPlatform> manager,
      std::shared_ptr<base::OCLOperationConfiguration> parameters, json::Node &configuration)
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

    // initialized in padDataset
    datasetSizeUnpadded = 0;
    datasetSizePadded = 0;
    datasetSizeBuffers = 0;

    // initialized in prepare
    gridSizeUnpadded = 0;
    gridSizePadded = 0;
    gridSizeBuffers = 0;

    this->padDataset(this->preparedDataset);
    this->preparedDataset.transpose();

    this->kernelDataset =
        std::vector<T>(this->preparedDataset.getNrows() * this->preparedDataset.getNcols());

    for (size_t i = 0; i < this->preparedDataset.getSize(); i++) {
      this->kernelDataset[i] = (T) this->preparedDataset[i];
    }

    for (size_t deviceIndex = 0; deviceIndex < devices.size(); deviceIndex++) {
      json::Node &platformConfiguration =
          (*parameters)["PLATFORMS"][devices[deviceIndex]->platformName];
      json::Node &deviceConfiguration =
          platformConfiguration["DEVICES"][devices[deviceIndex]->deviceName];
      json::Node &kernelConfiguration =
          deviceConfiguration["KERNELS"][StreamingOCLMultiPlatform::Configuration::getKernelName()];
      multKernels.emplace_back(devices[deviceIndex], dims, this->manager, kernelConfiguration,
                               queueLoadBalancerMult);
      multTransposeKernels.emplace_back(devices[deviceIndex], dims, this->manager,
                                        kernelConfiguration, queueLoadBalancerMultTranspose);
    }

    // create the kernel specific data structures and initialize gridSize and gridSizeExtra
    this->prepare();
  }

  ~OperationMultiEvalStreamingOCLMultiPlatform() {}

  void mult(base::DataVector &alpha, base::DataVector &result) override {
    this->myTimer.start();

    size_t gridFrom = 0;
    size_t gridTo = this->gridSizePadded;
    size_t datasetFrom = 0;
    size_t datasetTo = this->datasetSizePadded;

    queueLoadBalancerMult->initialize(datasetFrom, datasetTo);

    std::vector<T> alphaArray(this->gridSizePadded);

    for (size_t i = 0; i < alpha.getSize(); i++) {
      alphaArray[i] = (T)alpha[i];
    }

    for (size_t i = alpha.getSize(); i < this->gridSizePadded; i++) {
      alphaArray[i] = 0.0;
    }

    // additional padding to allow for devices with different block sizes
    std::vector<T> resultArray(this->datasetSizeBuffers);

    std::fill(resultArray.begin(), resultArray.end(), 0.0);

    omp_set_num_threads(static_cast<int>(devices.size()));

    for (size_t i = 0; i < devices.size(); i++) {
      std::cout << devices[i]->deviceName << std::endl;
    }

    std::once_flag onceFlag;
    std::exception_ptr exceptionPtr;

#pragma omp parallel
    {
      size_t threadId = omp_get_thread_num();
      try {
        this->multKernels[threadId].mult(this->level, this->index, this->kernelDataset, alphaArray,
                                         resultArray, gridFrom, gridTo, datasetFrom, datasetTo);
      } catch (...) {
        // store the first exception thrown for rethrow
        std::call_once(onceFlag, [&]() { exceptionPtr = std::current_exception(); });
      }
    }

    if (exceptionPtr) {
      std::rethrow_exception(exceptionPtr);
    }

    for (size_t i = 0; i < result.getSize(); i++) {
      result[i] = resultArray[i];
    }

    this->duration = this->myTimer.stop();

    for (StreamingOCLMultiPlatform::KernelMult<T> &kernel : multKernels) {
      this->duration -= kernel.getBuildDuration();
    }

    if (verbose) {
      std::cout << "duration mult ocl: " << this->duration << std::endl;
    }
  }

  void multTranspose(base::DataVector &source, base::DataVector &result) override {
    this->prepare();

    this->myTimer.start();

    size_t gridFrom = 0;
    size_t gridTo = this->gridSizePadded;
    size_t datasetFrom = 0;
    size_t datasetTo = this->datasetSizePadded;

    queueLoadBalancerMultTranspose->initialize(gridFrom, gridTo);

    std::vector<T> sourceArray(this->datasetSizePadded);

    for (size_t i = 0; i < source.getSize(); i++) {
      sourceArray[i] = (T)source[i];
    }

    for (size_t i = source.getSize(); i < this->datasetSizePadded; i++) {
      sourceArray[i] = 0.0;
    }

    std::vector<T> resultArray(this->gridSizeBuffers);

    std::fill(resultArray.begin(), resultArray.end(), 0.0);

    omp_set_num_threads(static_cast<int>(devices.size()));

    std::once_flag onceFlag;
    std::exception_ptr exceptionPtr;

#pragma omp parallel
    {
      size_t threadId = omp_get_thread_num();

      try {
        this->multTransposeKernels[threadId].multTranspose(
            this->level, this->index, this->kernelDataset, sourceArray, resultArray, gridFrom,
            gridTo, datasetFrom, datasetTo);
      } catch (...) {
        // store the first exception thrown for rethrow
        std::call_once(onceFlag, [&]() { exceptionPtr = std::current_exception(); });
      }
    }

    if (exceptionPtr) {
      std::rethrow_exception(exceptionPtr);
    }

    for (size_t i = 0; i < result.getSize(); i++) {
      result[i] = resultArray[i];
    }

    this->duration = this->myTimer.stop();

    for (StreamingOCLMultiPlatform::KernelMultTranspose<T> &kernelTranspose :
         multTransposeKernels) {
      this->duration -= kernelTranspose.getBuildDuration();
    }

    if (verbose) {
      std::cout << "duration multTranspose ocl: " << duration << std::endl;
    }
  }

  double getDuration() override { return this->duration; }

  void prepare() override {
    this->recalculateLevelAndIndex();

    for (size_t deviceIndex = 0; deviceIndex < devices.size(); deviceIndex++) {
      this->multKernels[deviceIndex].resetKernel();
      this->multTransposeKernels[deviceIndex].resetKernel();
    }
  }

 private:
  void padDataset(base::DataMatrix &dataset) {
    // Assure that data has a even number of instances -> padding might be
    // needed
    size_t remainder = dataset.getNrows() % commonDatasetPadding;
    size_t padding = commonDatasetPadding - remainder;
    datasetSizeUnpadded = dataset.getNrows();
    datasetSizePadded = dataset.getNrows() + padding;
    datasetSizeBuffers = dataset.getNrows() + commonDatasetPadding;

    // replicate last row for padding
    base::DataVector lastRow(dims);
    dataset.getRow(datasetSizeUnpadded - 1, lastRow);
    dataset.resize(datasetSizeBuffers);

    for (size_t i = datasetSizeUnpadded; i < datasetSizeBuffers; i++) {
      dataset.setRow(i, lastRow);
    }
  }

  void recalculateLevelAndIndex() {
    base::GridStorage &storage = grid.getStorage();

    size_t remainder = storage.getSize() % commonGridPadding;
    size_t padding = 0;

    if (remainder != 0) {
      padding = commonGridPadding - remainder;
    }

    gridSizeUnpadded = storage.getSize();
    gridSizePadded = storage.getSize() + padding;
    gridSizeBuffers = storage.getSize() + commonGridPadding;

    level = std::vector<T>(gridSizeBuffers * dims);
    index = std::vector<T>(gridSizeBuffers * dims);

    base::HashGridPoint::level_type curLevel;
    base::HashGridPoint::index_type curIndex;

    for (size_t i = 0; i < storage.getSize(); i++) {
      base::HashGridPoint &gridPoint = storage.getPoint(i);
      for (size_t dim = 0; dim < dims; dim++) {
        gridPoint.get(dim, curLevel, curIndex);
        level[i * dims + dim] = static_cast<T>(1 << curLevel);
        index[i * dims + dim] = static_cast<T>(curIndex);
      }
    }

    for (size_t i = storage.getSize(); i < gridSizeBuffers; i++) {
      for (size_t dim = 0; dim < storage.getDimension(); dim++) {
        level[i * dims + dim] = 1.0;
        index[i * dims + dim] = 1.0;
      }
    }
  }

  size_t calculateCommonDatasetPadding() {
    size_t commonPaddingRequiredment = 1;
    for (size_t deviceIndex = 0; deviceIndex < devices.size(); deviceIndex++) {
      json::Node &platformConfiguration =
          (*parameters)["PLATFORMS"][devices[deviceIndex]->platformName];
      json::Node &deviceConfiguration =
          platformConfiguration["DEVICES"][devices[deviceIndex]->deviceName];
      json::Node &kernelConfiguration =
          deviceConfiguration["KERNELS"][StreamingOCLMultiPlatform::Configuration::getKernelName()];

      commonPaddingRequiredment = std::max(commonPaddingRequiredment,
                                           kernelConfiguration["KERNEL_DATA_BLOCK_SIZE"].getUInt() *
                                               kernelConfiguration["LOCAL_SIZE"].getUInt());
    }
    return commonPaddingRequiredment;
  }

  size_t calculateCommonGridPadding() {
    size_t commonPaddingRequiredment = 1;
    for (size_t deviceIndex = 0; deviceIndex < devices.size(); deviceIndex++) {
      json::Node &platformConfiguration =
          (*parameters)["PLATFORMS"][devices[deviceIndex]->platformName];
      json::Node &deviceConfiguration =
          platformConfiguration["DEVICES"][devices[deviceIndex]->deviceName];
      json::Node &kernelConfiguration =
          deviceConfiguration["KERNELS"][StreamingOCLMultiPlatform::Configuration::getKernelName()];

      commonPaddingRequiredment = std::max(
          commonPaddingRequiredment, kernelConfiguration["KERNEL_TRANS_GRID_BLOCK_SIZE"].getUInt() *
                                         kernelConfiguration["LOCAL_SIZE"].getUInt());
    }
    return commonPaddingRequiredment;
  }
};

}  // namespace StreamingOCLMultiPlatform
}  // namespace datadriven
}  // namespace sgpp
