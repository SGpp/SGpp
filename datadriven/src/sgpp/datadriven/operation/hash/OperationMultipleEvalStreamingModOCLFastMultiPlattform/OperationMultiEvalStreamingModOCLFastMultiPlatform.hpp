// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <omp.h>

#include <algorithm>
#include <chrono>
#include <vector>

#include <sgpp/datadriven/operation/hash/OperationMultipleEvalStreamingModOCLFastMultiPlattform/Configuration.hpp>
#include <sgpp/datadriven/operation/hash/OperationMultipleEvalStreamingModOCLFastMultiPlattform/KernelMult.hpp>
#include <sgpp/datadriven/operation/hash/OperationMultipleEvalStreamingModOCLFastMultiPlattform/KernelMultTranspose.hpp>
#include <sgpp/base/exception/operation_exception.hpp>
#include <sgpp/base/opencl/OCLManager.hpp>
#include <sgpp/base/opencl/OCLOperationConfiguration.hpp>
#include <sgpp/base/tools/QueueLoadBalancerOpenMP.hpp>
#include <sgpp/base/operation/hash/OperationMultipleEval.hpp>
#include <sgpp/base/tools/SGppStopwatch.hpp>
#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace datadriven {

template <typename T>
class OperationMultiEvalStreamingModOCLFastMultiPlatform : public base::OperationMultipleEval {
 protected:
  size_t dims;

  sgpp::base::DataMatrix preparedDataset;

  std::shared_ptr<base::OCLOperationConfiguration> parameters;

  std::vector<T> kernelDataset;

  size_t datasetSizeUnpadded;
  size_t datasetSizePadded;

  //  size_t datasetSize = 0;
  /// Member to store the sparse grid's levels for better vectorization
  std::vector<T> level;
  /// Member to store the sparse grid's indices for better vectorization
  std::vector<T> index;

  size_t gridSizeUnpadded;
  size_t gridSizePadded;

  //  size_t gridSize = 0;

  /// Timer object to handle time measurements
  sgpp::base::SGppStopwatch myTimer;

  double duration;

  std::shared_ptr<base::QueueLoadBalancerMutex> queueLoadBalancerMult;
  std::shared_ptr<base::QueueLoadBalancerMutex> queueLoadBalancerMultTranspose;

  std::shared_ptr<base::OCLManagerMultiPlatform> manager;
  std::vector<std::shared_ptr<base::OCLDevice>> devices;

  std::vector<StreamingModOCLFastMultiPlatform::KernelMult<T>> multKernels;
  std::vector<StreamingModOCLFastMultiPlatform::KernelMultTranspose<T>> multTransposeKernels;

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
        devices(manager->getDevices()) {
    this->dims = dataset.getNcols();  // be aware of transpose!
    this->verbose = (*parameters)["VERBOSE"].getBool();

    this->commonDatasetPadding = calculateCommonDatasetPadding();
    this->commonGridPadding = calculateCommonGridPadding();

    queueLoadBalancerMult = std::make_shared<base::QueueLoadBalancerMutex>();
    queueLoadBalancerMultTranspose = std::make_shared<base::QueueLoadBalancerMutex>();

    // initialized in pad
    datasetSizeUnpadded = 0;
    datasetSizePadded = 0;

    // initialized in prepare
    gridSizeUnpadded = 0;
    gridSizePadded = 0;

    this->padDataset(this->preparedDataset);
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
      json::Node& kernelConfiguration =
          deviceConfiguration["KERNELS"]
                             [StreamingModOCLFastMultiPlatform::Configuration::getKernelName()];

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
    size_t gridTo = this->gridSizePadded;
    size_t datasetFrom = 0;
    size_t datasetTo = this->datasetSizePadded;

    queueLoadBalancerMult->initialize(datasetFrom, datasetTo, commonDatasetPadding);

    std::vector<T> alphaArray(this->gridSizePadded);

    for (size_t i = 0; i < alpha.getSize(); i++) {
      alphaArray[i] = (T)alpha[i];
    }

    for (size_t i = alpha.getSize(); i < this->gridSizePadded; i++) {
      alphaArray[i] = 0.0;
    }

    // additional padding to allow for devices with different block sizes
    std::vector<T> resultArray(this->datasetSizePadded);

    for (size_t i = 0; i < this->datasetSizePadded; i++) {
      resultArray[i] = 0.0;
    }

    std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();

    int oldThreads = omp_get_max_threads();
    omp_set_num_threads(static_cast<int>(devices.size()));

    for (size_t i = 0; i < devices.size(); i++) {
      std::cout << devices[i]->deviceName << std::endl;
    }

#pragma omp parallel
    {
      size_t threadId = omp_get_thread_num();
      this->multKernels[threadId].mult(this->level, this->index, this->kernelDataset, alphaArray,
                                       resultArray, gridFrom, gridTo, datasetFrom, datasetTo);
    }
    end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;

    if (verbose) {
      std::cout << "duration mult ocl: " << elapsed_seconds.count() << std::endl;
    }

    result.resize(this->datasetSizePadded);
    for (size_t i = 0; i < result.getSize(); i++) {
      result[i] = resultArray[i];
    }

    // restore old value of OMP_NUM_THREADS
    omp_set_num_threads(oldThreads);

    this->duration = this->myTimer.stop();
  }

  void multTranspose(base::DataVector& source, base::DataVector& result) override {
    this->prepare();

    this->myTimer.start();

    size_t gridFrom = 0;
    size_t gridTo = this->gridSizePadded;
    size_t datasetFrom = 0;
    size_t datasetTo = this->datasetSizePadded;

    queueLoadBalancerMultTranspose->initialize(gridFrom, gridTo, commonGridPadding);

    std::vector<T> sourceArray(this->datasetSizePadded);

    for (size_t i = 0; i < source.getSize(); i++) {
      sourceArray[i] = (T)source[i];
    }

    for (size_t i = source.getSize(); i < this->datasetSizePadded; i++) {
      sourceArray[i] = 0.0;
    }

    std::vector<T> resultArray(this->gridSizePadded);

    std::fill(resultArray.begin(), resultArray.end(), 0.0);

    std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();

    int oldThreads = omp_get_max_threads();
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

    result.resize(this->gridSizePadded);
    for (size_t i = 0; i < result.getSize(); i++) {
      result[i] = resultArray[i];
    }

    // restore old value of OMP_NUM_THREADS
    omp_set_num_threads(oldThreads);

    this->duration = this->myTimer.stop();
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
  void padDataset(base::DataMatrix& dataset) {
    datasetSizeUnpadded = dataset.getNrows();
    // Assure that data has a even number of instances -> padding might be
    // needed
    size_t remainder = datasetSizeUnpadded % commonDatasetPadding;
    size_t padding = commonDatasetPadding - remainder;

    datasetSizePadded = dataset.getNrows() + padding;

    if (padding != commonDatasetPadding) {
      base::DataVector lastRow(dims);
      dataset.getRow(datasetSizeUnpadded - 1, lastRow);
      dataset.resize(datasetSizePadded);

      for (size_t i = datasetSizeUnpadded; i < datasetSizePadded; i++) {
        dataset.setRow(i, lastRow);
      }
    }
  }

  void recalculateLevelAndIndex() {
    base::GridStorage& storage = grid.getStorage();

    size_t remainder = storage.getSize() % commonGridPadding;
    size_t padding = 0;

    if (remainder != 0) {
      padding = commonGridPadding - remainder;
    }

    gridSizeUnpadded = storage.getSize();
    gridSizePadded = storage.getSize() + padding;

    level = std::vector<T>(gridSizePadded * dims);
    index = std::vector<T>(gridSizePadded * dims);

    base::HashGridPoint::level_type curLevel;
    base::HashGridPoint::index_type curIndex;

    for (size_t i = 0; i < gridSizeUnpadded; i++) {
      base::HashGridPoint& gridPoint = storage.getPoint(i);
      for (size_t dim = 0; dim < dims; dim++) {
        gridPoint.get(dim, curLevel, curIndex);
        level[i * dims + dim] = static_cast<T>(1 << curLevel);
        index[i * dims + dim] = static_cast<T>(curIndex);
      }
    }

    for (size_t i = gridSizeUnpadded; i < this->gridSizePadded; i++) {
      for (size_t dim = 0; dim < storage.getDimension(); dim++) {
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
      json::Node& kernelConfiguration =
          deviceConfiguration["KERNELS"]
                             [StreamingModOCLFastMultiPlatform::Configuration::getKernelName()];

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
      json::Node& kernelConfiguration =
          deviceConfiguration["KERNELS"]
                             [StreamingModOCLFastMultiPlatform::Configuration::getKernelName()];

      commonPaddingRequiredment = std::max(
          commonPaddingRequiredment, kernelConfiguration["KERNEL_TRANS_GRID_BLOCK_SIZE"].getUInt() *
                                         kernelConfiguration["LOCAL_SIZE"].getUInt());
    }
    return commonPaddingRequiredment;
  }
};
}  // namespace datadriven
}  // namespace sgpp
