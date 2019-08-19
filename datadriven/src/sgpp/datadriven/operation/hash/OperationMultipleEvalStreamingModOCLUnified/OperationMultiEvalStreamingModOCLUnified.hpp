// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <omp.h>

#include <algorithm>
#include <chrono>
#include <vector>

#include <sgpp/datadriven/operation/hash/OperationMultipleEvalStreamingModOCLUnified/Configuration.hpp>
#include <sgpp/datadriven/operation/hash/OperationMultipleEvalStreamingModOCLUnified/KernelMult.hpp>
#include <sgpp/datadriven/operation/hash/OperationMultipleEvalStreamingModOCLUnified/KernelMultTranspose.hpp>
#include <sgpp/base/exception/operation_exception.hpp>
#include <sgpp/base/opencl/OCLManagerMultiPlatform.hpp>
#include <sgpp/base/opencl/OCLOperationConfiguration.hpp>
#include <sgpp/base/operation/hash/OperationMultipleEval.hpp>
#include <sgpp/base/tools/QueueLoadBalancerOpenMP.hpp>
#include <sgpp/base/tools/SGppStopwatch.hpp>
#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace datadriven {

template <typename T>
class OperationMultiEvalStreamingModOCLUnified
    : public base::OperationMultipleEval {
protected:
  bool verbose;
  size_t dims;
  sgpp::base::DataMatrix preparedDataset;
  std::shared_ptr<base::OCLOperationConfiguration> parameters;
  std::vector<T> kernelDataset;
  size_t datasetSizeUnpadded;
  size_t datasetSizePadded;

  std::vector<T> level;
  std::vector<T> index;
  std::vector<T> scaling;

  size_t gridSizeUnpadded;
  size_t gridSizePadded;

  /// Timer object to handle time measurements
  sgpp::base::SGppStopwatch myTimer;

  base::GridStorage &storage;

  double duration;

  std::shared_ptr<base::OCLManagerMultiPlatform> manager;
  std::vector<std::shared_ptr<base::OCLDevice>> devices;

  std::vector<StreamingModOCLUnified::KernelMult<T>> multKernels;
  std::vector<StreamingModOCLUnified::KernelMultTranspose<T>>
      multTransposeKernels;

  std::shared_ptr<sgpp::base::QueueLoadBalancerOpenMP> queueLoadBalancerMult;
  std::shared_ptr<sgpp::base::QueueLoadBalancerOpenMP>
      queueLoadBalancerMultTrans;

  size_t overallGridBlockingSize;
  size_t overallDataBlockingSize;

public:
  OperationMultiEvalStreamingModOCLUnified(
      base::Grid &grid, base::DataMatrix &dataset,
      std::shared_ptr<base::OCLManagerMultiPlatform> manager,
      std::shared_ptr<base::OCLOperationConfiguration> parameters)
      : OperationMultipleEval(grid, dataset), preparedDataset(dataset),
        parameters(parameters), myTimer(sgpp::base::SGppStopwatch()),
        storage(grid.getStorage()), duration(-1.0), manager(manager),
        devices(manager->getDevices()) {
    this->verbose = (*parameters)["VERBOSE"].getBool();

    this->dims = dataset.getNcols(); // be aware of transpose!

    overallGridBlockingSize = calculateCommonGridPadding();
    overallDataBlockingSize = calculateCommonDatasetPadding();

    // initialized in prepare
    this->gridSizeUnpadded = 0;
    this->gridSizePadded = 0;

    // initialize in pad
    this->datasetSizeUnpadded = 0;
    this->datasetSizePadded = 0;
    this->padDataset(this->preparedDataset);
    this->preparedDataset.transpose();

    queueLoadBalancerMult =
        std::make_shared<sgpp::base::QueueLoadBalancerOpenMP>();
    queueLoadBalancerMultTrans =
        std::make_shared<sgpp::base::QueueLoadBalancerOpenMP>();

    //    std::cout << "dims: " << this->dims << std::endl;
    //    std::cout << "padded instances: " << this->datasetSize << std::endl;

    // corresponds to size of dim * datasetSizeBuffers
    this->kernelDataset = std::vector<T>(this->preparedDataset.getNrows() *
                                         this->preparedDataset.getNcols());

    for (size_t i = 0; i < this->preparedDataset.getSize(); i++) {
      this->kernelDataset[i] = (T) this->preparedDataset[i];
    }

    for (size_t deviceIndex = 0; deviceIndex < devices.size(); deviceIndex++) {
      json::Node &platformConfiguration =
          (*parameters)["PLATFORMS"][devices[deviceIndex]->platformName];
      json::Node &deviceConfiguration =
          platformConfiguration["DEVICES"][devices[deviceIndex]->deviceName];
      json::Node &kernelConfiguration = deviceConfiguration
          ["KERNELS"][StreamingModOCLUnified::Configuration::getKernelName()];

      multKernels.emplace_back(devices[deviceIndex], dims, this->manager,
                               kernelConfiguration, queueLoadBalancerMult);

      multTransposeKernels.emplace_back(devices[deviceIndex], dims,
                                        this->manager, kernelConfiguration,
                                        queueLoadBalancerMultTrans);
    }

    // create the kernel specific data structures
    // also sets the correct padded grid size
    this->prepare();
  }

  ~OperationMultiEvalStreamingModOCLUnified() {}

  void mult(sgpp::base::DataVector &alpha,
            sgpp::base::DataVector &result) override {
    this->prepare();

    this->myTimer.start();

    size_t gridFrom = 0;
    size_t gridTo = this->gridSizePadded;
    size_t datasetFrom = 0;
    size_t datasetTo = this->datasetSizePadded;

    queueLoadBalancerMult->initialize(datasetFrom, datasetTo,
                                      overallDataBlockingSize);

    std::vector<T> alphaArray(this->gridSizePadded);

    for (size_t i = 0; i < alpha.getSize(); i++) {
      alphaArray[i] = (T)alpha[i] * scaling[i];
    }

    for (size_t i = alpha.getSize(); i < this->gridSizePadded; i++) {
      alphaArray[i] = static_cast<T>(0);
    }

    std::vector<T> resultArray(this->datasetSizePadded);
    std::fill(resultArray.begin(), resultArray.end(), 0.0);

    std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();

    int oldThreads = omp_get_max_threads();
    omp_set_num_threads(static_cast<int>(devices.size()));

#pragma omp parallel
    {
      size_t threadId = omp_get_thread_num();
      this->multKernels[threadId].mult(
          this->level, this->index, this->kernelDataset, alphaArray,
          resultArray, gridFrom, gridTo, datasetFrom, datasetTo);
    }
    end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;

    if (verbose) {
      std::cout << "duration mult ocl mod: " << elapsed_seconds.count()
                << std::endl;
    }

    result.resize(this->datasetSizePadded);
    for (size_t i = 0; i < result.getSize(); i++) {
      result[i] = resultArray[i];
    }

    // restore old value of OMP_NUM_THREADS
    omp_set_num_threads(oldThreads);

    this->duration = this->myTimer.stop();
  }

  void multTranspose(sgpp::base::DataVector &source,
                     sgpp::base::DataVector &result) override {
    this->prepare();

    this->myTimer.start();

    size_t gridFrom = 0;
    size_t gridTo = this->gridSizePadded;
    size_t datasetFrom = 0;
    size_t datasetTo = this->datasetSizePadded;

    queueLoadBalancerMultTrans->initialize(gridFrom, gridTo, overallGridBlockingSize);

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
          this->level, this->index, this->kernelDataset, sourceArray,
          resultArray, gridFrom, gridTo, datasetFrom, datasetTo);
    }
    end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;
    if (verbose) {
      std::cout << "duration multTranspose ocl mod: " << elapsed_seconds.count()
                << std::endl;
    }

    result.resize(this->gridSizePadded);
    for (size_t i = 0; i < result.getSize(); i++) {
      result[i] = resultArray[i] * scaling[i];
    }

    // restore old value of OMP_NUM_THREADS
    omp_set_num_threads(oldThreads);

    this->duration = this->myTimer.stop();
  }

  double getDuration() { return this->duration; }

  void prepare() override {
    this->recalculateLevelAndIndex();

    for (size_t deviceIndex = 0; deviceIndex < devices.size(); deviceIndex++) {
      this->multKernels[deviceIndex].resetKernel();
    }

    for (size_t deviceIndex = 0; deviceIndex < devices.size(); deviceIndex++) {
      this->multTransposeKernels[deviceIndex].resetKernel();
    }
  }

private:
  void padDataset(sgpp::base::DataMatrix &dataset) {
    datasetSizeUnpadded = dataset.getNrows();

    size_t commonDatasetPadding = calculateCommonDatasetPadding();

    // Assure that data has a even number of instances -> padding might be
    // needed
    size_t remainder = datasetSizeUnpadded % commonDatasetPadding;
    // round up to next number divisible by padding (distributable to threads)
    // then add another padding (for irregular schedules)
    size_t padding = commonDatasetPadding - remainder;

    // excluding the additional padding for irregular schedules
    datasetSizePadded = dataset.getNrows() + padding;

    sgpp::base::DataVector lastRow(dataset.getNcols());
    dataset.getRow(datasetSizeUnpadded - 1, lastRow);
    dataset.resize(datasetSizePadded);

    for (size_t i = datasetSizeUnpadded; i < datasetSizePadded; i++) {
      dataset.setRow(i, lastRow);
    }
  }

  void recalculateLevelAndIndex() {
    base::GridStorage &storage = grid.getStorage();

    size_t remainder = storage.getSize() % overallGridBlockingSize;
    size_t padding = 0;

    if (remainder != 0) {
      padding = overallGridBlockingSize - remainder;
    }

    gridSizeUnpadded = storage.getSize();
    gridSizePadded = storage.getSize() + padding;

    level = std::vector<T>(gridSizePadded * dims);
    index = std::vector<T>(gridSizePadded * dims);
    scaling = std::vector<T>(gridSizePadded);

    base::HashGridPoint::level_type curLevel;
    base::HashGridPoint::index_type curIndex;

    for (size_t i = 0; i < storage.getSize(); i++) {
      base::HashGridPoint &gridPoint = storage.getPoint(i);
      //      std::cout << "--------------" << std::endl;

      T scalingFactor = 1.0;
      for (size_t d = 0; d < dims; d++) {
        gridPoint.get(d, curLevel, curIndex);
        //        std::cout << "level: " << curLevel << " index: " << curIndex
        //        << std::endl;

        // handle the special cases for the extrapolating grid points first
        if (curLevel == 1 && curIndex == 1) {
          level[i * dims + d] = static_cast<T>(0); // special value to
                                                   // differentiate between
                                                   // border huts on level 1
          index[i * dims + d] = static_cast<T>(0);
          //          std::cout << "mod level: " << level[i * dims + d] <<
          //          std::endl;
          //          std::cout << "mod index: " << index[i * dims + d] <<
          //          std::endl;
          scalingFactor *= static_cast<T>(1.0);
        } else if (curLevel > 1 && curIndex == 1) {
          level[i * dims + d] = static_cast<T>(1 << (curLevel - 1));
          index[i * dims + d] = static_cast<T>(0);
          //          std::cout << "mod level: " << level[i * dims + d] <<
          //          std::endl;
          //          std::cout << "mod index: " << index[i * dims + d] <<
          //          std::endl;
          scalingFactor *= static_cast<T>(2.0);
        } else if (curLevel > 1 &&
                   curIndex == static_cast<uint32_t>(1 << curLevel) - 1) {
          level[i * dims + d] = static_cast<T>(1 << (curLevel - 1));
          index[i * dims + d] = static_cast<T>(1 << (curLevel - 1));
          //          std::cout << "mod level: " << level[i * dims + d] <<
          //          std::endl;
          //          std::cout << "mod index: " << index[i * dims + d] <<
          //          std::endl;
          scalingFactor *= static_cast<T>(2.0);
        } else {
          level[i * dims + d] = static_cast<T>(1 << curLevel);
          index[i * dims + d] = static_cast<T>(curIndex);
          //          std::cout << "mod level: " << level[i * dims + d] <<
          //          std::endl;
          //          std::cout << "mod index: " << index[i * dims + d] <<
          //          std::endl;
          scalingFactor *= static_cast<T>(1.0);
        }
      }
      scaling[i] = scalingFactor;
    }

    // grid points are disabled via surplus set to zero
    for (size_t i = storage.getSize(); i < gridSizePadded; i++) {
      for (size_t dim = 0; dim < storage.getDimension(); dim++) {
        level[i * dims + dim] =
            1.0; // same treatment as constant basis function
        index[i * dims + dim] = 0.0;
      }
      scaling[i] = 1.0; // don't care
    }
  }

  size_t calculateCommonDatasetPadding() {
    size_t commonPaddingRequiredment = 1;
    for (size_t deviceIndex = 0; deviceIndex < devices.size(); deviceIndex++) {
      json::Node &platformConfiguration =
          (*parameters)["PLATFORMS"][devices[deviceIndex]->platformName];
      json::Node &deviceConfiguration =
          platformConfiguration["DEVICES"][devices[deviceIndex]->deviceName];
      json::Node &kernelConfiguration = deviceConfiguration
          ["KERNELS"][StreamingModOCLUnified::Configuration::getKernelName()];

      commonPaddingRequiredment =
          std::max(commonPaddingRequiredment,
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
      json::Node &kernelConfiguration = deviceConfiguration
          ["KERNELS"][StreamingModOCLUnified::Configuration::getKernelName()];

      commonPaddingRequiredment = std::max(
          commonPaddingRequiredment,
          kernelConfiguration["KERNEL_TRANS_GRID_BLOCK_SIZE"].getUInt() *
              kernelConfiguration["LOCAL_SIZE"].getUInt());
    }
    return commonPaddingRequiredment;
  }
};

} // namespace datadriven
} // namespace sgpp
