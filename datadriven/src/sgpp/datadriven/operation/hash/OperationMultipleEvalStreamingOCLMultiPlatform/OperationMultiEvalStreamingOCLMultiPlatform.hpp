// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <omp.h>

#include <algorithm>
#include <chrono>
#include <mutex>  // NOLINT(build/c++11)
#include <vector>

#include <sgpp/datadriven/operation/hash/OperationMultipleEvalStreamingOCLMultiPlatform/Configuration.hpp>
#include <sgpp/datadriven/operation/hash/OperationMultipleEvalStreamingOCLMultiPlatform/KernelMult.hpp>
#include <sgpp/datadriven/operation/hash/OperationMultipleEvalStreamingOCLMultiPlatform/KernelMultTranspose.hpp>
#include <sgpp/base/exception/operation_exception.hpp>
#include <sgpp/base/opencl/OCLManager.hpp>
#include <sgpp/base/opencl/OCLOperationConfiguration.hpp>
#include <sgpp/base/operation/hash/OperationMultipleEval.hpp>
#include <sgpp/base/tools/QueueLoadBalancerOpenMP.hpp>
#include <sgpp/base/tools/SGppStopwatch.hpp>
#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace datadriven {
namespace StreamingOCLMultiPlatform {

/**
 * This class provides an operation for evaluating multiple grid points in the domain and doing
 * least squares data mining.
 * This algorithmic variant uses the streaming algorithm for evaluation.
 * It uses high performance OpenCL kernels and is well-suited for large irregular datasets and
 * grids.
 * This class manages one OpenCL kernel for each devices configured using the
 * OCLOperationConfiguration.
 * When a operation is called it triggers the device work by using OpenMP and delegating the work to
 * instances of the kernels.
 * Furthermore, this class converts the received grid and dataset into a representation that is
 * suited for the streaming algorithm.
 *
 * @see base::OperationMultipleEval
 * @see StreamingOCLMultiPlatform::KernelMult
 * @see StreamingOCLMultiPlatform::KernelMultTranspose
 */
template <typename T>
class OperationMultiEvalStreamingOCLMultiPlatform : public base::OperationMultipleEval {
 protected:
  size_t dims;

  base::DataMatrix preparedDataset;

  std::shared_ptr<base::OCLOperationConfiguration> parameters;
  std::vector<T> kernelDataset;

  size_t datasetSizeUnpadded;
  size_t datasetSizePadded;

  //  // includes padding
  //  size_t datasetSize;

  // Member to store the sparse grid's levels for better vectorization
  std::vector<T> level;
  // Member to store the sparse grid's indices for better vectorization
  std::vector<T> index;

  size_t gridSizeUnpadded;
  size_t gridSizePadded;

  //  // includes padding
  //  size_t gridSize;

  // Timer object to handle time measurements
  base::SGppStopwatch myTimer;

  double duration;

  std::shared_ptr<base::QueueLoadBalancerOpenMP> queueLoadBalancerMult;
  std::shared_ptr<base::QueueLoadBalancerOpenMP> queueLoadBalancerMultTranspose;

  std::shared_ptr<base::OCLManagerMultiPlatform> manager;
  std::vector<std::shared_ptr<base::OCLDevice>> devices;

  std::vector<StreamingOCLMultiPlatform::KernelMult<T>> multKernels;
  std::vector<StreamingOCLMultiPlatform::KernelMultTranspose<T>> multTransposeKernels;

  //  json::Node &configuration;

  bool verbose;

  size_t commonDatasetPadding;
  size_t commonGridPadding;

 public:
  /**
   * Creates a new instance of the OperationMultiEvalStreamingOCLMultiPlatform class.
   * This class should not be created directly, instead the datadriven operator factory should be
   * used or at least the factory method.
   *
   * @see createStreamingOCLMultiPlatformConfigured
   *
   * @param grid The grid to evaluate
   * @param dataset The datapoints to evaluate
   * @param manager The OpenCL manager that manages OpenCL internels for this kernel
   * @param parameters The configuration of the kernel leading to different compute kernels
   */
  OperationMultiEvalStreamingOCLMultiPlatform(
      base::Grid &grid, base::DataMatrix &dataset,
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

    queueLoadBalancerMult = std::make_shared<base::QueueLoadBalancerOpenMP>();
    queueLoadBalancerMultTranspose = std::make_shared<base::QueueLoadBalancerOpenMP>();

    // initialized in padDataset
    datasetSizeUnpadded = 0;
    datasetSizePadded = 0;

    // initialized in prepare
    gridSizeUnpadded = 0;
    gridSizePadded = 0;

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

  /**
   * Destructor
   */
  ~OperationMultiEvalStreamingOCLMultiPlatform() {}

  /**
   * Performs the MultiEval operation \f$v:= B^T \alpha\f$.
   *
   * @param alpha The surpluses of the grid
   * @param result A vector that contains the result in the order of the dataset
   */
  void mult(base::DataVector &alpha, base::DataVector &result) override {
    this->mult(alpha, result, 0, this->datasetSizePadded);
  }

  void mult(base::DataVector &alpha, base::DataVector &result, size_t startIndexData,
            size_t endIndexData) override {
    // ensure padding requirements are fulfilled by start and end index
    size_t startIndexDataPadded = this->padIndexDown(startIndexData, commonDatasetPadding);
    size_t endIndexDataPadded = this->padIndexUp(endIndexData, commonDatasetPadding);

    this->myTimer.start();

    if (startIndexData > endIndexData) {
      std::stringstream errorString;
      errorString << "Error: cannot process negative data set input range" << std::endl;
      throw base::operation_exception(errorString.str());
    }

    this->prepare();

    size_t gridFrom = 0;
    size_t gridTo = this->gridSizePadded;
    size_t datasetFrom = startIndexDataPadded;
    size_t datasetTo = endIndexDataPadded;

    queueLoadBalancerMult->initialize(datasetFrom, datasetTo, commonDatasetPadding);

    std::vector<T> alphaArray(this->gridSizePadded);

    for (size_t i = 0; i < alpha.getSize(); i++) {
      alphaArray[i] = (T)alpha[i];
    }

    for (size_t i = alpha.getSize(); i < this->gridSizePadded; i++) {
      alphaArray[i] = 0.0;
    }

    // additional padding to allow for devices with different block sizes
    std::vector<T> resultArray(datasetTo - datasetFrom);

    std::fill(resultArray.begin(), resultArray.end(), 0.0);

    int oldThreads = omp_get_max_threads();
    omp_set_num_threads(static_cast<int>(devices.size()));

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
        std::call_once(onceFlag,
                       [&]() { exceptionPtr = std::current_exception(); });  // NOLINT(build/c++11)
      }
    }

    if (exceptionPtr) {
      std::rethrow_exception(exceptionPtr);
    }

    result.resize(endIndexData - startIndexData);
    for (size_t i = 0; i < result.getSize(); i++) {
      result[i] = resultArray[(startIndexData - startIndexDataPadded) + i];
    }

    // restore old value of OMP_NUM_THREADS
    omp_set_num_threads(oldThreads);

    this->duration = this->myTimer.stop();

    for (StreamingOCLMultiPlatform::KernelMult<T> &kernel : multKernels) {
      this->duration -= kernel.getBuildDuration();
    }

    if (verbose) {
      std::cout << "duration mult ocl: " << this->duration << std::endl;
    }
  }

  /**
   * Performs the transposed MultiEval operation  \f$v':= B v\f$.
   *
   * @param source The vector \f$v\f$
   * @param result The result of the matrix vector multiplication in the order of grid (of the alpha
   * vector)
   */
  void multTranspose(base::DataVector &source, base::DataVector &result) override {
    this->multTranspose(source, result, 0, this->gridSizePadded);
  }

  void multTranspose(base::DataVector &source, base::DataVector &result, size_t startIndexGrid,
                     size_t endIndexGrid) override {
    this->myTimer.start();

    // ensure padding requirements are fulfilled by start and end index
    size_t startIndexGridPadded = this->padIndexDown(startIndexGrid, commonGridPadding);
    size_t endIndexGridPadded = this->padIndexUp(endIndexGrid, commonGridPadding);
    std::cout << "endIndexGrid: " << endIndexGrid << std::endl;
    std::cout << "commonGridPadding: " << commonGridPadding << std::endl;
    std::cout << "endIndexGridPadded: " << endIndexGridPadded << std::endl;

    if (startIndexGrid > endIndexGrid) {
      std::stringstream errorString;
      errorString << "Error: cannot process negative grid input range" << std::endl;
      throw base::operation_exception(errorString.str());
    }

    this->prepare();

    size_t gridFrom = startIndexGridPadded;
    size_t gridTo = endIndexGridPadded;
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

    std::vector<T> resultArray(gridTo - gridFrom);

    std::fill(resultArray.begin(), resultArray.end(), 0.0);

    int oldThreads = omp_get_max_threads();
    omp_set_num_threads(static_cast<int>(devices.size()));

    std::once_flag onceFlag;
    std::exception_ptr exceptionPtr;

#pragma omp parallel
    {
      size_t threadId = omp_get_thread_num();

      try {
        this->multTransposeKernels[threadId].multTranspose(this->level, this->index,
                                                           this->kernelDataset, sourceArray,
                                                           resultArray, datasetFrom, datasetTo);
      } catch (...) {
        // store the first exception thrown for rethrow
        std::call_once(onceFlag,
                       [&]() { exceptionPtr = std::current_exception(); });  // NOLINT(build/c++11)
      }
    }

    if (exceptionPtr) {
      std::rethrow_exception(exceptionPtr);
    }

    result.resize(endIndexGrid - startIndexGrid);
    for (size_t i = 0; i < result.getSize(); i++) {
      result[i] = resultArray[(startIndexGrid - startIndexGridPadded) + i];
    }

    // restore old value of OMP_NUM_THREADS
    omp_set_num_threads(oldThreads);

    this->duration = this->myTimer.stop();

    for (StreamingOCLMultiPlatform::KernelMultTranspose<T> &kernelTranspose :
         multTransposeKernels) {
      this->duration -= kernelTranspose.getBuildDuration();
    }

    if (verbose) {
      std::cout << "duration multTranspose ocl: " << duration << std::endl;
    }
  }

  /**
   * @return The duration of the last call to mult or multTranspose
   */
  double getDuration() override { return this->duration; }

  /**
   * Creates the internal data structures used by the algorithm. Needs to be called every time the
   * grid changes e.g., due to refinement.
   */
  void prepare() override { this->recalculateLevelAndIndex(); }

 private:
  /**
   * Pads the dataset according to the requirement of the parameter configuration
   * @param dataset The dataset to be padded
   */
  void padDataset(base::DataMatrix &dataset) {
    // Assure that data has a even number of instances -> padding might be
    // needed
    size_t remainder = dataset.getNrows() % commonDatasetPadding;
    size_t padding = commonDatasetPadding - remainder;
    datasetSizeUnpadded = dataset.getNrows();
    datasetSizePadded = dataset.getNrows() + padding;

    // replicate last row for padding
    base::DataVector lastRow(dims);
    dataset.getRow(datasetSizeUnpadded - 1, lastRow);
    dataset.resize(datasetSizePadded);

    for (size_t i = datasetSizeUnpadded; i < datasetSizePadded; i++) {
      dataset.setRow(i, lastRow);
    }
  }

  /**
   * Creates the internal grid data structure which consists of two lists, one for the level and one
   * for the index values.
   * Has to be recalculated via prepare whenever the grid is changed.
   */
  void recalculateLevelAndIndex() {
    base::GridStorage &storage = grid.getStorage();

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

    for (size_t i = 0; i < storage.getSize(); i++) {
      base::HashGridPoint &gridPoint = storage.getPoint(i);
      for (size_t dim = 0; dim < dims; dim++) {
        gridPoint.get(dim, curLevel, curIndex);
        level[i * dims + dim] = static_cast<T>(1 << curLevel);
        index[i * dims + dim] = static_cast<T>(curIndex);
      }
    }

    for (size_t i = storage.getSize(); i < gridSizePadded; i++) {
      for (size_t dim = 0; dim < storage.getDimension(); dim++) {
        level[i * dims + dim] = 1.0;
        index[i * dims + dim] = 1.0;
      }
    }
  }

  /**
   * Calculates the padding requirements of the dataset according to the configuration.
   * Takes into account the configuration of individual devices.
   */
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

  /**
   * Calculates the padding requirements of the internal grid data structure according to the
   * configuration.
   * Takes into account the configuration of individual devices.
   */
  size_t calculateCommonGridPadding() {
    size_t commonPaddingRequirement = 1;
    for (size_t deviceIndex = 0; deviceIndex < devices.size(); deviceIndex++) {
      json::Node &platformConfiguration =
          (*parameters)["PLATFORMS"][devices[deviceIndex]->platformName];
      json::Node &deviceConfiguration =
          platformConfiguration["DEVICES"][devices[deviceIndex]->deviceName];
      json::Node &kernelConfiguration =
          deviceConfiguration["KERNELS"][StreamingOCLMultiPlatform::Configuration::getKernelName()];

      commonPaddingRequirement = std::max(
          commonPaddingRequirement, kernelConfiguration["KERNEL_TRANS_GRID_BLOCK_SIZE"].getUInt() *
                                        kernelConfiguration["LOCAL_SIZE"].getUInt());
    }
    return commonPaddingRequirement;
  }

  // pads an index to make it divisible by the provided blocksize ()
  size_t padIndexUp(size_t index, size_t blockSize) {
    size_t remainder = index % blockSize;
    if (remainder == 0) {
      return index;
    }
    return index + (blockSize - remainder);
  }

  // pads an index to make it divisible by the provided blocksize
  size_t padIndexDown(size_t index, size_t blockSize) {
    size_t remainder = index % blockSize;
    if (remainder == 0) {
      return index;
    }    
    return index - remainder;
  }
};

}  // namespace StreamingOCLMultiPlatform
}  // namespace datadriven
}  // namespace sgpp
