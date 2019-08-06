// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <omp.h>

#include <chrono>

#include <sgpp/base/operation/hash/OperationMultipleEval.hpp>
#include <sgpp/base/tools/SGppStopwatch.hpp>
#include <sgpp/base/exception/operation_exception.hpp>
#include <sgpp/base/grid/type/BsplineGrid.hpp>
#include <sgpp/globaldef.hpp>
#include <sgpp/base/opencl/OCLOperationConfiguration.hpp>
#include <sgpp/base/opencl/OCLManager.hpp>
#include <sgpp/datadriven/operation/hash/OperationMultipleEvalStreamingBSplineOCL/StreamingBSplineOCLKernelImpl.hpp>

namespace sgpp {
namespace datadriven {

template <typename T>
class OperationMultiEvalStreamingBSplineOCL : public base::OperationMultipleEval {
 protected:
  size_t dims;
  sgpp::base::DataMatrix preparedDataset;
  std::shared_ptr<base::OCLOperationConfiguration> parameters;
  T* kernelDataset = nullptr;
  size_t datasetSize = 0;
  /// Member to store the sparse grid's levels for better vectorization
  T* level = nullptr;
  /// Member to store the sparse grid's indices for better vectorization
  T* index = nullptr;
  size_t gridSize = 0;
  /// Timer object to handle time measurements
  sgpp::base::SGppStopwatch myTimer;

  base::GridStorage* storage;

  double duration;

  std::shared_ptr<base::OCLManager> manager;
  std::unique_ptr<StreamingBSplineOCLKernelImpl<T>> kernel;

 public:
  OperationMultiEvalStreamingBSplineOCL(base::Grid& grid, base::DataMatrix& dataset,
                                        std::shared_ptr<base::OCLOperationConfiguration> parameters)
      : OperationMultipleEval(grid, dataset),
        preparedDataset(dataset),
        parameters(parameters),
        myTimer(sgpp::base::SGppStopwatch()),
        duration(-1.0) {
    this->manager = std::make_shared<base::OCLManager>(parameters);

    this->dims = dataset.getNcols();  // be aware of transpose!
    this->kernel =
        std::unique_ptr<StreamingBSplineOCLKernelImpl<T>>(new StreamingBSplineOCLKernelImpl<T>(
            dims, dynamic_cast<base::BsplineGrid&>(grid).getDegree(), this->manager, parameters));

    this->storage = &grid.getStorage();
    this->padDataset(this->preparedDataset);
    this->preparedDataset.transpose();
    this->datasetSize = this->preparedDataset.getNcols();

    //    std::cout << "dims: " << this->dims << std::endl;
    //    std::cout << "padded instances: " << this->datasetSize << std::endl;

    this->kernelDataset =
        new T[this->preparedDataset.getNrows() * this->preparedDataset.getNcols()];

    for (size_t i = 0; i < this->preparedDataset.getSize(); i++) {
      this->kernelDataset[i] = (T) this->preparedDataset[i];
    }

    // create the kernel specific data structures
    this->prepare();
  }

  ~OperationMultiEvalStreamingBSplineOCL() {
    if (this->level != nullptr) {
      delete this->level;
    }

    if (this->index != nullptr) {
      delete this->index;
    }

    if (this->kernelDataset != nullptr) {
      delete this->kernelDataset;
    }
  }

  void mult(sgpp::base::DataVector& alpha, sgpp::base::DataVector& result) override {
    this->myTimer.start();

    size_t gridFrom = 0;
    size_t gridTo = this->gridSize;
    size_t datasetFrom = 0;
    size_t datasetTo = this->datasetSize;

    T* alphaArray = new T[this->gridSize];

    for (size_t i = 0; i < alpha.getSize(); i++) {
      alphaArray[i] = (T)alpha[i];
    }

    for (size_t i = alpha.getSize(); i < this->gridSize; i++) {
      alphaArray[i] = 0.0;
    }

    T* resultArray = new T[this->datasetSize];

    for (size_t i = 0; i < this->datasetSize; i++) {
      resultArray[i] = 0.0;
    }

    std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();
    this->kernel->mult(this->level, this->index, this->gridSize, this->kernelDataset,
                       this->datasetSize, alphaArray, resultArray, gridFrom, gridTo, datasetFrom,
                       datasetTo);
    end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;
    std::cout << "duration mult ocl B-spline: " << elapsed_seconds.count() << std::endl;

    result.resize(this->datasetSize);
    for (size_t i = 0; i < result.getSize(); i++) {
      result[i] = resultArray[i];
    }

    delete alphaArray;
    delete resultArray;
    this->duration = this->myTimer.stop();
  }

  void multTranspose(sgpp::base::DataVector& source, sgpp::base::DataVector& result) override {
    this->myTimer.start();

    size_t gridFrom = 0;
    size_t gridTo = this->gridSize;
    size_t datasetFrom = 0;
    size_t datasetTo = this->datasetSize;

    T* sourceArray = new T[this->datasetSize];

    for (size_t i = 0; i < source.getSize(); i++) {
      sourceArray[i] = (T)source[i];
    }

    for (size_t i = source.getSize(); i < this->datasetSize; i++) {
      sourceArray[i] = 0.0;
    }

    T* resultArray = new T[this->gridSize];

    for (size_t i = 0; i < this->gridSize; i++) {
      resultArray[i] = 0.0;
    }

    std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();
    this->kernel->multTranspose(this->level, this->index, this->gridSize, this->kernelDataset,
                                this->preparedDataset.getNcols(), sourceArray, resultArray,
                                gridFrom, gridTo, datasetFrom, datasetTo);
    end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;
    std::cout << "duration multTranspose ocl B-spline: " << elapsed_seconds.count() << std::endl;

    result.resize(this->gridSize);
    for (size_t i = 0; i < result.getSize(); i++) {
      result[i] = resultArray[i];
    }

    delete sourceArray;
    delete resultArray;
    this->duration = this->myTimer.stop();
  }

  double getDuration() { return this->duration; }

  void prepare() override {
    this->recalculateLevelAndIndex();

    this->kernel->resetKernel();

    //    std::cout << "gridSize: " << this->gridSize << std::endl;
  }

 private:
  size_t padDataset(sgpp::base::DataMatrix& dataset) {
    size_t vecWidth =
        (*parameters)["LOCAL_SIZE"].getUInt() * (*parameters)["KERNEL_DATA_BLOCK_SIZE"].getUInt();

    // Assure that data has a even number of instances -> padding might be needed
    size_t remainder = dataset.getNrows() % vecWidth;
    size_t loopCount = vecWidth - remainder;

    if (loopCount != vecWidth) {
      sgpp::base::DataVector lastRow(dataset.getNcols());
      size_t oldSize = dataset.getNrows();
      dataset.getRow(dataset.getNrows() - 1, lastRow);
      dataset.resize(dataset.getNrows() + loopCount);

      for (size_t i = 0; i < loopCount; i++) {
        dataset.setRow(oldSize + i, lastRow);
      }
    }

    return dataset.getNrows();
  }

  void recalculateLevelAndIndex() {
    if (this->level != nullptr) delete this->level;

    if (this->index != nullptr) delete this->index;

    uint32_t localWorkSize = (uint32_t)(*parameters)["LOCAL_SIZE"].getUInt();

    size_t remainder = this->storage->getSize() % localWorkSize;
    size_t padding = 0;

    if (remainder != 0) {
      padding = localWorkSize - remainder;
    }

    this->gridSize = this->storage->getSize() + padding;

    sgpp::base::DataMatrix* levelMatrix = new sgpp::base::DataMatrix(this->gridSize, this->dims);
    sgpp::base::DataMatrix* indexMatrix = new sgpp::base::DataMatrix(this->gridSize, this->dims);

    this->storage->getLevelIndexArraysForEval(*levelMatrix, *indexMatrix);

    for (size_t i = this->storage->getSize(); i < this->gridSize; i++) {
      for (size_t j = 0; j < this->storage->getDimension(); j++) {
        levelMatrix->set(i, j, 1.0);
        indexMatrix->set(i, j, 1.0);
      }
    }

    this->level = new T[this->gridSize * this->dims];
    this->index = new T[this->gridSize * this->dims];

    for (size_t i = 0; i < this->gridSize * this->dims; i++) {
      this->level[i] = (T)(*levelMatrix)[i];
      this->index[i] = (T)(*indexMatrix)[i];
    }

    delete levelMatrix;
    delete indexMatrix;
  }
};
}  // namespace datadriven
}  // namespace sgpp
