// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <chrono>
#include <omp.h>

#include <sgpp/base/operation/hash/OperationMultipleEval.hpp>
#include <sgpp/base/tools/SGppStopwatch.hpp>
#include <sgpp/base/exception/operation_exception.hpp>
#include <sgpp/globaldef.hpp>
#include <sgpp/base/opencl/OCLOperationConfiguration.hpp>
#include <sgpp/base/opencl/OCLManager.hpp>
#include <sgpp/datadriven/operation/hash/OperationMultipleEvalStreamingModOCLMask/StreamingModOCLMaskKernelImpl.hpp>

namespace SGPP {
  namespace datadriven {

    template<typename T>
    class OperationMultiEvalStreamingModOCLMask: public base::OperationMultipleEval {
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
        /// Member to store the sparse grid's mask for better vectorization
        std::vector<T> mask;
        /// Member to store the sparse grid's offset for better vectorization
        std::vector<T> offset;
        size_t gridSize;
        /// Timer object to handle time measurements
        SGPP::base::SGppStopwatch myTimer;

        base::GridStorage* storage;

        float_t duration;

        std::shared_ptr<base::OCLManager> manager;
        std::unique_ptr<StreamingModOCLMaskKernelImpl<T>> kernel;
      public:

        OperationMultiEvalStreamingModOCLMask(base::Grid& grid, base::DataMatrix& dataset,
                                              std::shared_ptr<base::OCLOperationConfiguration> parameters) :
          OperationMultipleEval(grid, dataset), preparedDataset(dataset), parameters(parameters), myTimer(
            SGPP::base::SGppStopwatch()), duration(-1.0) {

          this->manager = std::make_shared<base::OCLManager>(parameters);

          this->dims = dataset.getNcols(); //be aware of transpose!
          this->kernel = std::unique_ptr<StreamingModOCLMaskKernelImpl<T>>(
                           new StreamingModOCLMaskKernelImpl<T>(dims, this->manager, parameters));

          this->storage = grid.getStorage();
          this->gridSize = this->storage->size();
          this->padDataset(this->preparedDataset);
          this->preparedDataset.transpose();
          this->datasetSize = this->preparedDataset.getNcols();

          //    std::cout << "dims: " << this->dims << std::endl;
          //    std::cout << "padded instances: " << this->datasetSize << std::endl;

          this->kernelDataset = std::vector<T>(this->preparedDataset.getNrows() * this->preparedDataset.getNcols());
          //        this->kernelDataset = new T[this->preparedDataset.getNrows() * this->preparedDataset.getNcols()];

          for (size_t i = 0; i < this->preparedDataset.getSize(); i++) {
            this->kernelDataset[i] = (T) this->preparedDataset[i];
          }

          //create the kernel specific data structures
          this->prepare();
        }

        ~OperationMultiEvalStreamingModOCLMask() {
          //        if (this->level != nullptr) {
          //            delete this->level;
          //        }
          //
          //        if (this->index != nullptr) {
          //            delete this->index;
          //        }
          //
          //        if (this->mask != nullptr) {
          //            delete this->mask;
          //        }
          //
          //        if (this->offset != nullptr) {
          //            delete this->offset;
          //        }
          //
          //        if (this->kernelDataset != nullptr) {
          //            delete this->kernelDataset;
          //        }
        }

        void mult(SGPP::base::DataVector& alpha,
                  SGPP::base::DataVector& result) override {
          this->myTimer.start();

          size_t gridFrom = 0;
          size_t gridTo = this->gridSize;
          size_t datasetFrom = 0;
          size_t datasetTo = this->datasetSize;

          std::vector<T> alphaArray(this->gridSize);

          for (size_t i = 0; i < alpha.getSize(); i++) {
            alphaArray[i] = (T) alpha[i];
          }

          for (size_t i = alpha.getSize(); i < this->gridSize; i++) {
            alphaArray[i] = 0.0;
          }

          std::vector<T> resultArray(this->datasetSize);

          for (size_t i = 0; i < this->datasetSize; i++) {
            resultArray[i] = 0.0;
          }

          std::chrono::time_point<std::chrono::system_clock> start, end;
          start = std::chrono::system_clock::now();
          this->kernel->mult(this->level, this->index, this->mask, this->offset, this->gridSize, this->kernelDataset, this->datasetSize, alphaArray,
                             resultArray, gridFrom, gridTo, datasetFrom, datasetTo);
          end = std::chrono::system_clock::now();
          std::chrono::duration<double> elapsed_seconds = end - start;
          std::cout << "duration mult ocl mod: " << elapsed_seconds.count() << std::endl;

          for (size_t i = 0; i < result.getSize(); i++) {
            result[i] = resultArray[i];
          }

          this->duration = this->myTimer.stop();
        }

        void multTranspose(
          SGPP::base::DataVector& source,
          SGPP::base::DataVector& result) override {
          this->myTimer.start();

          size_t gridFrom = 0;
          size_t gridTo = this->gridSize;
          size_t datasetFrom = 0;
          size_t datasetTo = this->datasetSize;

          std::vector<T> sourceArray(this->datasetSize);

          for (size_t i = 0; i < source.getSize(); i++) {
            sourceArray[i] = (T) source[i];
          }

          for (size_t i = source.getSize(); i < this->datasetSize; i++) {
            sourceArray[i] = 0.0;
          }

          std::vector<T> resultArray(this->gridSize);

          for (size_t i = 0; i < this->gridSize; i++) {
            resultArray[i] = 0.0;
          }

          std::chrono::time_point<std::chrono::system_clock> start, end;
          start = std::chrono::system_clock::now();
          this->kernel->multTranspose(this->level, this->index, this->mask, this->offset, this->gridSize, this->kernelDataset,
                                      this->preparedDataset.getNcols(), sourceArray, resultArray, gridFrom, gridTo, datasetFrom, datasetTo);
          end = std::chrono::system_clock::now();
          std::chrono::duration<double> elapsed_seconds = end - start;
          std::cout << "duration multTranspose ocl mod: " << elapsed_seconds.count() << std::endl;

          for (size_t i = 0; i < result.getSize(); i++) {
            result[i] = resultArray[i];
          }

          this->duration = this->myTimer.stop();
        }

        float_t getDuration() {
          return this->duration;
        }

        void prepare() override {
          this->recalculateLevelIndexMask();

          this->kernel->resetKernel();

          //    std::cout << "gridSize: " << this->gridSize << std::endl;
        }

      private:

        size_t padDataset(
          SGPP::base::DataMatrix& dataset) {

          size_t vecWidth = (*parameters)["LOCAL_SIZE"].getUInt();
          //                * (*parameters)["KERNEL_DATA_BLOCKING_SIZE"].getUInt();

          // Assure that data has a even number of instances -> padding might be needed
          size_t remainder = dataset.getNrows() % vecWidth;
          size_t loopCount = vecWidth - remainder;

          if (loopCount != vecWidth) {
            SGPP::base::DataVector lastRow(dataset.getNcols());
            size_t oldSize = dataset.getNrows();
            dataset.getRow(dataset.getNrows() - 1, lastRow);
            dataset.resize(dataset.getNrows() + loopCount);

            for (size_t i = 0; i < loopCount; i++) {
              dataset.setRow(oldSize + i, lastRow);
            }
          }

          return dataset.getNrows();
        }

        /**
         * Converts this storage from AOS (array of structures) to SOA (structure of array)
         * with modification to speed up iterative function evaluation. The Level
         * array won't contain the levels, it contains the level to the power of two.
         *
         * The returned format is only useful for a multi-evaluation of modlinear grids
         */
        void recalculateLevelIndexMask() {

          uint32_t localWorkSize = (uint32_t) (*parameters)["LOCAL_SIZE"].getUInt();

          size_t remainder = this->storage->size() % localWorkSize;
          size_t padding = 0;

          if (remainder != 0) {
            padding = localWorkSize - remainder;
          }

          this->gridSize = this->storage->size() + padding;

          SGPP::base::HashGridIndex::level_type curLevel;
          SGPP::base::HashGridIndex::index_type curIndex;

          //TODO: update the other kernels with this style

          this->level = std::vector<T>(this->gridSize * this->dims);
          this->index = std::vector<T>(this->gridSize * this->dims);
          this->mask = std::vector<T>(this->gridSize * this->dims);
          this->offset = std::vector<T>(this->gridSize * this->dims);

          for (size_t i = 0; i < this->storage->size(); i++) {
            for (size_t dim = 0; dim < this->dims; dim++) {
              storage->get(i)->get(dim, curLevel, curIndex);

              if (curLevel == 1) {
                this->level[i * this->dims + dim] = 0.0;
                this->index[i * this->dims + dim] = 0.0;

                if (std::is_same<T, double>::value) {
                  uint64_t intmask = 0x0000000000000000;
                  this->mask[i * this->dims + dim] = *reinterpret_cast<T*>(&intmask);
                } else {
                  uint32_t intmask = 0x00000000;
                  this->mask[i * this->dims + dim] = *reinterpret_cast<T*>(&intmask);
                }

                this->offset[i * this->dims + dim] = 1.0;
              } else if (curIndex == 1) {
                this->level[i * this->dims + dim] = static_cast<T>(-1.0) * static_cast<T>(1 << curLevel);
                this->index[i * this->dims + dim] = 0.0;

                if (std::is_same<T, double>::value) {
                  uint64_t intmask = 0x0000000000000000;
                  this->mask[i * this->dims + dim] = *reinterpret_cast<T*>(&intmask);
                } else {
                  uint32_t intmask = 0x00000000;
                  this->mask[i * this->dims + dim] = *reinterpret_cast<T*>(&intmask);
                }

                this->offset[i * this->dims + dim] = 2.0;
              } else if (curIndex == static_cast<SGPP::base::HashGridIndex::level_type>(((1 << curLevel) - 1))) {
                this->level[i * this->dims + dim] = static_cast<T>(1 << curLevel);
                this->index[i * this->dims + dim] = static_cast<T>(curIndex);

                if (std::is_same<T, double>::value) {
                  uint64_t intmask = 0x0000000000000000;
                  this->mask[i * this->dims + dim] = *reinterpret_cast<T*>(&intmask);
                } else {
                  uint32_t intmask = 0x00000000;
                  this->mask[i * this->dims + dim] = *reinterpret_cast<T*>(&intmask);
                }

                this->offset[i * this->dims + dim] = 1.0;
              } else {
                this->level[i * this->dims + dim] = static_cast<T>(1 << curLevel);
                this->index[i * this->dims + dim] = static_cast<T>(curIndex);

                if (std::is_same<T, double>::value) {
                  uint64_t intmask = 0x8000000000000000;
                  this->mask[i * this->dims + dim] = *reinterpret_cast<T*>(&intmask);
                } else {
                  uint32_t intmask = 0x80000000;
                  this->mask[i * this->dims + dim] = *reinterpret_cast<T*>(&intmask);
                }

                this->offset[i * this->dims + dim] = 1.0;
              }
            }
          }

          for (size_t i = this->storage->size(); i < this->gridSize; i++) {
            for (size_t dim = 0; dim < this->dims; dim++) {
              this->level[i * this->dims + dim] = 0;
              this->index[i * this->dims + dim] = 0;

              if (std::is_same<T, double>::value) {
                uint64_t intmask = 0x0000000000000000;
                this->mask[i * this->dims + dim] = *reinterpret_cast<T*>(&intmask);
              } else {
                uint32_t intmask = 0x00000000;
                this->mask[i * this->dims + dim] = *reinterpret_cast<T*>(&intmask);
              }

              this->offset[i * this->dims + dim] = 1.0;
            }
          }
        }
    };

  }
}
