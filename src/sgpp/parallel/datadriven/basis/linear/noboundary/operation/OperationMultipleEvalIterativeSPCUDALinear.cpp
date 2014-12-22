/******************************************************************************
* Copyright (C) 2011 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "base/exception/operation_exception.hpp"

#include "parallel/datadriven/basis/linear/noboundary/operation/OperationMultipleEvalIterativeSPCUDALinear.hpp"
#include "parallel/datadriven/basis/common/CUDAKernels.hpp"

namespace sg {

  namespace parallel {

    OperationMultipleEvalIterativeSPCUDALinear::OperationMultipleEvalIterativeSPCUDALinear(sg::base::GridStorage* storage, sg::base::DataMatrixSP* dataset) : sg::parallel::OperationMultipleEvalVectorizedSP(storage, dataset), old_storage_size(0) {
      this->storage_ = storage;

      this->level_ = new sg::base::DataMatrixSP(storage->size(), storage->dim());
      this->index_ = new sg::base::DataMatrixSP(storage->size(), storage->dim());

      storage->getLevelIndexArraysForEval(*(this->level_), *(this->index_));

      myTimer = new sg::base::SGppStopwatch();

      this->old_storage_size = this->storage_->size();

      uploadGridSPCUDA(this->level_->getPointer(), this->index_->getPointer(), this->storage_->size(), this->storage_->dim());
      uploadDataSPCUDA(this->dataset_->getPointer(), this->dataset_->getNcols(), this->dataset_->getNrows());
    }

    OperationMultipleEvalIterativeSPCUDALinear::~OperationMultipleEvalIterativeSPCUDALinear() {
      deleteGridSPCUDA();
      deleteDataSPCUDA();
      delete myTimer;
    }

    void OperationMultipleEvalIterativeSPCUDALinear::rebuildLevelAndIndex(size_t gridFrom, size_t gridTo) {
      deleteGridSPCUDA();
      gridFrom = gridTo;
      gridTo = gridFrom;

      delete this->level_;
      delete this->index_;

      this->level_ = new sg::base::DataMatrixSP(this->storage_->size(), this->storage_->dim());
      this->index_ = new sg::base::DataMatrixSP(this->storage_->size(), this->storage_->dim());

      this->storage_->getLevelIndexArraysForEval(*(this->level_), *(this->index_));

      uploadGridSPCUDA(this->level_->getPointer(), this->index_->getPointer(), this->storage_->size(), this->storage_->dim());
    }
    
    double OperationMultipleEvalIterativeSPCUDALinear::multTransposeVectorized(sg::base::DataVectorSP& source, sg::base::DataVectorSP& result) {
      size_t source_size = source.getSize();
      size_t dims = this->storage_->dim();
      size_t storageSize = this->storage_->size();
      double time = 0.0;

      result.setAll(0.0f);

      float* ptrSource = source.getPointer();
      float* ptrGlobalResult = result.getPointer();

      //std::cout << this->dataset_->getNrows() << std::endl;
      //std::cout << this->dataset_->getNcols() << std::endl;

      if (this->dataset_->getNrows() % 64 != 0 || source_size != this->dataset_->getNrows()) {
        throw sg::base::operation_exception("For iterative mult transpose an even number of instances is required and result vector length must fit to data!");
      }

      size_t tmp = storageSize / 64;
      size_t myStorageSize = tmp * 64;

      if (myStorageSize > 0) {
        time = multTransSPCUDA(ptrSource, ptrGlobalResult, source_size, myStorageSize, dims);
      }

      // do padding
      float* ptrData = this->dataset_->getPointer();
      float* ptrLevel = this->level_->getPointer();
      float* ptrIndex = this->index_->getPointer();

      #pragma omp parallel for

      for (size_t j = myStorageSize; j < storageSize; j++) {
        ptrGlobalResult[j] = 0.0f;

        for (size_t i = 0; i < source_size; i++) {
          float curSupport = ptrSource[i];

          for (size_t d = 0; d < dims; d++) {
            float eval = ((ptrLevel[(j * dims) + d]) * (ptrData[(i * dims) + d]));
            float index_calc = eval - (ptrIndex[(j * dims) + d]);
            float abs = (float)fabs(index_calc);
            float last = 1.0f - abs;
            float localSupport = std::max<float>(last, 0.0f);
            curSupport *= localSupport;
          }

          ptrGlobalResult[j] += curSupport;
        }
      }

      return time;
    }
    
    double OperationMultipleEvalIterativeSPCUDALinear::multVectorized(sg::base::DataVectorSP& alpha, sg::base::DataVectorSP& result) {
      size_t result_size = result.getSize();
      size_t dims = this->storage_->dim();
      size_t storageSize = this->storage_->size();
      double time = 0.0;

      result.setAll(0.0f);

      float* ptrAlpha = alpha.getPointer();
      float* ptrResult = result.getPointer();

      //std::cout << this->dataset_->getNrows() << std::endl;
      //std::cout << this->dataset_->getNcols() << std::endl;
      
      if (this->dataset_->getNrows() % 64 != 0 || result_size != this->dataset_->getNrows()) {
        throw sg::base::operation_exception("For iterative mult an even number of instances is required and result vector length must fit to data!");
      }

      time = multSPCUDA(ptrAlpha, ptrResult, result_size, storageSize, dims);

      return time;
    }
  }

}
