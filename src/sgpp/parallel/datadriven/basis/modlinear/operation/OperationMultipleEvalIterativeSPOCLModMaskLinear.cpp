/******************************************************************************
* Copyright (C) 2013 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "parallel/datadriven/basis/modlinear/operation/OperationMultipleEvalIterativeSPOCLModMaskLinear.hpp"
#include "base/exception/operation_exception.hpp"

namespace sg {

  namespace parallel {

    OperationMultipleEvalIterativeSPOCLModMaskLinear::OperationMultipleEvalIterativeSPOCLModMaskLinear(
      base::GridStorage* storage, base::DataMatrixSP* dataset) :
      sg::parallel::OperationMultipleEvalVectorizedSP(storage, dataset) {
      this->level_ = new base::DataMatrixSP(storage->size(), storage->dim());
      this->index_ = new base::DataMatrixSP(storage->size(), storage->dim());
      this->mask_ = new base::DataMatrixSP(storage->size(), storage->dim());
      this->offset_ = new base::DataMatrixSP(storage->size(), storage->dim());

      storage->getLevelIndexMaskArraysForModEval(*(this->level_), *(this->index_), *(this->mask_), *(this->offset_));

      myOCLKernels = new OCLKernels();
    }

    OperationMultipleEvalIterativeSPOCLModMaskLinear::~OperationMultipleEvalIterativeSPOCLModMaskLinear() {
      delete myOCLKernels;
    }

    void OperationMultipleEvalIterativeSPOCLModMaskLinear::rebuildLevelAndIndex() {
      delete this->level_;
      delete this->index_;
      delete this->mask_;
      delete this->offset_;

      this->level_ = new base::DataMatrixSP(storage_->size(), storage_->dim());
      this->index_ = new base::DataMatrixSP(storage_->size(), storage_->dim());
      this->mask_ = new base::DataMatrixSP(storage_->size(), storage_->dim());
      this->offset_ = new base::DataMatrixSP(storage_->size(), storage_->dim());

      storage_->getLevelIndexMaskArraysForModEval(*(this->level_), *(this->index_), *(this->mask_), *(this->offset_));

      myOCLKernels->resetKernels();
    }

    double OperationMultipleEvalIterativeSPOCLModMaskLinear::multTransposeVectorized(base::DataVectorSP& source, base::DataVectorSP& result) {
      size_t source_size = source.getSize();
      size_t dims = storage_->dim();
      size_t storageSize = storage_->size();

      result.setAll(0.0f);

      float* ptrSource = source.getPointer();
      float* ptrData = this->dataset_->getPointer();
      float* ptrLevel = this->level_->getPointer();
      float* ptrIndex = this->index_->getPointer();
      float* ptrMask = this->mask_->getPointer();
      float* ptrOffset = this->offset_->getPointer();
      float* ptrGlobalResult = result.getPointer();

      if (this->dataset_->getNcols() % OCL_SGPP_LOCAL_WORKGROUP_SIZE != 0 || source_size != this->dataset_->getNcols()) {
        throw base::operation_exception("For iterative mult an even number of instances is required and result vector length must fit to data!");
      }

      size_t numWGs = storageSize / OCL_SGPP_LOCAL_WORKGROUP_SIZE;
      size_t global = numWGs * OCL_SGPP_LOCAL_WORKGROUP_SIZE;

      double time = 0.0;

      if (global > 0)
        time = myOCLKernels->multTransModMaskSPOCL(ptrSource, ptrData, ptrLevel, ptrIndex, ptrMask, ptrOffset, ptrGlobalResult, source_size, storageSize, dims, global);

      #pragma omp parallel for

      for (size_t j = global; j < storageSize; j++) {
        ptrGlobalResult[j] = 0.0f;

        for (size_t i = 0; i < source_size; i++) {
          float curSupport = ptrSource[i];

          for (size_t d = 0; d < dims; d++) {
            float eval = ((ptrLevel[(j * dims) + d]) * (ptrData[(d * source_size) + i])) - (ptrIndex[(j * dims) + d]);
            uint64_t maskresult = *reinterpret_cast<unsigned int*>(&eval) | *reinterpret_cast<unsigned int*>(&(ptrMask[(j * dims) + d]));
            float masking = *reinterpret_cast<float*>( &maskresult );
            float last = masking + ptrOffset[(j * dims) + d];
            float localSupport = std::max<float>(last, 0.0f);
            curSupport *= localSupport;
          }

          ptrGlobalResult[j] += curSupport;
        }
      }

      return time;
    }

    double OperationMultipleEvalIterativeSPOCLModMaskLinear::multVectorized(base::DataVectorSP& alpha, base::DataVectorSP& result) {
      size_t result_size = result.getSize();
      size_t dims = storage_->dim();
      size_t storageSize = storage_->size();

      result.setAll(0.0f);

      float* ptrAlpha = alpha.getPointer();
      float* ptrData = this->dataset_->getPointer();
      float* ptrResult = result.getPointer();
      float* ptrLevel = this->level_->getPointer();
      float* ptrIndex = this->index_->getPointer();
      float* ptrMask = this->mask_->getPointer();
      float* ptrOffset = this->offset_->getPointer();

      if (this->dataset_->getNcols() % OCL_SGPP_LOCAL_WORKGROUP_SIZE != 0 || result_size != this->dataset_->getNcols()) {
        throw base::operation_exception("For iterative mult transpose an even number of instances is required and result vector length must fit to data!");
      }

      double time = myOCLKernels->multModMaskSPOCL(ptrAlpha, ptrData, ptrLevel, ptrIndex, ptrMask, ptrOffset, ptrResult, result_size, storageSize, dims, result_size);

      return time;
    }

  }

}
