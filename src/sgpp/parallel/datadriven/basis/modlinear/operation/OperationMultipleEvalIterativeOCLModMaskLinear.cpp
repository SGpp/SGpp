/******************************************************************************
* Copyright (C) 2013 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "parallel/datadriven/basis/modlinear/operation/OperationMultipleEvalIterativeOCLModMaskLinear.hpp"
#include "base/exception/operation_exception.hpp"

namespace sg {

  namespace parallel {

    OperationMultipleEvalIterativeOCLModMaskLinear::OperationMultipleEvalIterativeOCLModMaskLinear(
      base::GridStorage* storage, base::DataMatrix* dataset) :
      sg::parallel::OperationMultipleEvalVectorized(storage, dataset) {
      this->level_ = new base::DataMatrix(storage->size(), storage->dim());
      this->index_ = new base::DataMatrix(storage->size(), storage->dim());
      this->mask_ = new base::DataMatrix(storage->size(), storage->dim());
      this->offset_ = new base::DataMatrix(storage->size(), storage->dim());

      storage->getLevelIndexMaskArraysForModEval(*(this->level_), *(this->index_), *(this->mask_), *(this->offset_));

      myOCLKernels = new OCLKernels();
    }

    OperationMultipleEvalIterativeOCLModMaskLinear::~OperationMultipleEvalIterativeOCLModMaskLinear() {
      delete myOCLKernels;
    }

    void OperationMultipleEvalIterativeOCLModMaskLinear::rebuildLevelAndIndex() {
      delete this->level_;
      delete this->index_;
      delete this->mask_;
      delete this->offset_;

      this->level_ = new base::DataMatrix(storage_->size(), storage_->dim());
      this->index_ = new base::DataMatrix(storage_->size(), storage_->dim());
      this->mask_ = new base::DataMatrix(storage_->size(), storage_->dim());
      this->offset_ = new base::DataMatrix(storage_->size(), storage_->dim());

      storage_->getLevelIndexMaskArraysForModEval(*(this->level_), *(this->index_), *(this->mask_), *(this->offset_));

      myOCLKernels->resetKernels();
    }

    double OperationMultipleEvalIterativeOCLModMaskLinear::multTransposeVectorized(base::DataVector& source, base::DataVector& result) {
      size_t source_size = source.getSize();
      size_t dims = storage_->dim();
      size_t storageSize = storage_->size();

      result.setAll(0.0);

      double* ptrSource = source.getPointer();
      double* ptrData = this->dataset_->getPointer();
      double* ptrLevel = this->level_->getPointer();
      double* ptrIndex = this->index_->getPointer();
      double* ptrMask = this->mask_->getPointer();
      double* ptrOffset = this->offset_->getPointer();
      double* ptrGlobalResult = result.getPointer();

      if (this->dataset_->getNrows() % OCL_SGPP_LOCAL_WORKGROUP_SIZE != 0 || source_size != this->dataset_->getNrows()) {
        throw base::operation_exception("For iterative mult an even number of instances is required and result vector length must fit to data!");
      }

      size_t numWGs = storageSize / OCL_SGPP_LOCAL_WORKGROUP_SIZE;
      size_t global = numWGs * OCL_SGPP_LOCAL_WORKGROUP_SIZE;

      double time = 0.0;

      if (global > 0)
        time = myOCLKernels->multTransModMaskOCL(ptrSource, ptrData, ptrLevel, ptrIndex, ptrMask, ptrOffset, ptrGlobalResult, source_size, storageSize, dims, global);

      #pragma omp parallel for

      for (size_t j = global; j < storageSize; j++) {
        ptrGlobalResult[j] = 0.0f;

        for (size_t i = 0; i < source_size; i++) {
          double curSupport = ptrSource[i];

          for (size_t d = 0; d < dims; d++) {
            double eval = ((ptrLevel[(j * dims) + d]) * (ptrData[(i * dims) + d])) - (ptrIndex[(j * dims) + d]);
            uint64_t maskresult = *reinterpret_cast<uint64_t*>(&eval) | *reinterpret_cast<uint64_t*>(&(ptrMask[(j * dims) + d]));
            double masking = *reinterpret_cast<double*>( &maskresult );
            double last = masking + ptrOffset[(j * dims) + d];
            double localSupport = std::max<double>(last, 0.0);
            curSupport *= localSupport;
          }

          ptrGlobalResult[j] += curSupport;
        }
      }

      return time;
    }

    double OperationMultipleEvalIterativeOCLModMaskLinear::multVectorized(base::DataVector& alpha, base::DataVector& result) {
      size_t result_size = result.getSize();
      size_t dims = storage_->dim();
      size_t storageSize = storage_->size();

      result.setAll(0.0);

      double* ptrAlpha = alpha.getPointer();
      double* ptrData = this->dataset_->getPointer();
      double* ptrResult = result.getPointer();
      double* ptrLevel = this->level_->getPointer();
      double* ptrIndex = this->index_->getPointer();
      double* ptrMask = this->mask_->getPointer();
      double* ptrOffset = this->offset_->getPointer();

      if (this->dataset_->getNrows() % OCL_SGPP_LOCAL_WORKGROUP_SIZE != 0 || result_size != this->dataset_->getNrows()) {
        throw base::operation_exception("For iterative mult transpose an even number of instances is required and result vector length must fit to data!");
      }

      double time = myOCLKernels->multModMaskOCL(ptrAlpha, ptrData, ptrLevel, ptrIndex, ptrMask, ptrOffset, ptrResult, result_size, storageSize, dims, result_size);

      return time;
    }

  }

}
