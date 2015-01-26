/******************************************************************************
* Copyright (C) 2010 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "parallel/datadriven/basis/linear/noboundary/operation/OperationMultipleEvalIterativeSPArBBLinear.hpp"
#include "base/exception/operation_exception.hpp"

#define ARBB_ARRAY

namespace sg {
  namespace parallel {

    OperationMultipleEvalIterativeSPArBBLinear::OperationMultipleEvalIterativeSPArBBLinear(sg::base::GridStorage* storage, sg::base::DataMatrixSP* dataset) : sg::parallel::OperationMultipleEvalVectorizedSP(dataset) {
      this->storage = storage;

      this->level_ = new sg::base::DataMatrixSP(storage->size(), storage->dim());
      this->index_ = new sg::base::DataMatrixSP(storage->size(), storage->dim());

      storage->getLevelIndexArraysForEval(*(this->level_), *(this->index_));

      myTimer = new sg::base::SGppStopwatch();
      myArBBKernels = new ArBBKernels();
      myArBBKernels2D = new ArBBKernels2D();
      myArBBKernels4D = new ArBBKernels4D();
      myArBBKernels5D = new ArBBKernels5D();
      myArBBKernels10D = new ArBBKernels10D();
    }

    OperationMultipleEvalIterativeSPArBBLinear::~OperationMultipleEvalIterativeSPArBBLinear() {
      delete myTimer;
      delete myArBBKernels;
      delete myArBBKernels2D;
      delete myArBBKernels4D;
      delete myArBBKernels5D;
      delete myArBBKernels10D;
    }

    void OperationMultipleEvalIterativeSPArBBLinear::rebuildLevelAndIndex() {
      delete this->level_;
      delete this->index_;

      this->level_ = new sg::base::DataMatrixSP(storage->size(), storage->dim());
      this->index_ = new sg::base::DataMatrixSP(storage->size(), storage->dim());

      storage->getLevelIndexArraysForEval(*(this->level_), *(this->index_));

      myArBBKernels->resetKernels();
      myArBBKernels2D->resetKernels();
      myArBBKernels4D->resetKernels();
      myArBBKernels5D->resetKernels();
      myArBBKernels10D->resetKernels();
    }

    double OperationMultipleEvalIterativeSPArBBLinear::multVectorized(sg::base::DataVectorSP& alpha, sg::base::DataVectorSP& result) {
      size_t result_size = result.getSize();
      size_t dims = storage->dim();
      size_t storageSize = storage->size();

      result.setAll(0.0f);

      float* ptrAlpha = alpha.getPointer();
      float* ptrData = this->dataset_->getPointer();
      float* ptrResult = result.getPointer();
      float* ptrLevel = this->level_->getPointer();
      float* ptrIndex = this->index_->getPointer();

      if (this->dataset_->getNrows() % 16 != 0 || result_size != this->dataset_->getNrows()) {
        throw sg::base::operation_exception("For iterative mult transpose an even number of instances is required and result vector length must fit to data!");
      }

      double time = 0.0;

      if (this->dataset_->getNcols() == 2) {
#ifdef ARBB_ARRAY
        time = myArBBKernels2D->multSPArBB(ptrAlpha, ptrData, ptrLevel, ptrIndex, ptrResult, result_size, storageSize, dims);
#else
        time = myArBBKernels->multSPArBB(ptrAlpha, ptrData, ptrLevel, ptrIndex, ptrResult, result_size, storageSize, dims);
#endif
      } else if (this->dataset_->getNcols() == 4) {
#ifdef ARBB_ARRAY
        time = myArBBKernels4D->multSPArBB(ptrAlpha, ptrData, ptrLevel, ptrIndex, ptrResult, result_size, storageSize, dims);
#else
        time = myArBBKernels->multSPArBB(ptrAlpha, ptrData, ptrLevel, ptrIndex, ptrResult, result_size, storageSize, dims);
#endif
      } else if (this->dataset_->getNcols() == 5) {
#ifdef ARBB_ARRAY
        time = myArBBKernels5D->multSPArBB(ptrAlpha, ptrData, ptrLevel, ptrIndex, ptrResult, result_size, storageSize, dims);
#else
        time = myArBBKernels->multSPArBB(ptrAlpha, ptrData, ptrLevel, ptrIndex, ptrResult, result_size, storageSize, dims);
#endif
      } else if (this->dataset_->getNcols() == 10) {
#ifdef ARBB_ARRAY
        time = myArBBKernels10D->multSPArBB(ptrAlpha, ptrData, ptrLevel, ptrIndex, ptrResult, result_size, storageSize, dims);
#else
        time = myArBBKernels->multSPArBB(ptrAlpha, ptrData, ptrLevel, ptrIndex, ptrResult, result_size, storageSize, dims);
#endif
      } else {
        time = myArBBKernels->multSPArBB(ptrAlpha, ptrData, ptrLevel, ptrIndex, ptrResult, result_size, storageSize, dims);
      }

      return time;
    }

    double OperationMultipleEvalIterativeSPArBBLinear::multTransposeVectorized(sg::base::DataVectorSP& source, sg::base::DataVectorSP& result) {
      size_t source_size = source.getSize();
      size_t dims = storage->dim();
      size_t storageSize = storage->size();

      result.setAll(0.0f);

      float* ptrSource = source.getPointer();
      float* ptrData = this->dataset_->getPointer();
      float* ptrLevel = this->level_->getPointer();
      float* ptrIndex = this->index_->getPointer();
      float* ptrGlobalResult = result.getPointer();

      if (this->dataset_->getNrows() % 16 != 0 || source_size != this->dataset_->getNrows()) {
        throw sg::base::operation_exception("For iterative mult an even number of instances is required and result vector length must fit to data!");
      }

      double time = 0.0;

      if (this->dataset_->getNcols() == 2) {
#ifdef ARBB_ARRAY
        time = myArBBKernels2D->multTransSPArBB(ptrSource, ptrData, ptrLevel, ptrIndex, ptrGlobalResult, source_size, storageSize, dims);
#else
        time = myArBBKernels->multTransSPArBB(ptrSource, ptrData, ptrLevel, ptrIndex, ptrGlobalResult, source_size, storageSize, dims);
#endif
      } else if (this->dataset_->getNcols() == 4) {
#ifdef ARBB_ARRAY
        time = myArBBKernels4D->multTransSPArBB(ptrSource, ptrData, ptrLevel, ptrIndex, ptrGlobalResult, source_size, storageSize, dims);
#else
        time = myArBBKernels->multTransSPArBB(ptrSource, ptrData, ptrLevel, ptrIndex, ptrGlobalResult, source_size, storageSize, dims);
#endif
      } else if (this->dataset_->getNcols() == 5) {
#ifdef ARBB_ARRAY
        time = myArBBKernels5D->multTransSPArBB(ptrSource, ptrData, ptrLevel, ptrIndex, ptrGlobalResult, source_size, storageSize, dims);
#else
        time = myArBBKernels->multTransSPArBB(ptrSource, ptrData, ptrLevel, ptrIndex, ptrGlobalResult, source_size, storageSize, dims);
#endif
      } else if (this->dataset_->getNcols() == 10) {
#ifdef ARBB_ARRAY
        time = myArBBKernels10D->multTransSPArBB(ptrSource, ptrData, ptrLevel, ptrIndex, ptrGlobalResult, source_size, storageSize, dims);
#else
        time = myArBBKernels->multTransSPArBB(ptrSource, ptrData, ptrLevel, ptrIndex, ptrGlobalResult, source_size, storageSize, dims);
#endif
      } else {
        time = myArBBKernels->multTransSPArBB(ptrSource, ptrData, ptrLevel, ptrIndex, ptrGlobalResult, source_size, storageSize, dims);
      }

      return time;
    }

  }

}
