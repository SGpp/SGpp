// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/parallel/datadriven/basis/modlinear/operation/OperationMultipleEvalIterativeSPArBBModLinear.hpp>
#include <sgpp/base/exception/operation_exception.hpp>

//#define ARBB_ARRAY

#include <sgpp/globaldef.hpp>

#if USE_DOUBLE_PRECISION==0

namespace SGPP {
  namespace parallel {

    OperationMultipleEvalIterativeSPArBBModLinear::OperationMultipleEvalIterativeSPArBBModLinear(SGPP::base::GridStorage* storage, SGPP::base::DataMatrixSP* dataset) : SGPP::parallel::OperationMultipleEvalVectorizedSP(dataset) {
      this->storage = storage;

      this->level_ = new SGPP::base::DataMatrixSP(storage->size(), storage->dim());
      this->index_ = new SGPP::base::DataMatrixSP(storage->size(), storage->dim());

      storage->getLevelIndexArraysForEval(*(this->level_), *(this->index_));

      myTimer = new SGPP::base::SGppStopwatch();
      myArBBKernels = new ArBBKernels();
      myArBBKernels2D = new ArBBKernels2D();
      //  myArBBKernels5D = new ArBBKernels5D();
    }

    OperationMultipleEvalIterativeSPArBBModLinear::~OperationMultipleEvalIterativeSPArBBModLinear() {
      delete myTimer;
      delete myArBBKernels;
      delete myArBBKernels2D;
      //  delete myArBBKernels5D;
    }

    void OperationMultipleEvalIterativeSPArBBModLinear::rebuildLevelAndIndex() {
      delete this->level_;
      delete this->index_;

      this->level_ = new SGPP::base::DataMatrixSP(storage->size(), storage->dim());
      this->index_ = new SGPP::base::DataMatrixSP(storage->size(), storage->dim());

      storage->getLevelIndexArraysForEval(*(this->level_), *(this->index_));

      myArBBKernels->resetKernels();
      myArBBKernels2D->resetKernels();
      //  myArBBKernels5D->resetKernels();
    }

    double OperationMultipleEvalIterativeSPArBBModLinear::multVectorized(SGPP::base::DataVectorSP& alpha, SGPP::base::DataVectorSP& result) {
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
        throw SGPP::base::operation_exception("For iterative mult transpose an even number of instances is required and result vector length must fit to data!");
      }

      double time = 0.0;

      if (this->dataset_->getNcols() == 2) {
#ifdef ARBB_ARRAY
        time = myArBBKernels2D->multModSPArBB(ptrAlpha, ptrData, ptrLevel, ptrIndex, ptrResult, result_size, storageSize, dims);
#else
        time = myArBBKernels->multModSPArBB(ptrAlpha, ptrData, ptrLevel, ptrIndex, ptrResult, result_size, storageSize, dims);
#endif
      }
      //    else if (this->dataset_->getNcols() == 5)
      //    {
      //#ifdef ARBB_ARRAY
      //      time = myArBBKernels5D->multModSPArBB(ptrAlpha, ptrData, ptrLevel, ptrIndex, ptrResult, result_size, storageSize, dims);
      //#else
      //      time = myArBBKernels->multModSPArBB(ptrAlpha, ptrData, ptrLevel, ptrIndex, ptrResult, result_size, storageSize, dims);
      //#endif
      //    }
      else {
        time = myArBBKernels->multModSPArBB(ptrAlpha, ptrData, ptrLevel, ptrIndex, ptrResult, result_size, storageSize, dims);
      }

      return time;
    }

    double OperationMultipleEvalIterativeSPArBBModLinear::multTransposeVectorized(SGPP::base::DataVectorSP& source, SGPP::base::DataVectorSP& result) {
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
        throw SGPP::base::operation_exception("For iterative mult an even number of instances is required and result vector length must fit to data!");
      }

      double time = 0.0;

      if (this->dataset_->getNcols() == 2) {
#ifdef ARBB_ARRAY
        time = myArBBKernels2D->multModTransSPArBB(ptrSource, ptrData, ptrLevel, ptrIndex, ptrGlobalResult, source_size, storageSize, dims);
#else
        time = myArBBKernels->multModTransSPArBB(ptrSource, ptrData, ptrLevel, ptrIndex, ptrGlobalResult, source_size, storageSize, dims);
#endif
      }
      //    else if (this->dataset_->getNcols() == 5)
      //    {
      //#ifdef ARBB_ARRAY
      //      time = myArBBKernels5D->multModTransSPArBB(ptrSource, ptrData, ptrLevel, ptrIndex, ptrGlobalResult, source_size, storageSize, dims);
      //#else
      //      time = myArBBKernels->multModTransSPArBB(ptrSource, ptrData, ptrLevel, ptrIndex, ptrGlobalResult, source_size, storageSize, dims);
      //#endif
      //    }
      else {
        time = myArBBKernels->multModTransSPArBB(ptrSource, ptrData, ptrLevel, ptrIndex, ptrGlobalResult, source_size, storageSize, dims);
      }

      return time;
    }

  }

}

#endif
