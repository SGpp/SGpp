/******************************************************************************
* Copyright (C) 2010 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include <sgpp/parallel/datadriven/basis/linear/noboundary/operation/OperationMultipleEvalIterativeArBBLinear.hpp>
#include <sgpp/base/exception/operation_exception.hpp>

#define ARBB_ARRAY

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace parallel {

    OperationMultipleEvalIterativeArBBLinear::OperationMultipleEvalIterativeArBBLinear(SGPP::base::GridStorage* storage, SGPP::base::DataMatrix* dataset) : SGPP::parallel::OperationMultipleEvalVectorized(dataset) {
      this->storage = storage;

      this->level_ = new SGPP::base::DataMatrix(storage->size(), storage->dim());
      this->index_ = new SGPP::base::DataMatrix(storage->size(), storage->dim());

      storage->getLevelIndexArraysForEval(*(this->level_), *(this->index_));

      myTimer = new SGPP::base::SGppStopwatch();
      myArBBKernels = new ArBBKernels();
      myArBBKernels2D = new ArBBKernels2D();
      myArBBKernels4D = new ArBBKernels4D();
      myArBBKernels5D = new ArBBKernels5D();
      myArBBKernels10D = new ArBBKernels10D();
    }

    OperationMultipleEvalIterativeArBBLinear::~OperationMultipleEvalIterativeArBBLinear() {
      delete myTimer;
      delete myArBBKernels;
      delete myArBBKernels2D;
      delete myArBBKernels4D;
      delete myArBBKernels5D;
      delete myArBBKernels10D;
    }

    void OperationMultipleEvalIterativeArBBLinear::rebuildLevelAndIndex() {
      delete this->level_;
      delete this->index_;

      this->level_ = new SGPP::base::DataMatrix(storage->size(), storage->dim());
      this->index_ = new SGPP::base::DataMatrix(storage->size(), storage->dim());

      storage->getLevelIndexArraysForEval(*(this->level_), *(this->index_));

      myArBBKernels->resetKernels();
      myArBBKernels2D->resetKernels();
      myArBBKernels4D->resetKernels();
      myArBBKernels5D->resetKernels();
      myArBBKernels10D->resetKernels();
    }

    double OperationMultipleEvalIterativeArBBLinear::multVectorized(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result) {
      size_t result_size = result.getSize();
      size_t dims = storage->dim();
      size_t storageSize = storage->size();

      result.setAll(0.0);

      double* ptrAlpha = alpha.getPointer();
      double* ptrData =  this->dataset_->getPointer();
      double* ptrLevel = this->level_->getPointer();
      double* ptrIndex = this->index_->getPointer();
      double* ptrResult = result.getPointer();

      if (this->dataset_->getNrows() % 16 != 0 || result_size !=  this->dataset_->getNrows()) {
        throw SGPP::base::operation_exception("For iterative mult an even number of instances is required and result vector length must fit to data!");
      }

      double time = 0.0;

      if (this->dataset_->getNcols() == 2) {
#ifdef ARBB_ARRAY
        time = myArBBKernels2D->multArBB(ptrAlpha, ptrData, ptrLevel, ptrIndex, ptrResult, result_size, storageSize, dims);
#else
        time = myArBBKernels->multArBB(ptrAlpha, ptrData, ptrLevel, ptrIndex, ptrResult, result_size, storageSize, dims);
#endif
      } else if (this->dataset_->getNcols() == 4) {
#ifdef ARBB_ARRAY
        time = myArBBKernels4D->multArBB(ptrAlpha, ptrData, ptrLevel, ptrIndex, ptrResult, result_size, storageSize, dims);
#else
        time = myArBBKernels->multArBB(ptrAlpha, ptrData, ptrLevel, ptrIndex, ptrResult, result_size, storageSize, dims);
#endif
      } else if (this->dataset_->getNcols() == 5) {
#ifdef ARBB_ARRAY
        time = myArBBKernels5D->multArBB(ptrAlpha, ptrData, ptrLevel, ptrIndex, ptrResult, result_size, storageSize, dims);
#else
        time = myArBBKernels->multArBB(ptrAlpha, ptrData, ptrLevel, ptrIndex, ptrResult, result_size, storageSize, dims);
#endif
      } else if (this->dataset_->getNcols() == 10) {
#ifdef ARBB_ARRAY
        time = myArBBKernels10D->multArBB(ptrAlpha, ptrData, ptrLevel, ptrIndex, ptrResult, result_size, storageSize, dims);
#else
        time = myArBBKernels->multArBB(ptrAlpha, ptrData, ptrLevel, ptrIndex, ptrResult, result_size, storageSize, dims);
#endif
      } else {
        time = myArBBKernels->multArBB(ptrAlpha, ptrData, ptrLevel, ptrIndex, ptrResult, result_size, storageSize, dims);
      }

      return time;
    }

    double OperationMultipleEvalIterativeArBBLinear::multTransposeVectorized(SGPP::base::DataVector& source, SGPP::base::DataVector& result) {
      size_t soruceSize = source.getSize();
      size_t dims = storage->dim();
      size_t storageSize = storage->size();

      result.setAll(0.0);

      double* ptrSource = source.getPointer();
      double* ptrData = this->dataset_->getPointer();
      double* ptrLevel = this->level_->getPointer();
      double* ptrIndex = this->index_->getPointer();
      double* ptrGlobalResult = result.getPointer();

      if (this->dataset_->getNrows() % 16 != 0 || soruceSize != this->dataset_->getNrows()) {
        throw SGPP::base::operation_exception("For iterative mult transpose an even number of instances is required and result vector length must fit to data!");
      }

      double time = 0.0;

      if (this->dataset_->getNcols() == 2) {
#ifdef ARBB_ARRAY
        time = myArBBKernels2D->multTransArBB(ptrSource, ptrData, ptrLevel, ptrIndex, ptrGlobalResult, soruceSize, storageSize, dims);
#else
        time = myArBBKernels->multTransArBB(ptrSource, ptrData, ptrLevel, ptrIndex, ptrGlobalResult, soruceSize, storageSize, dims);
#endif
      } else if (this->dataset_->getNcols() == 4) {
#ifdef ARBB_ARRAY
        time = myArBBKernels4D->multTransArBB(ptrSource, ptrData, ptrLevel, ptrIndex, ptrGlobalResult, soruceSize, storageSize, dims);
#else
        time = myArBBKernels->multTransArBB(ptrSource, ptrData, ptrLevel, ptrIndex, ptrGlobalResult, soruceSize, storageSize, dims);
#endif
      } else if (this->dataset_->getNcols() == 5) {
#ifdef ARBB_ARRAY
        time = myArBBKernels5D->multTransArBB(ptrSource, ptrData, ptrLevel, ptrIndex, ptrGlobalResult, soruceSize, storageSize, dims);
#else
        time = myArBBKernels->multTransArBB(ptrSource, ptrData, ptrLevel, ptrIndex, ptrGlobalResult, soruceSize, storageSize, dims);
#endif
      } else if (this->dataset_->getNcols() == 10) {
#ifdef ARBB_ARRAY
        time = myArBBKernels10D->multTransArBB(ptrSource, ptrData, ptrLevel, ptrIndex, ptrGlobalResult, soruceSize, storageSize, dims);
#else
        time = myArBBKernels->multTransArBB(ptrSource, ptrData, ptrLevel, ptrIndex, ptrGlobalResult, soruceSize, storageSize, dims);
#endif
      } else {
        time = myArBBKernels->multTransArBB(ptrSource, ptrData, ptrLevel, ptrIndex, ptrGlobalResult, soruceSize, storageSize, dims);
      }

      return time;
    }

  }

}
