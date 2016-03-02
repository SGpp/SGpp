// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/parallel/datadriven/basis/modlinear/operation/OperationMultipleEvalIterativeArBBModLinear.hpp>
#include <sgpp/base/exception/operation_exception.hpp>

#define ARBB_ARRAY

#include <sgpp/globaldef.hpp>

namespace SGPP {
namespace parallel {

OperationMultipleEvalIterativeArBBModLinear::OperationMultipleEvalIterativeArBBModLinear(
    SGPP::base::GridStorage* storage, SGPP::base::DataMatrix* dataset)
    : SGPP::parallel::OperationMultipleEvalVectorized(storage, dataset) {
  this->storage = storage;

  this->level_ = new SGPP::base::DataMatrix(storage->getSize(), storage->getDimension());
  this->index_ = new SGPP::base::DataMatrix(storage->getSize(), storage->getDimension());

  storage->getLevelIndexArraysForEval(*(this->level_), *(this->index_));

  myTimer = new SGPP::base::SGppStopwatch();
  myArBBKernels = new ArBBKernels();
  myArBBKernels2D = new ArBBKernels2D();
  //  myArBBKernels5D = new ArBBKernels5D();
}

OperationMultipleEvalIterativeArBBModLinear::~OperationMultipleEvalIterativeArBBModLinear() {
  delete myTimer;
  delete myArBBKernels;
  delete myArBBKernels2D;
  //  delete myArBBKernels5D;
}

void OperationMultipleEvalIterativeArBBModLinear::rebuildLevelAndIndex() {
  delete this->level_;
  delete this->index_;

  this->level_ = new SGPP::base::DataMatrix(storage->getSize(), storage->getDimension());
  this->index_ = new SGPP::base::DataMatrix(storage->getSize(), storage->getDimension());

  storage->getLevelIndexArraysForEval(*(this->level_), *(this->index_));

  myArBBKernels->resetKernels();
  myArBBKernels2D->resetKernels();
  //  myArBBKernels5D->resetKernels();
}

double OperationMultipleEvalIterativeArBBModLinear::multVectorized(SGPP::base::DataVector& alpha,
                                                                   SGPP::base::DataVector& result) {
  size_t result_size = result.getSize();
  size_t dims = storage->getDimension();
  size_t storageSize = storage->getSize();

  result.setAll(0.0);

  double* ptrAlpha = alpha.getPointer();
  double* ptrData = this->dataset_->getPointer();
  double* ptrLevel = this->level_->getPointer();
  double* ptrIndex = this->index_->getPointer();
  double* ptrResult = result.getPointer();

  if (this->dataset_->getNrows() % 16 != 0 || result_size != this->dataset_->getNrows()) {
    throw SGPP::base::operation_exception(
        "For iterative mult an even number of instances is required and result vector length must "
        "fit to data!");
  }

  double time = 0.0;

  if (this->dataset_->getNcols() == 2) {
#ifdef ARBB_ARRAY
    time = myArBBKernels2D->multModArBB(ptrAlpha, ptrData, ptrLevel, ptrIndex, ptrResult,
                                        result_size, storageSize, dims);
#else
    time = myArBBKernels->multModArBB(ptrAlpha, ptrData, ptrLevel, ptrIndex, ptrResult, result_size,
                                      storageSize, dims);
#endif
  }
  //    else if (this->dataset_->getNcols() == 5)
  //    {
  //#ifdef ARBB_ARRAY
  //      time = myArBBKernels5D->multModArBB(ptrAlpha, ptrData, ptrLevel, ptrIndex, ptrResult,
  //      result_size, storageSize, dims);
  //#else
  //      time = myArBBKernels->multModArBB(ptrAlpha, ptrData, ptrLevel, ptrIndex, ptrResult,
  //      result_size, storageSize, dims);
  //#endif
  //    }
  else {
    time = myArBBKernels->multModArBB(ptrAlpha, ptrData, ptrLevel, ptrIndex, ptrResult, result_size,
                                      storageSize, dims);
  }

  return time;
}

double OperationMultipleEvalIterativeArBBModLinear::multTransposeVectorized(
    SGPP::base::DataVector& source, SGPP::base::DataVector& result) {
  size_t soruceSize = source.getSize();
  size_t dims = storage->getDimension();
  size_t storageSize = storage->getSize();

  result.setAll(0.0);

  double* ptrSource = source.getPointer();
  double* ptrData = this->dataset_->getPointer();
  double* ptrLevel = this->level_->getPointer();
  double* ptrIndex = this->index_->getPointer();
  double* ptrGlobalResult = result.getPointer();

  if (this->dataset_->getNrows() % 16 != 0 || soruceSize != this->dataset_->getNrows()) {
    throw SGPP::base::operation_exception(
        "For iterative mult transpose an even number of instances is required and result vector "
        "length must fit to data!");
  }

  double time = 0.0;

  if (this->dataset_->getNcols() == 2) {
#ifdef ARBB_ARRAY
    time = myArBBKernels2D->multModTransArBB(ptrSource, ptrData, ptrLevel, ptrIndex,
                                             ptrGlobalResult, soruceSize, storageSize, dims);
#else
    time = myArBBKernels->multModTransArBB(ptrSource, ptrData, ptrLevel, ptrIndex, ptrGlobalResult,
                                           soruceSize, storageSize, dims);
#endif
  }
  //    else if (this->dataset_->getNcols() == 5)
  //    {
  //#ifdef ARBB_ARRAY
  //      time = myArBBKernels5D->multModTransArBB(ptrSource, ptrData, ptrLevel, ptrIndex,
  //      ptrGlobalResult, soruceSize, storageSize, dims);
  //#else
  //      time = myArBBKernels->multModTransArBB(ptrSource, ptrData, ptrLevel, ptrIndex,
  //      ptrGlobalResult, soruceSize, storageSize, dims);
  //#endif
  //    }
  else {
    time = myArBBKernels->multModTransArBB(ptrSource, ptrData, ptrLevel, ptrIndex, ptrGlobalResult,
                                           soruceSize, storageSize, dims);
  }

  return time;
}
}
}
