// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/grid/type/LinearGrid.hpp>
#include <sgpp/base/grid/generation/GridGenerator.hpp>

#include <sgpp/parallel/pde/basis/linear/noboundary/operation/OperationLTwoDotProductVectorizedLinearOCL.hpp>

#include <cmath>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace parallel {

OperationLTwoDotProductVectorizedLinearOCL::OperationLTwoDotProductVectorizedLinearOCL(
    sgpp::base::GridStorage* storage)
    : storage(storage) {
  this->lcl_q = new double[this->storage->getDimension()];
  this->OCLPDEKernelsHandle = OCLPDEKernels();

  this->level_ = new sgpp::base::DataMatrix(storage->getSize(), storage->getDimension());
  this->level_int_ = new sgpp::base::DataMatrix(storage->getSize(), storage->getDimension());
  this->index_ = new sgpp::base::DataMatrix(storage->getSize(), storage->getDimension());

  storage->getLevelIndexArraysForEval(*(this->level_), *(this->index_));
  storage->getLevelForIntegral(*(this->level_int_));
}

OperationLTwoDotProductVectorizedLinearOCL::~OperationLTwoDotProductVectorizedLinearOCL() {
  delete this->level_;
  delete this->level_int_;
  delete this->index_;
  delete[] lcl_q;
  this->OCLPDEKernelsHandle.CleanUpGPU();
}

void OperationLTwoDotProductVectorizedLinearOCL::mult(sgpp::base::DataVector& alpha,
                                                      sgpp::base::DataVector& result) {
  result.setAll(0.0);

  for (size_t d = 0; d < this->storage->getDimension(); d++) {
    sgpp::base::BoundingBox* boundingBox = this->storage->getBoundingBox();
    this->lcl_q[d] = boundingBox->getIntervalWidth(d);
  }

  this->OCLPDEKernelsHandle.RunOCLKernelLTwoDotInner(
      alpha, result, this->lcl_q, this->level_->getPointer(), this->index_->getPointer(),
      this->level_int_->getPointer(), storage->getSize(), storage->getDimension(), storage);
}
}
}