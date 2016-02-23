// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/grid/type/LinearGrid.hpp>
#include <sgpp/base/grid/generation/GridGenerator.hpp>
#include <sgpp/base/exception/operation_exception.hpp>

#include <sgpp/parallel/pde/basis/linear/boundary/operation/OperationLaplaceVectorizedLinearBoundaryOCL.hpp>

#include <cmath>
#include <assert.h>

#include <sgpp/globaldef.hpp>


namespace SGPP {
namespace parallel {

OperationLaplaceVectorizedLinearBoundaryOCL::OperationLaplaceVectorizedLinearBoundaryOCL(
  SGPP::base::GridStorage* storage,
  SGPP::base::DataVector& lambda) : storage(storage) {

  this->lambda = new SGPP::base::DataVector(lambda);
  this->OCLPDEKernelsHandle = OCLPDEKernels();
  this->level_ = new SGPP::base::DataMatrix(storage->getSize(), storage->getDimension());
  this->level_int_ = new SGPP::base::DataMatrix(storage->getSize(), storage->getDimension());
  this->index_ = new SGPP::base::DataMatrix(storage->getSize(), storage->getDimension());
  lcl_q = new double[this->storage->getDimension()];
  lcl_q_inv = new double[this->storage->getDimension()];
  storage->getLevelIndexArraysForEval(*(this->level_), *(this->index_));
  storage->getLevelForIntegral(*(this->level_int_));

}

OperationLaplaceVectorizedLinearBoundaryOCL::OperationLaplaceVectorizedLinearBoundaryOCL(
  SGPP::base::GridStorage* storage) : storage(storage) {


  this->lambda = new base::DataVector(storage->getDimension());
  this->lambda->setAll(1.0);
  this->OCLPDEKernelsHandle = OCLPDEKernels();
  this->level_ = new SGPP::base::DataMatrix(storage->getSize(), storage->getDimension());
  this->level_int_ = new SGPP::base::DataMatrix(storage->getSize(), storage->getDimension());
  this->index_ = new SGPP::base::DataMatrix(storage->getSize(), storage->getDimension());
  lcl_q = new double[this->storage->getDimension()];
  lcl_q_inv = new double[this->storage->getDimension()];
  storage->getLevelIndexArraysForEval(*(this->level_), *(this->index_));
  storage->getLevelForIntegral(*(this->level_int_));

}


OperationLaplaceVectorizedLinearBoundaryOCL::~OperationLaplaceVectorizedLinearBoundaryOCL() {
  delete this->level_;
  delete this->level_int_;
  delete this->index_;
  delete[] lcl_q;
  delete[] lcl_q_inv;

  this->OCLPDEKernelsHandle.CleanUpGPU();
}

void OperationLaplaceVectorizedLinearBoundaryOCL::mult_dirichlet(
  SGPP::base::DataVector& alpha, SGPP::base::DataVector& result) {
  result.setAll(0.0);

  this->OCLPDEKernelsHandle.RunOCLKernelLaplaceBound(alpha, result, lcl_q,
      lcl_q_inv,
      this->level_->getPointer(),
      this->index_->getPointer(),
      this->level_int_->getPointer(),
      lambda->getPointer(),
      storage->getSize(),
      storage->getDimension(),
      storage);


}

void OperationLaplaceVectorizedLinearBoundaryOCL::mult(SGPP::base::DataVector&
    alpha, SGPP::base::DataVector& result) {
  result.setAll(0.0);
  bool dirichlet = true;

  // fill q array
  for (size_t d = 0; d < this->storage->getDimension(); d++) {
    SGPP::base::BoundingBox* boundingBox = this->storage->getBoundingBox();
    lcl_q[d] = boundingBox->getIntervalWidth(d);
    lcl_q_inv[d] = 1.0 / boundingBox->getIntervalWidth(d);
    dirichlet = dirichlet && boundingBox->hasDirichletBoundaryLeft(d);
    dirichlet = dirichlet && boundingBox->hasDirichletBoundaryRight(d);
  }

  if (dirichlet) {
    mult_dirichlet(alpha, result);
  } else {
    throw SGPP::base::operation_exception("OperationLaplaceVectorizedLinearBoundaryOCL::mult : This method is only available on grids with Dirichlet boundaries in all dimensions!");
  }
}
}
}