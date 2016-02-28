// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include "sgpp/base/exception/operation_exception.hpp"
#include "sgpp/datadriven/DatadrivenOpFactory.hpp"
#include "SystemMatrixLeastSquaresIdentity.hpp"
#include "sgpp/globaldef.hpp"

// #include <iostream>

// using namespace std;

namespace SGPP {
namespace datadriven {

SystemMatrixLeastSquaresIdentity::SystemMatrixLeastSquaresIdentity(base::Grid& grid,
                                                                   base::DataMatrix& trainData,
                                                                   float_t lambda)
    : DMSystemMatrixBase(trainData, lambda), instances(0), paddedInstances(0), grid(grid) {
  this->dataset_ = new base::DataMatrix(trainData);
  this->instances = this->dataset_->getNrows();
  // this->paddedInstances = PaddingAssistant::padDataset(*(this->dataset_));
  // datadriven::OperationMultipleEvalType type =
  // datadriven::OperationMultipleEvalType::SUBSPACELINEAR;
  this->B = op_factory::createOperationMultipleEval(grid, *(this->dataset_),
                                                    this->implementationConfiguration).release();
  // padded during Operator construction, fetch new size
  this->paddedInstances = this->dataset_->getNrows();
}

SystemMatrixLeastSquaresIdentity::~SystemMatrixLeastSquaresIdentity() {
  delete this->B;
  delete this->dataset_;
}

void SystemMatrixLeastSquaresIdentity::mult(base::DataVector& alpha, base::DataVector& result) {
  base::DataVector temp(this->paddedInstances);

  // Operation B
  this->myTimer_->start();
  this->B->mult(alpha, temp);
  this->completeTimeMult_ += this->myTimer_->stop();

  this->myTimer_->start();
  this->B->multTranspose(temp, result);
  this->completeTimeMultTrans_ += this->myTimer_->stop();

  result.axpy(static_cast<float_t>(this->instances) * this->lambda_, alpha);
}

void SystemMatrixLeastSquaresIdentity::generateb(base::DataVector& classes, base::DataVector& b) {
  base::DataVector myClasses(classes);

  this->myTimer_->start();
  this->B->multTranspose(myClasses, b);
  this->completeTimeMultTrans_ += this->myTimer_->stop();
}

void SystemMatrixLeastSquaresIdentity::prepareGrid() { this->B->prepare(); }

}  // namespace datadriven
}  // namespace SGPP
