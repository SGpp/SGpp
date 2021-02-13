// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/exception/operation_exception.hpp>
#include <sgpp/datadriven/DatadrivenOpFactory.hpp>
#include <sgpp/datadriven/algorithm/SystemMatrixLeastSquaresIdentity.hpp>
#include <sgpp/globaldef.hpp>

// #include <iostream>

// using namespace std;

namespace sgpp {
namespace datadriven {

SystemMatrixLeastSquaresIdentity::SystemMatrixLeastSquaresIdentity(base::Grid& grid,
                                                                   base::DataMatrix& trainData,
                                                                   double lambda)
    : DMSystemMatrixBase(trainData, lambda), instances(0), paddedInstances(0), grid(grid) {
  this->instances = this->dataset_.getNrows();
  this->B.reset(op_factory::createOperationMultipleEval(grid, this->dataset_,
                                                        this->implementationConfiguration));

  // padded during Operator construction, fetch new size
  this->paddedInstances = this->dataset_.getNrows();
}

SystemMatrixLeastSquaresIdentity::~SystemMatrixLeastSquaresIdentity() {}

void SystemMatrixLeastSquaresIdentity::mult(base::DataVector& alpha, base::DataVector& result) {
  base::DataVector temp(this->paddedInstances);

  // Operation B
  this->myTimer_->start();
  this->B->mult(alpha, temp);
  this->completeTimeMult_ += this->myTimer_->stop();
  this->computeTimeMult_ += this->B->getDuration();

  this->myTimer_->start();
  this->B->multTranspose(temp, result);
  this->completeTimeMultTrans_ += this->myTimer_->stop();
  this->computeTimeMultTrans_ += this->B->getDuration();

  result.axpy(static_cast<double>(this->instances) * this->lambda_, alpha);
}

void SystemMatrixLeastSquaresIdentity::generateb(base::DataVector& classes, base::DataVector& b) {
  base::DataVector myClasses(classes);

  this->myTimer_->start();
  this->B->multTranspose(myClasses, b);
  this->completeTimeMultTrans_ += this->myTimer_->stop();
  this->computeTimeMultTrans_ += this->B->getDuration();
}

void SystemMatrixLeastSquaresIdentity::prepareGrid() { this->B->prepare(); }

void SystemMatrixLeastSquaresIdentity::setImplementation(
    datadriven::OperationMultipleEvalConfiguration operationConfiguration) {
  this->implementationConfiguration = operationConfiguration;
  this->B.reset(op_factory::createOperationMultipleEval(this->grid, this->dataset_,
                                                        this->implementationConfiguration));
  // padded during Operator construction, fetch new size
  this->paddedInstances = this->dataset_.getNrows();
}

}  // namespace datadriven
}  // namespace sgpp
