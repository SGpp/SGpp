// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/exception/operation_exception.hpp>
#include <sgpp/datadriven/DatadrivenOpFactory.hpp>
#include <sgpp/base/operation/BaseOpFactory.hpp>
#include <sgpp/datadriven/algorithm/SystemMatrixDensityDerivativeRatioEstimation.hpp>
#include <sgpp/globaldef.hpp>

// #include <iostream>

// using namespace std;

namespace sgpp {
namespace datadriven {

SystemMatrixDensityDerivativeRatioEstimation::SystemMatrixDensityDerivativeRatioEstimation(
    base::Grid& grid, base::DataMatrix& trainData, double lambda, size_t derivDim)
    : DMSystemMatrixBase(trainData, lambda),
      instances(0),
      paddedInstances(0),
      grid(grid),
      derivDim(derivDim) {
  this->instances = this->dataset_.getNrows();
  this->B.reset(op_factory::createOperationMultipleEval(grid, this->dataset_,
                                                        this->implementationConfiguration));
  this->Bderiv.reset(op_factory::createOperationMultipleEvalPartialDerivativeNaive(
      grid, this->dataset_, derivDim));
  // padded during Operator construction, fetch new size
  this->paddedInstances = this->dataset_.getNrows();
}

SystemMatrixDensityDerivativeRatioEstimation::~SystemMatrixDensityDerivativeRatioEstimation() {}

void SystemMatrixDensityDerivativeRatioEstimation::mult(base::DataVector& alpha,
                                                        base::DataVector& result) {
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

void SystemMatrixDensityDerivativeRatioEstimation::generateb(base::DataVector& classes,
                                                             base::DataVector& b) {
  // classes are irrelevant;
  // Compute d_B * (-1)
  base::DataVector y(this->instances, -1.0);
  this->myTimer_->start();
  this->Bderiv->multTranspose(y, b);
  this->completeTimeMultTrans_ += this->myTimer_->stop();
  this->computeTimeMultTrans_ += this->Bderiv->getDuration();
}

void SystemMatrixDensityDerivativeRatioEstimation::prepareGrid() {
  this->B->prepare();
  this->Bderiv->prepare();
}

void SystemMatrixDensityDerivativeRatioEstimation::setImplementation(
    datadriven::OperationMultipleEvalConfiguration operationConfiguration) {
  this->implementationConfiguration = operationConfiguration;
  this->B.reset(op_factory::createOperationMultipleEval(this->grid, this->dataset_,
                                                        this->implementationConfiguration));
  this->Bderiv.reset(op_factory::createOperationMultipleEvalPartialDerivativeNaive(
      grid, this->dataset_, derivDim));
  // padded during Operator construction, fetch new size
  this->paddedInstances = this->dataset_.getNrows();
}

}  // namespace datadriven
}  // namespace sgpp
