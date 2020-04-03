// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/datadriven/algorithm/SystemMatrixDensityRatioEstimation.hpp>

#include <sgpp/base/exception/operation_exception.hpp>
#include <sgpp/datadriven/DatadrivenOpFactory.hpp>
#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace datadriven {

SystemMatrixDensityRatioEstimation::SystemMatrixDensityRatioEstimation(base::Grid& grid,
                                                                       base::DataMatrix& trainDataP,
                                                                       base::DataMatrix& trainDataQ,
                                                                       double lambda)
    : DMSystemMatrixDRE(trainDataP, trainDataQ, lambda),
      instancesP(0),
      paddedInstancesP(0),
      instancesQ(0),
      paddedInstancesQ(0),
      grid(grid) {
  this->instancesP = this->datasetP_.getNrows();
  this->instancesQ = this->datasetQ_.getNrows();
  this->B_p.reset(op_factory::createOperationMultipleEval(grid, this->datasetP_,
                                                          this->implementationConfiguration));
  this->B_q.reset(op_factory::createOperationMultipleEval(grid, this->datasetQ_,
                                                          this->implementationConfiguration));

  // padded during Operator construction, fetch new size
  this->paddedInstancesP = this->datasetP_.getNrows();
  // padded during Operator construction, fetch new size
  this->paddedInstancesQ = this->datasetQ_.getNrows();
}

SystemMatrixDensityRatioEstimation::~SystemMatrixDensityRatioEstimation() {}

void SystemMatrixDensityRatioEstimation::mult(base::DataVector& alpha, base::DataVector& result) {
  base::DataVector tempQ(this->paddedInstancesP);

  // Compute (B_q^T * alpha)
  this->myTimer_->start();
  this->B_q->mult(alpha, tempQ);
  this->completeTimeMult_ += this->myTimer_->stop();
  this->computeTimeMult_ += this->B_q->getDuration();

  // Compute (B_q * B_q^T * alpha)
  this->myTimer_->start();
  this->B_q->multTranspose(tempQ, result);
  this->completeTimeMultTrans_ += this->myTimer_->stop();
  this->computeTimeMultTrans_ += this->B_q->getDuration();

  // Compute (B_q * B_q^T + lambda * nq * I) * alpha
  result.axpy(static_cast<double>(this->instancesQ) * this->lambda_, alpha);

  // Compute (np * B_q * B_q^T + lambda * np * nq * I) * alpha
  result.mult(static_cast<double>(this->instancesP));
}

void SystemMatrixDensityRatioEstimation::generateb(base::DataVector& b) {
  base::DataVector y(this->instancesP, static_cast<double>(this->instancesQ));

  // Compute nq * B_p * 1
  this->myTimer_->start();
  this->B_p->multTranspose(y, b);
  this->completeTimeMultTrans_ += this->myTimer_->stop();
  this->computeTimeMultTrans_ += this->B_p->getDuration();
}

void SystemMatrixDensityRatioEstimation::prepareGrid() {
  this->B_p->prepare();
  this->B_q->prepare();
}

void SystemMatrixDensityRatioEstimation::setImplementation(
    datadriven::OperationMultipleEvalConfiguration operationConfiguration) {
  this->implementationConfiguration = operationConfiguration;
  this->B_p.reset(op_factory::createOperationMultipleEval(this->grid, this->datasetP_,
                                                          this->implementationConfiguration));
  this->B_q.reset(op_factory::createOperationMultipleEval(this->grid, this->datasetQ_,
                                                          this->implementationConfiguration));
}

}  // namespace datadriven
}  // namespace sgpp
