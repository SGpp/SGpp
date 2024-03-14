// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/datadriven/algorithm/SystemMatrixRelativeDensityRatioEstimation.hpp>

#include <sgpp/base/exception/operation_exception.hpp>
#include <sgpp/datadriven/DatadrivenOpFactory.hpp>
#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace datadriven {

SystemMatrixRelativeDensityRatioEstimation::SystemMatrixRelativeDensityRatioEstimation(
    base::Grid& grid, base::DataMatrix& trainDataP, base::DataMatrix& trainDataQ, double lambda,
    double omega)
    : DMSystemMatrixTwoDatasets(trainDataP, trainDataQ, lambda),
      instancesP(0),
      paddedInstancesP(0),
      instancesQ(0),
      paddedInstancesQ(0),
      grid(grid),
      omega(omega) {
  this->instancesP = this->datasetP.getNrows();
  this->instancesQ = this->datasetQ.getNrows();
  this->B_p.reset(op_factory::createOperationMultipleEval(grid, this->datasetP,
                                                          this->implementationConfiguration));
  this->B_q.reset(op_factory::createOperationMultipleEval(grid, this->datasetQ,
                                                          this->implementationConfiguration));

  // padded during Operator construction, fetch new size
  this->paddedInstancesP = this->datasetP.getNrows();
  // padded during Operator construction, fetch new size
  this->paddedInstancesQ = this->datasetQ.getNrows();
}

SystemMatrixRelativeDensityRatioEstimation::~SystemMatrixRelativeDensityRatioEstimation() {}

void SystemMatrixRelativeDensityRatioEstimation::mult(base::DataVector& alpha,
                                                      base::DataVector& result) {
  base::DataVector tempQ(this->paddedInstancesQ);
  base::DataVector tempP(this->paddedInstancesP);
  base::DataVector result2(result.getSize());

  // Compute (B_q^T * alpha)
  this->myTimer->start();
  this->B_q->mult(alpha, tempQ);
  this->completeTimeMult += this->myTimer->stop();
  this->computeTimeMult += this->B_q->getDuration();

  // Compute (B_q * B_q^T * alpha)
  this->myTimer->start();
  this->B_q->multTranspose(tempQ, result);
  this->completeTimeMultTrans += this->myTimer->stop();
  this->computeTimeMultTrans += this->B_q->getDuration();

  // Compute ((1-omega) * np * B_q * B_q^T * alpha)
  result.mult((1. - this->omega) * static_cast<double>(this->instancesP));

  // Compute (B_p^T * alpha)
  this->myTimer->start();
  this->B_p->mult(alpha, tempP);
  this->completeTimeMult += this->myTimer->stop();
  this->computeTimeMult += this->B_p->getDuration();

  // Compute (B_p * B_p^T * alpha)
  this->myTimer->start();
  this->B_p->multTranspose(tempP, result2);
  this->completeTimeMultTrans += this->myTimer->stop();
  this->computeTimeMultTrans += this->B_p->getDuration();

  // Compute ((1 - omega) * np * B_q * B_q^T + omega * nq * B_p * B_p^T) * alpha
  result.axpy(this->omega * static_cast<double>(this->instancesQ), result2);

  // ((1 - omega) * np * B_q * B_q^T + omega * nq * B_p * B_p^T + lambda * np * nq * I) * alpha
  result.axpy(
      static_cast<double>(this->instancesP) * static_cast<double>(this->instancesQ) * this->lambda,
      alpha);
}

void SystemMatrixRelativeDensityRatioEstimation::generateb(base::DataVector& b) {
  base::DataVector y(this->paddedInstancesP, static_cast<double>(this->instancesQ));

  // Compute nq * B_p * 1
  this->myTimer->start();
  this->B_p->multTranspose(y, b);
  this->completeTimeMultTrans += this->myTimer->stop();
  this->computeTimeMultTrans += this->B_p->getDuration();
}

void SystemMatrixRelativeDensityRatioEstimation::prepareGrid() {
  this->B_p->prepare();
  this->B_q->prepare();
}

void SystemMatrixRelativeDensityRatioEstimation::setImplementation(
    datadriven::OperationMultipleEvalConfiguration operationConfiguration) {
  this->implementationConfiguration = operationConfiguration;
  this->B_p.reset(op_factory::createOperationMultipleEval(this->grid, this->datasetP,
                                                          this->implementationConfiguration));
  this->B_q.reset(op_factory::createOperationMultipleEval(this->grid, this->datasetQ,
                                                          this->implementationConfiguration));
  // padded during Operator construction, fetch new size
  this->paddedInstancesP = this->datasetP.getNrows();
  // padded during Operator construction, fetch new size
  this->paddedInstancesQ = this->datasetQ.getNrows();
}

}  // namespace datadriven
}  // namespace sgpp
