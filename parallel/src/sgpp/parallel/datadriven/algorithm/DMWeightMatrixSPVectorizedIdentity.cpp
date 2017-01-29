// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/parallel/datadriven/algorithm/DMWeightMatrixSPVectorizedIdentity.hpp>
#include <sgpp/parallel/operation/SPParallelOpFactory.hpp>
#include <sgpp/base/exception/operation_exception.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace parallel {

DMWeightMatrixSPVectorizedIdentity::DMWeightMatrixSPVectorizedIdentity(
    sgpp::base::Grid& SparseGrid, sgpp::base::DataMatrixSP& trainData, float lambda,
    sgpp::base::DataVectorSP& w, VectorizationType vecMode) {
  // handle unsupported vector extensions
  // one may find a better way to determine vectorwidth
  if (this->vecMode == parallel::X86SIMD) {
    this->vecWidth = 48;
  } else if (this->vecMode == parallel::OpenCL) {
    this->vecWidth = 128;
  } else if (this->vecMode == parallel::Hybrid_X86SIMD_OpenCL) {
    this->vecWidth = 128;
  } else {
    throw sgpp::base::operation_exception(
        "DMWeightMatrixVectorizedIdentity : Only X86SIMD or OCL or HYBRID_X86SIMD_OCL are "
        "supported vector extensions!");
  }

  resetTimers();

  // create the operations needed in ApplyMatrix
  this->vecMode = vecMode;
  this->lamb = lambda;
  this->data = new sgpp::base::DataMatrixSP(trainData);
  this->weight = new sgpp::base::DataVectorSP(w);

  numTrainingInstances = data->getNrows();

  // Assure that data has a even number of instances -> padding might be needed
  size_t remainder = data->getNrows() % this->vecWidth;
  size_t loopCount = this->vecWidth - remainder;

  if (loopCount != this->vecWidth) {
    sgpp::base::DataVectorSP lastRow(data->getNcols());

    for (size_t i = 0; i < loopCount; i++) {
      // resize the data matrix
      data->getRow(data->getNrows() - 1, lastRow);
      data->resize(data->getNrows() + 1);
      data->setRow(data->getNrows() - 1, lastRow);
      // resize the weight vector
      weight->resize(weight->getSize() + 1);
      weight->set(weight->getSize() - 1, 0.0f);
    }
  }

  numPatchedTrainingInstances = data->getNrows();

  if (this->vecMode != OpenCL && this->vecMode != Hybrid_X86SIMD_OpenCL) {
    data->transpose();
  }

  this->myTimer = new sgpp::base::SGppStopwatch();

  this->B = sgpp::op_factory::createOperationMultipleEvalVectorizedSP(SparseGrid, this->vecMode,
                                                                      this->data);
}

DMWeightMatrixSPVectorizedIdentity::~DMWeightMatrixSPVectorizedIdentity() {
  delete this->B;
  delete this->data;
  delete this->weight;
  delete this->myTimer;
}

void DMWeightMatrixSPVectorizedIdentity::mult(sgpp::base::DataVectorSP& alpha,
                                              sgpp::base::DataVectorSP& result) {
  sgpp::base::DataVectorSP temp(numPatchedTrainingInstances);

  // Operation B
  this->myTimer->start();
  this->computeTimeMult += this->B->multVectorized(alpha, temp);
  this->completeTimeMult += this->myTimer->stop();

  // auto set additional entries zero (upward weight vector set the additional entries to zero)
  temp.componentwise_mult(*weight);

  this->myTimer->start();
  this->computeTimeMultTrans += this->B->multTransposeVectorized(temp, result);
  this->completeTimeMultTrans += this->myTimer->stop();

  result.axpy(this->lamb, alpha);
}

void DMWeightMatrixSPVectorizedIdentity::generateb(sgpp::base::DataVectorSP& classes,
                                                   sgpp::base::DataVectorSP& b) {
  sgpp::base::DataVectorSP myClassesWithWeights(classes);
  size_t loopcount = numPatchedTrainingInstances - numTrainingInstances;

  // Apply padding
  // resize the data class vector and weight
  if (numPatchedTrainingInstances != numTrainingInstances) {
    float lastClass;

    for (size_t i = 0; i < loopcount; i++) {
      lastClass = myClassesWithWeights.get(myClassesWithWeights.getSize() - 1);
      myClassesWithWeights.resize(myClassesWithWeights.getSize() + 1);
      myClassesWithWeights.set(myClassesWithWeights.getSize() - 1, lastClass);
    }
  }

  myClassesWithWeights.componentwise_mult(*weight);

  this->myTimer->start();
  this->computeTimeMultTrans += this->B->multTransposeVectorized(myClassesWithWeights, b);
  this->completeTimeMultTrans += this->myTimer->stop();
}

void DMWeightMatrixSPVectorizedIdentity::rebuildLevelAndIndex() { this->B->rebuildLevelAndIndex(); }

void DMWeightMatrixSPVectorizedIdentity::resetTimers() {
  this->completeTimeMult = 0.0;
  this->computeTimeMult = 0.0;
  this->completeTimeMultTrans = 0.0;
  this->computeTimeMultTrans = 0.0;
}

void DMWeightMatrixSPVectorizedIdentity::getTimers(double& timeMult, double& computeMult,
                                                   double& timeMultTrans,
                                                   double& computeMultTrans) {
  timeMult = this->completeTimeMult;
  computeMult = this->computeTimeMult;
  timeMultTrans = this->completeTimeMultTrans;
  computeMultTrans = this->computeTimeMultTrans;
}

}  // namespace parallel
}  // namespace sgpp
