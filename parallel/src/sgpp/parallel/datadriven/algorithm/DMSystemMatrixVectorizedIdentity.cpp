// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/exception/operation_exception.hpp>

#include <sgpp/parallel/datadriven/tools/DMVectorizationPaddingAssistant.hpp>
#include <sgpp/parallel/datadriven/algorithm/DMSystemMatrixVectorizedIdentity.hpp>
#include <sgpp/parallel/operation/ParallelOpFactory.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace parallel {

DMSystemMatrixVectorizedIdentity::DMSystemMatrixVectorizedIdentity(
    sgpp::base::Grid& SparseGrid, sgpp::base::DataMatrix& trainData, double lambda,
    VectorizationType vecMode)
    : DMSystemMatrixBase(trainData, lambda),
      vecMode_(vecMode),
      numTrainingInstances_(0),
      numPatchedTrainingInstances_(0) {
  // handle unsupported vector extensions
  if (this->vecMode_ != X86SIMD && this->vecMode_ != MIC && this->vecMode_ != Hybrid_X86SIMD_MIC &&
      this->vecMode_ != OpenCL &&
      this->vecMode_ != Hybrid_X86SIMD_OpenCL) {
    throw sgpp::base::operation_exception(
        "DMSystemMatrixSPVectorizedIdentity : un-supported vector extension!");
  }

  // already initialized in constructor DMSystemMatrixBase
  //  this->dataset_ = new sgpp::base::DataMatrix(trainData);
  this->numTrainingInstances_ = this->dataset_.getNrows();
  this->numPatchedTrainingInstances_ =
      sgpp::parallel::DMVectorizationPaddingAssistant::padDataset(this->dataset_, vecMode_);

  this->dataset_.transpose();

  this->B_.reset(sgpp::op_factory::createOperationMultipleEvalVectorized(SparseGrid, this->vecMode_,
                                                                         &this->dataset_));
}

DMSystemMatrixVectorizedIdentity::~DMSystemMatrixVectorizedIdentity() {
  //    delete this->dataset_;
}

void DMSystemMatrixVectorizedIdentity::mult(sgpp::base::DataVector& alpha,
                                            sgpp::base::DataVector& result) {
  sgpp::base::DataVector temp(this->numPatchedTrainingInstances_);

  // Operation B
  this->myTimer_->start();
  this->computeTimeMult_ += this->B_->multVectorized(alpha, temp);
  this->completeTimeMult_ += this->myTimer_->stop();

  // patch result -> set additional entries zero
  if (this->numTrainingInstances_ != temp.getSize()) {
    for (size_t i = 0; i < (temp.getSize() - this->numTrainingInstances_); i++) {
      temp.set(temp.getSize() - (i + 1), 0.0f);
    }
  }

  this->myTimer_->start();
  this->computeTimeMultTrans_ += this->B_->multTransposeVectorized(temp, result);
  this->completeTimeMultTrans_ += this->myTimer_->stop();

  result.axpy(static_cast<double>(this->numTrainingInstances_) * this->lambda_, alpha);
}

void DMSystemMatrixVectorizedIdentity::generateb(sgpp::base::DataVector& classes,
                                                 sgpp::base::DataVector& b) {
  sgpp::base::DataVector myClasses(classes);

  // Apply padding
  if (this->numPatchedTrainingInstances_ != myClasses.getSize()) {
    myClasses.resizeZero(this->numPatchedTrainingInstances_);
  }

  this->myTimer_->start();
  this->computeTimeMultTrans_ += this->B_->multTransposeVectorized(myClasses, b);
  this->completeTimeMultTrans_ += this->myTimer_->stop();
}

void DMSystemMatrixVectorizedIdentity::rebuildLevelAndIndex() { this->B_->rebuildLevelAndIndex(); }

}  // namespace parallel
}  // namespace sgpp
