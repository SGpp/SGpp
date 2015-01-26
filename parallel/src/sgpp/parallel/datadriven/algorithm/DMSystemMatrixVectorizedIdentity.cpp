/******************************************************************************
* Copyright (C) 2010 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include <sgpp/base/exception/operation_exception.hpp>

#include <sgpp/parallel/datadriven/tools/DMVectorizationPaddingAssistant.hpp>
#include <sgpp/parallel/datadriven/algorithm/DMSystemMatrixVectorizedIdentity.hpp>
#include <sgpp/parallel/operation/ParallelOpFactory.hpp>

namespace sg {
  namespace parallel {

    DMSystemMatrixVectorizedIdentity::DMSystemMatrixVectorizedIdentity(sg::base::Grid& SparseGrid, sg::base::DataMatrix& trainData, double lambda, VectorizationType vecMode)
      : DMSystemMatrixBase(trainData, lambda), vecMode_(vecMode), numTrainingInstances_(0), numPatchedTrainingInstances_(0) {
      // handle unsupported vector extensions
      if (this->vecMode_ != X86SIMD && this->vecMode_ != MIC && this->vecMode_ != Hybrid_X86SIMD_MIC && this->vecMode_ != OpenCL && this->vecMode_ != ArBB && this->vecMode_ != Hybrid_X86SIMD_OpenCL) {
        throw sg::base::operation_exception("DMSystemMatrixSPVectorizedIdentity : un-supported vector extension!");
      }

      this->dataset_ = new sg::base::DataMatrix(trainData);
      this->numTrainingInstances_ = this->dataset_->getNrows();
      this->numPatchedTrainingInstances_ = sg::parallel::DMVectorizationPaddingAssistant::padDataset(*(this->dataset_), vecMode_);

      if (this->vecMode_ != ArBB) {
        this->dataset_->transpose();
      }

      this->B_ = sg::op_factory::createOperationMultipleEvalVectorized(SparseGrid, this->vecMode_, this->dataset_);
    }

    DMSystemMatrixVectorizedIdentity::~DMSystemMatrixVectorizedIdentity() {
      delete this->B_;
      delete this->dataset_;
    }

    void DMSystemMatrixVectorizedIdentity::mult(sg::base::DataVector& alpha, sg::base::DataVector& result) {
      sg::base::DataVector temp(this->numPatchedTrainingInstances_);

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

      //@TODO make MPI version of this

      this->myTimer_->start();
      this->computeTimeMultTrans_ += this->B_->multTransposeVectorized(temp, result);
      this->completeTimeMultTrans_ += this->myTimer_->stop();

      result.axpy(static_cast<double>(this->numTrainingInstances_)*this->lambda_, alpha);
    }

    void DMSystemMatrixVectorizedIdentity::generateb(sg::base::DataVector& classes, sg::base::DataVector& b) {
      sg::base::DataVector myClasses(classes);

      // Apply padding
      if (this->numPatchedTrainingInstances_ != myClasses.getSize()) {
        myClasses.resizeZero(this->numPatchedTrainingInstances_);
      }

      this->myTimer_->start();
      this->computeTimeMultTrans_ += this->B_->multTransposeVectorized(myClasses, b);
      this->completeTimeMultTrans_ += this->myTimer_->stop();
    }

    void DMSystemMatrixVectorizedIdentity::rebuildLevelAndIndex() {
      this->B_->rebuildLevelAndIndex();
    }

  }
}
