/******************************************************************************
* Copyright (C) 2010 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "base/exception/operation_exception.hpp"

#include "parallel/datadriven/tools/DMVectorizationPaddingAssistant.hpp"
#include "parallel/datadriven/algorithm/DMSystemMatrixSPVectorizedIdentity.hpp"
#include "parallel/operation/SPParallelOpFactory.hpp"

namespace sg {
  namespace parallel {

    DMSystemMatrixSPVectorizedIdentity::DMSystemMatrixSPVectorizedIdentity(sg::base::Grid& SparseGrid, sg::base::DataMatrixSP& trainData, float lambda, VectorizationType vecMode)
      :  DMSystemMatrixBaseSP(trainData, lambda), vecMode_(vecMode), numTrainingInstances_(0), numPatchedTrainingInstances_(0) {
      // handle unsupported vector extensions
      if (this->vecMode_ != X86SIMD && this->vecMode_ != MIC && this->vecMode_ != Hybrid_X86SIMD_MIC && this->vecMode_ != OpenCL && this->vecMode_ != ArBB && this->vecMode_ != Hybrid_X86SIMD_OpenCL && this->vecMode_ != CUDA) {
        throw sg::base::operation_exception("DMSystemMatrixSPVectorizedIdentity : un-supported vector extension!");
      }

      // create the operations needed in ApplyMatrix
      this->dataset_ = new sg::base::DataMatrixSP(trainData);
      this->numTrainingInstances_ = this->dataset_->getNrows();
      this->numPatchedTrainingInstances_ = sg::parallel::DMVectorizationPaddingAssistant::padDataset(*(this->dataset_), vecMode_);

      if (this->vecMode_ != ArBB && this->vecMode_ != CUDA) {
        this->dataset_->transpose();
      }

      this->B_ = sg::op_factory::createOperationMultipleEvalVectorizedSP(SparseGrid, this->vecMode_, this->dataset_);
    }

    DMSystemMatrixSPVectorizedIdentity::~DMSystemMatrixSPVectorizedIdentity() {
      delete this->B_;
      delete this->dataset_;
    }

    void DMSystemMatrixSPVectorizedIdentity::mult(sg::base::DataVectorSP& alpha, sg::base::DataVectorSP& result) {
      sg::base::DataVectorSP temp(this->numPatchedTrainingInstances_);

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

      result.axpy(static_cast<float>(this->numTrainingInstances_)*this->lambda_, alpha);
    }

    void DMSystemMatrixSPVectorizedIdentity::generateb(sg::base::DataVectorSP& classes, sg::base::DataVectorSP& b) {
      sg::base::DataVectorSP myClasses(classes);

      // Apply padding
      if (this->numPatchedTrainingInstances_ != myClasses.getSize()) {
        myClasses.resizeZero(this->numPatchedTrainingInstances_);
      }

      this->myTimer_->start();
      this->computeTimeMultTrans_ += this->B_->multTransposeVectorized(myClasses, b);
      this->completeTimeMultTrans_ += this->myTimer_->stop();
    }

    void DMSystemMatrixSPVectorizedIdentity::rebuildLevelAndIndex() {
      this->B_->rebuildLevelAndIndex();
    }

  }
}
