// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/exception/operation_exception.hpp>

#include <sgpp/parallel/datadriven/tools/DMVectorizationPaddingAssistant.hpp>
#include <sgpp/parallel/datadriven/algorithm/DMSystemMatrixSPVectorizedIdentity.hpp>
#include <sgpp/parallel/operation/SPParallelOpFactory.hpp>

#include <sgpp/globaldef.hpp>

// TODO David
#if USE_DOUBLE_PRECISION==0


namespace SGPP {
  namespace parallel {

    DMSystemMatrixSPVectorizedIdentity::DMSystemMatrixSPVectorizedIdentity(SGPP::base::Grid& SparseGrid, SGPP::base::DataMatrixSP& trainData, float lambda, VectorizationType vecMode)
      :  DMSystemMatrixBaseSP(trainData, lambda), vecMode_(vecMode), numTrainingInstances_(0), numPatchedTrainingInstances_(0) {
      // handle unsupported vector extensions
      if (this->vecMode_ != X86SIMD && this->vecMode_ != MIC && this->vecMode_ != Hybrid_X86SIMD_MIC && this->vecMode_ != OpenCL && this->vecMode_ != ArBB && this->vecMode_ != Hybrid_X86SIMD_OpenCL && this->vecMode_ != CUDA) {
        throw SGPP::base::operation_exception("DMSystemMatrixSPVectorizedIdentity : un-supported vector extension!");
      }

      // create the operations needed in ApplyMatrix
      this->dataset_ = new SGPP::base::DataMatrixSP(trainData);
      this->numTrainingInstances_ = this->dataset_->getNrows();
      this->numPatchedTrainingInstances_ = SGPP::parallel::DMVectorizationPaddingAssistant::padDataset(*(this->dataset_), vecMode_);

      if (this->vecMode_ != ArBB) {
        this->dataset_->transpose();
      }

      this->B_ = SGPP::op_factory::createOperationMultipleEvalVectorizedSP(SparseGrid, this->vecMode_, this->dataset_);
    }

    DMSystemMatrixSPVectorizedIdentity::~DMSystemMatrixSPVectorizedIdentity() {
      delete this->B_;
      delete this->dataset_;
    }

    void DMSystemMatrixSPVectorizedIdentity::mult(SGPP::base::DataVectorSP& alpha, SGPP::base::DataVectorSP& result) {
      SGPP::base::DataVectorSP temp(this->numPatchedTrainingInstances_);

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

    void DMSystemMatrixSPVectorizedIdentity::generateb(SGPP::base::DataVectorSP& classes, SGPP::base::DataVectorSP& b) {
      SGPP::base::DataVectorSP myClasses(classes);

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

#endif
