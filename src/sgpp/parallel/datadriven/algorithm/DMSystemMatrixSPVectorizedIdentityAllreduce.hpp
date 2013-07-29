/* ****************************************************************************
* Copyright (C) 2010 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)
// @author Roman Karlstetter (karlstetter@mytum.de)

#ifndef DMSYSTEMMATRIXSPVECTORIZEDIDENTITYALLREDUCE_H
#define DMSYSTEMMATRIXSPVECTORIZEDIDENTITYALLREDUCE_H

#include "base/exception/operation_exception.hpp"
#include "base/grid/Grid.hpp"
#include "parallel/datadriven/algorithm/DMSystemMatrixSPVectorizedIdentityMPIBase.hpp"
#include "parallel/datadriven/tools/DMVectorizationPaddingAssistant.hpp"
#include "parallel/datadriven/tools/LevelIndexMaskOffsetHelper.hpp"
#include "parallel/tools/MPI/SGppMPITools.hpp"
#include "parallel/tools/PartitioningTool.hpp"
#include "parallel/tools/TypesParallel.hpp"

namespace sg {
  namespace parallel {
    template<typename KernelImplementation>
    class DMSystemMatrixSPVectorizedIdentityAllreduce : public sg::parallel::DMSystemMatrixSPVectorizedIdentityMPIBase<KernelImplementation> {
      private:
        // which part of the dataset this process handles, it always handles the complete grid
        size_t data_size;
        size_t data_offset;

      public:
        DMSystemMatrixSPVectorizedIdentityAllreduce(sg::base::Grid& SparseGrid, sg::base::DataMatrixSP& trainData, float lambda, VectorizationType vecMode)
          : DMSystemMatrixSPVectorizedIdentityMPIBase<KernelImplementation>(SparseGrid, trainData, lambda, vecMode) {
          sg::parallel::PartitioningTool::getMPIPartitionSegment(this->numPatchedTrainingInstances_, &data_size, &data_offset, sg::parallel::DMVectorizationPaddingAssistant::getVecWidthSP(this->vecMode_));
          rebuildLevelAndIndex();
        }

        virtual void mult(base::DataVectorSP& alpha, base::DataVectorSP& result) {
#ifdef X86_MIC_SYMMETRIC
          myGlobalMPIComm->broadcastSPGridCoefficientsFromRank0(alpha);
#endif
          this->tempData->setAll(0.0f);
          this->result_tmp->setAll(0.0f);
          result.setAll(0.0f);

          this->myTimer_->start();
          #pragma omp parallel
          {
            this->kernel_.mult(
              this->level_,
              this->index_,
              this->mask_,
              this->offset_,
              this->dataset_,
              alpha,
              *(this->tempData),
              0,
              this->storage_->size(),
              data_offset,
              data_offset + data_size);

            #pragma omp barrier
            // make sure that all threads finished their part, so that we can
            // safely overwrite the padded range
            #pragma omp single
            {
              // patch result -> set additional entries zero
              for (size_t i = this->numTrainingInstances_; i < this->numPatchedTrainingInstances_; i++) {
                this->tempData->set(i, 0.0f);
              }

              // the time measured here does not represent the complete
              // time spent computing mult, it's probably the first threads
              // that finished mult() that enters here while the other
              // threads might still be busy with mult
              double timeMult = this->myTimer_->stop();
              this->computeTimeMult_ += timeMult;
              this->completeTimeMult_ += timeMult;
              this->myTimer_->start();
            }
            // #pragma omp barrier
            // implicit openMP barrier here (after omp single), which is
            // needed, as multTranspose works on the full data part of this
            // process, so threads might work on unfinished results of mult

            this->kernel_.multTranspose(
              this->level_,
              this->index_,
              this->mask_,
              this->offset_,
              this->dataset_,
              *(this->tempData),
              *(this->result_tmp),
              0,
              this->storage_->size(),
              data_offset,
              data_offset + data_size);
          }
          myGlobalMPIComm->Barrier();
          this->computeTimeMultTrans_ += this->myTimer_->stop();
          myGlobalMPIComm->allreduceSumSP(*(this->result_tmp), result);
          this->completeTimeMultTrans_ += this->myTimer_->stop();

          result.axpy(static_cast<float>(this->numTrainingInstances_)*this->lambda_, alpha);
        }

        virtual void generateb(base::DataVectorSP& classes, base::DataVectorSP& b) {
          this->tempData->setAll(0.0f);
          this->tempData->copyFrom(classes);
          this->result_tmp->setAll(0.0f);

          #pragma omp parallel
          {
            this->kernel_.multTranspose(
              this->level_,
              this->index_,
              this->mask_,
              this->offset_,
              this->dataset_,
              *(this->tempData),
              *(this->result_tmp),
              0,
              this->storage_->size(),
              data_offset,
              data_offset + data_size);
          }
          myGlobalMPIComm->allreduceSumSP(*(this->result_tmp), b);
        }

        virtual void rebuildLevelAndIndex() {
          DMSystemMatrixSPVectorizedIdentityMPIBase<KernelImplementation>::rebuildLevelAndIndex();
        }

    };

  }
}

#endif // DMSYSTEMMATRIXSPVECTORIZEDIDENTITYALLREDUCE_H
