/* ****************************************************************************
* Copyright (C) 2010 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)
// @author Roman Karlstetter (karlstetter@mytum.de)

#ifndef DMSYSTEMMATRIXVECTORIZEDIDENTITYALLREDUCE_H
#define DMSYSTEMMATRIXVECTORIZEDIDENTITYALLREDUCE_H

#include "base/datatypes/DataVector.hpp"
#include "base/exception/operation_exception.hpp"
#include "base/grid/Grid.hpp"
#include "parallel/datadriven/algorithm/DMSystemMatrixVectorizedIdentityMPIBase.hpp"
#include "parallel/datadriven/tools/DMVectorizationPaddingAssistant.hpp"
#include "parallel/datadriven/tools/LevelIndexMaskOffsetHelper.hpp"
#include "parallel/tools/MPI/SGppMPITools.hpp"
#include "parallel/tools/PartitioningTool.hpp"
#include "parallel/tools/TypesParallel.hpp"

namespace sg {
  namespace parallel {
    template<typename KernelImplementation>
    class DMSystemMatrixVectorizedIdentityAllreduce : public sg::parallel::DMSystemMatrixVectorizedIdentityMPIBase<KernelImplementation> {
      private:
        // which part of the dataset this process handles, it always handles the complete grid
        size_t data_size;
        size_t data_offset;

      public:
        DMSystemMatrixVectorizedIdentityAllreduce(sg::base::Grid& SparseGrid, sg::base::DataMatrix& trainData, double lambda, VectorizationType vecMode)
          : DMSystemMatrixVectorizedIdentityMPIBase<KernelImplementation>(SparseGrid, trainData, lambda, vecMode) {
          sg::parallel::PartitioningTool::getMPIPartitionSegment(this->numPatchedTrainingInstances_, &data_size, &data_offset,                sg::parallel::DMVectorizationPaddingAssistant::getVecWidth(this->vecMode_));
          rebuildLevelAndIndex();
        }

        virtual void mult(base::DataVector& alpha, base::DataVector& result) {
#ifdef X86_MIC_SYMMETRIC
          myGlobalMPIComm->broadcastGridCoefficientsFromRank0(alpha);
#endif
          this->tempData->setAll(0.0);
          this->result_tmp->setAll(0.0);
          result.setAll(0.0);

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
                this->tempData->set(i, 0.0);
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
          myGlobalMPIComm->allreduceSum(*(this->result_tmp), result);
          this->completeTimeMultTrans_ += this->myTimer_->stop();

          result.axpy(static_cast<double>(this->numTrainingInstances_)*this->lambda_, alpha);
        }

        virtual void generateb(base::DataVector& classes, base::DataVector& b) {
          this->tempData->setAll(0.0);
          this->tempData->copyFrom(classes);
          this->result_tmp->setAll(0.0);

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
          myGlobalMPIComm->allreduceSum(*(this->result_tmp), b);
        }

        virtual void rebuildLevelAndIndex() {
          DMSystemMatrixVectorizedIdentityMPIBase<KernelImplementation>::rebuildLevelAndIndex();
        }

    };

  }
}

#endif // DMSYSTEMMATRIXVECTORIZEDIDENTITYALLREDUCE_H
