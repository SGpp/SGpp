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
    class DMSystemMatrixVectorizedIdentityAllreduce : public sg::parallel::DMSystemMatrixVectorizedIdentityMPIBase<KernelImplementation::kernelType> {
      private:
        // which part of the dataset this process handles
        size_t data_size;
        size_t data_offset;

      public:
        DMSystemMatrixVectorizedIdentityAllreduce(sg::base::Grid& SparseGrid, sg::base::DataMatrix& trainData, double lambda, VectorizationType vecMode)
          : DMSystemMatrixVectorizedIdentityMPIBase<KernelImplementation::kernelType>(SparseGrid, trainData, lambda, vecMode) {
          sg::parallel::PartitioningTool::getMPIPartitionSegment(this->numPatchedTrainingInstances_, &data_size, &data_offset, sg::parallel::DMVectorizationPaddingAssistant::getVecWidth(this->vecMode_));
          rebuildLevelAndIndex();
        }

        virtual void mult(base::DataVector& alpha, base::DataVector& result) {
          this->tempData->setAll(0.0);
          this->result_tmp->setAll(0.0);
          result.setAll(0.0);

          this->myTimer_->start();
          #pragma omp parallel
          {
            size_t threadChunkStart, threadChunkEnd;
            sg::parallel::PartitioningTool::getOpenMPPartitionSegment(
              data_offset, data_offset + data_size,
              &threadChunkStart, &threadChunkEnd,
              sg::parallel::DMVectorizationPaddingAssistant::getVecWidth(this->vecMode_));

            KernelImplementation::mult(
              this->level_,
              this->index_,
              this->mask_,
              this->offset_,
              this->dataset_,
              alpha,
              *(this->tempData),
              0,
              this->storage_->size(),
              threadChunkStart,
              threadChunkEnd);

            // patch result -> set additional entries zero
            // only done for processes that need this part of the temp data for multTrans
            for (size_t i = std::max<size_t>(this->numTrainingInstances_, threadChunkStart); i < threadChunkEnd; i++) {
              this->tempData->set(i, 0.0f);
            }

            #pragma omp single
            {
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

            sg::parallel::PartitioningTool::getOpenMPPartitionSegment(
              0, this->storage_->size(),
              &threadChunkStart, &threadChunkEnd, 1);
            KernelImplementation::multTranspose(
              this->level_,
              this->index_,
              this->mask_,
              this->offset_,
              this->dataset_,
              *(this->tempData),
              *(this->result_tmp),
              threadChunkStart,
              threadChunkEnd,
              data_offset,
              data_offset + data_size);
          }
          this->computeTimeMultTrans_ += this->myTimer_->stop();
          MPI_Allreduce(this->result_tmp->getPointer(), result.getPointer(), result.getSize(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
          this->completeTimeMultTrans_ += this->myTimer_->stop();
        }

        virtual void generateb(base::DataVector& classes, base::DataVector& b) {
          this->tempData->setAll(0.0);
          this->tempData->copyFrom(classes);
          this->result_tmp->setAll(0.0);

          #pragma omp parallel
          {
            size_t threadChunkStart, threadChunkEnd;
            sg::parallel::PartitioningTool::getOpenMPPartitionSegment(
              0, this->storage_->size(),
              &threadChunkStart, &threadChunkEnd, 1);
            KernelImplementation::multTranspose(
              this->level_,
              this->index_,
              this->mask_,
              this->offset_,
              this->dataset_,
              *(this->tempData),
              *(this->result_tmp),
              threadChunkStart,
              threadChunkEnd,
              data_offset,
              data_offset + data_size);
          }
          MPI_Allreduce(this->result_tmp->getPointer(), b.getPointer(), b.getSize(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        }

        virtual void rebuildLevelAndIndex() {
          DMSystemMatrixVectorizedIdentityMPIBase<KernelImplementation::kernelType>::rebuildLevelAndIndex();
        }

    };

  }
}

#endif // DMSYSTEMMATRIXVECTORIZEDIDENTITYALLREDUCE_H
