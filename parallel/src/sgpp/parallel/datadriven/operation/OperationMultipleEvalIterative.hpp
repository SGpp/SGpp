// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at 
// sgpp.sparsegrids.org

#ifndef OPERATIONMULTIPLEEVALITERATIVE_H
#define OPERATIONMULTIPLEEVALITERATIVE_H

#include <sgpp/parallel/datadriven/operation/OperationMultipleEvalVectorized.hpp>
#include <sgpp/base/tools/AlignedMemory.hpp>
#include <sgpp/parallel/tools/PartitioningTool.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace parallel {

    template<typename KernelImplementation>
    class OperationMultipleEvalIterative : public OperationMultipleEvalVectorized {
      public:
        /**
         * Constructor of OperationMultipleEvalIterativeSPX86Simd
         *
         * Within the constructor SGPP::base::DataMatrixSP Level and SGPP::base::DataMatrixSP Index are set up.
         * If the grid changes during your calculations and you don't want to create
         * a new instance of this class, you have to call rebuildLevelAndIndex before
         * doing any further mult or multTranspose calls.
         *
         * @param storage Pointer to the grid's gridstorage object
         * @param dataset dataset that should be evaluated
         * @param gridFrom local part of grid (start)
         * @param gridTo local part of grid (end)
         * @param datasetFrom local part of dataset (start)
         * @param datasetTo local part of dataset (end)
         */
        OperationMultipleEvalIterative(base::GridStorage* storage, base::DataMatrix* dataset,
                                       size_t gridFrom, size_t gridTo, size_t datasetFrom, size_t datasetTo):
          OperationMultipleEvalVectorized(storage, dataset) {
          m_gridFrom = gridFrom;
          m_gridTo = gridTo;
          m_datasetFrom = datasetFrom;
          m_datasetTo = datasetTo;

          rebuildLevelAndIndex(m_gridFrom, m_gridTo);
        }

        virtual double multVectorized(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result) {
          myTimer_->start();
          result.setAll(0.0);

          #pragma omp parallel
          {
            m_kernel.mult(
              level_,
              index_,
              mask_,
              offset_,
              dataset_,
              alpha,
              result,
              0,
              alpha.getSize(),
              m_datasetFrom,
              m_datasetTo);
          }

          return myTimer_->stop();
        }

        virtual double multTransposeVectorized(SGPP::base::DataVector& source, SGPP::base::DataVector& result) {
          myTimer_->start();
          result.setAll(0.0);

          #pragma omp parallel
          {
            m_kernel.multTranspose(
              level_,
              index_,
              mask_,
              offset_,
              dataset_,
              source,
              result,
              m_gridFrom,
              m_gridTo,
              0,
              dataset_->getNcols());
          }

          return myTimer_->stop();
        }

        virtual void rebuildLevelAndIndex(size_t gridFrom = 0, size_t gridTo = std::numeric_limits<size_t>::max()) {
          LevelIndexMaskOffsetHelper::rebuild<KernelImplementation::kernelType, OperationMultipleEvalVectorized >(this);

          if (gridTo == std::numeric_limits<size_t>::max()) {
            gridTo = this->storage_->size();
          }

          m_gridFrom = gridFrom;
          m_gridTo = gridTo;
          m_kernel.resetKernel();
        }
      private:
        KernelImplementation m_kernel;
    };

  }
}
#endif // OPERATIONMULTIPLEEVALITERATIVE_H