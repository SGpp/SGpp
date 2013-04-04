/* ****************************************************************************
* Copyright (C) 2010 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)
// @author Roman Karlstetter (karlstetter@mytum.de)

#ifndef OPERATIONMULTIPLEEVALITERATIVESP_H
#define OPERATIONMULTIPLEEVALITERATIVESP_H

#include "parallel/datadriven/operation/OperationMultipleEvalVectorizedSP.hpp"
#include "parallel/tools/PartitioningTool.hpp"

namespace sg {
  namespace parallel {

    template<typename KernelImplementation>
    class OperationMultipleEvalIterativeSP : public OperationMultipleEvalVectorizedSP {
      public:
        /**
         * Constructor of OperationMultipleEvalIterativeSPX86Simd
         *
         * Within the constructor sg::base::DataMatrixSP Level and sg::base::DataMatrixSP Index are set up.
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
        OperationMultipleEvalIterativeSP(base::GridStorage* storage, base::DataMatrixSP* dataset,
                                         int gridFrom, int gridTo, int datasetFrom, int datasetTo):
          OperationMultipleEvalVectorizedSP(storage, dataset) {
          m_gridFrom = gridFrom;
          m_gridTo = gridTo;
          m_datasetFrom = datasetFrom;
          m_datasetTo = datasetTo;

          rebuildLevelAndIndex();
        }

        virtual double multVectorized(sg::base::DataVectorSP& alpha, sg::base::DataVectorSP& result) {
          myTimer_->start();
          result.setAll(0.0f);

          #pragma omp parallel
          {
            size_t start;
            size_t end;
            PartitioningTool::getOpenMPPartitionSegment(m_datasetFrom, m_datasetTo, &start, &end, KernelImplementation::getChunkDataPoints());

            KernelImplementation::mult(
              level_,
              index_,
              mask_,
              offset_,
              dataset_,
              alpha,
              result,
              0,
              alpha.getSize(),
              start,
              end);
          }
          return myTimer_->stop();
        }

        virtual double multTransposeVectorized(sg::base::DataVectorSP& source, sg::base::DataVectorSP& result) {
          myTimer_->start();
          result.setAll(0.0);

          #pragma omp parallel
          {
            size_t start;
            size_t end;
            PartitioningTool::getOpenMPPartitionSegment(m_gridFrom, m_gridTo, &start, &end, 1);

            KernelImplementation::multTranspose(
              level_,
              index_,
              mask_,
              offset_,
              dataset_,
              source,
              result,
              start,
              end,
              0,
              dataset_->getNcols());
          }

          return myTimer_->stop();
        }

        virtual void rebuildLevelAndIndex() {
          LevelIndexMaskOffsetHelperSP::rebuild<KernelImplementation::kernelType, OperationMultipleEvalVectorizedSP >(this);
        }

        virtual void updateGridComputeBoundaries(int gridFrom, int gridTo) {
          m_gridFrom = gridFrom;
          m_gridTo = gridTo;
        }
    };

  }
}
#endif // OPERATIONMULTIPLEEVALITERATIVESP_H
