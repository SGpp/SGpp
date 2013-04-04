/******************************************************************************
* Copyright (C) 2013 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "parallel/datadriven/basis/modlinear/operation/OperationMultipleEvalIterativeSPX86SimdModLinearMask.hpp"
#include "parallel/datadriven/basis/modlinear/operation/impl/SPX86SimdModLinearMask.hpp"
#include "parallel/tools/PartitioningTool.hpp"

namespace sg {
  namespace parallel {

    OperationMultipleEvalIterativeSPX86SimdModLinearMask::OperationMultipleEvalIterativeSPX86SimdModLinearMask(
      sg::base::GridStorage* storage, sg::base::DataMatrixSP* dataset) :
      sg::parallel::OperationMultipleEvalVectorizedSP(storage, dataset) {
      rebuildLevelAndIndex();
    }

    void OperationMultipleEvalIterativeSPX86SimdModLinearMask::rebuildLevelAndIndex() {
      LevelIndexMaskOffsetHelperSP::rebuild<SPX86SimdModLinearMask::kernelType, OperationMultipleEvalVectorizedSP>(this);
    }

    double OperationMultipleEvalIterativeSPX86SimdModLinearMask::multTransposeVectorized(sg::base::DataVectorSP& source, sg::base::DataVectorSP& result) {
      myTimer_->start();
      result.setAll(0.0);

      #pragma omp parallel
      {
        size_t start;
        size_t end;
        PartitioningTool::getOpenMPPartitionSegment(m_gridFrom, m_gridTo, &start, &end, 1);

        SPX86SimdModLinearMask::multTranspose(
          level_, index_, mask_, offset_, dataset_, source, result, start, end, 0, this->dataset_->getNcols());
      }

      return myTimer_->stop();
    }

    double OperationMultipleEvalIterativeSPX86SimdModLinearMask::multVectorized(sg::base::DataVectorSP& alpha, sg::base::DataVectorSP& result) {
      myTimer_->start();
      result.setAll(0.0);

      #pragma omp parallel
      {
        size_t start;
        size_t end;
        PartitioningTool::getOpenMPPartitionSegment(m_datasetFrom, m_datasetTo, &start, &end, SPX86SimdModLinearMask::getChunkDataPoints());

        SPX86SimdModLinearMask::mult(
          level_, index_, mask_, offset_, dataset_, alpha, result, 0, alpha.getSize(), start, end);
      }

      return myTimer_->stop();
    }

  }
}
