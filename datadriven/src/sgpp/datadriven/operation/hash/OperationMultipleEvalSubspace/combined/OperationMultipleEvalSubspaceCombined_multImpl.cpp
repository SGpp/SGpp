// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/datadriven/operation/hash/OperationMultipleEvalSubspace/combined/OperationMultipleEvalSubspaceCombined.hpp>
#include <algorithm>
#include <iomanip>
#include <limits>
#include <utility>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace datadriven {

/**
 * Internal mult operator, should not be called directly.
 *
 * @see OperationMultipleEval
 *
 * @param source source operand for the operator
 * @param result stores the result
 * @param start_index_data beginning of the range to process
 * @param end_index_data end of the range to process
 */
void OperationMultipleEvalSubspaceCombined::multImpl(sgpp::base::DataVector& source,
                                                     sgpp::base::DataVector& result,
                                                     const size_t start_index_data,
                                                     const size_t end_index_data) {
  size_t tid = omp_get_thread_num();

  if (tid == 0) {
    this->setCoefficients(source);
  }

#pragma omp barrier

  size_t dim = this->paddedDataset->getNcols();
  double* datasetPtr = this->paddedDataset->getPointer();

  size_t totalThreadNumber = X86COMBINED_PARALLEL_DATA_POINTS + X86COMBINED_VEC_PADDING;

  double* evalIndexValuesAll = new double[(dim + 1) * totalThreadNumber];

  for (size_t i = 0; i < (dim + 1) * totalThreadNumber; i++) {
    evalIndexValuesAll[i] = 1.0;
  }

  // for faster index flattening, last element is for padding
  uint32_t* intermediatesAll = new uint32_t[(dim + 1) * totalThreadNumber];

  for (size_t i = 0; i < (dim + 1) * totalThreadNumber; i++) {
    intermediatesAll[i] = 0.0;
  }

  size_t validIndices[X86COMBINED_PARALLEL_DATA_POINTS + X86COMBINED_VEC_PADDING];
  size_t validIndicesCount;

  double componentResults[X86COMBINED_PARALLEL_DATA_POINTS + X86COMBINED_VEC_PADDING];
  size_t levelIndices[X86COMBINED_PARALLEL_DATA_POINTS + X86COMBINED_VEC_PADDING];
  // size_t nextIterationToRecalcReferences[X86COMBINED_PARALLEL_DATA_POINTS +
  // X86COMBINED_VEC_PADDING];

  double* listSubspace = new double[this->maxGridPointsOnLevel];

  for (size_t i = 0; i < this->maxGridPointsOnLevel; i++) {
    listSubspace[i] = std::numeric_limits<double>::quiet_NaN();
  }

  // process the next chunk of data tuples in parallel
  for (size_t dataIndexBase = start_index_data; dataIndexBase < end_index_data;
       dataIndexBase += X86COMBINED_PARALLEL_DATA_POINTS) {
    for (size_t i = 0; i < totalThreadNumber; i++) {
      levelIndices[i] = 0.0;
      componentResults[i] = 0.0;
      // nextIterationToRecalcReferences[i] = 0;
    }

    for (size_t subspaceIndex = 0; subspaceIndex < subspaceCount; subspaceIndex++) {
      SubspaceNodeCombined& subspace = this->allSubspaceNodes[subspaceIndex];

      double* levelArrayContinuous = nullptr;

      // prepare the subspace array for a list type subspace
      if (subspace.type == SubspaceNodeCombined::SubspaceType::LIST) {
        // fill with surplusses
        for (std::pair<uint32_t, double> tuple : subspace.indexFlatSurplusPairs) {
          // actual values are utilized, but only read
          listSubspace[tuple.first] = tuple.second;
        }

        levelArrayContinuous = listSubspace;
      } else {
        levelArrayContinuous = subspace.subspaceArray.data();
      }

      validIndicesCount = 0;

      for (size_t parallelIndex = 0; parallelIndex < X86COMBINED_PARALLEL_DATA_POINTS;
           parallelIndex++) {
        size_t parallelLevelIndex = levelIndices[parallelIndex];

        if (parallelLevelIndex == subspaceIndex) {
          validIndices[validIndicesCount] = parallelIndex;
          validIndicesCount += 1;
        }
      }

      size_t paddingSize = std::min(static_cast<int>(validIndicesCount + X86COMBINED_VEC_PADDING),
                                    X86COMBINED_PARALLEL_DATA_POINTS + X86COMBINED_VEC_PADDING);

      for (size_t i = validIndicesCount; i < paddingSize; i++) {
        size_t threadId = X86COMBINED_PARALLEL_DATA_POINTS + (i - validIndicesCount);
        validIndices[i] = threadId;
        componentResults[threadId] = 0.0;
        levelIndices[threadId] = 0;
        // nextIterationToRecalcReferences[threadId] = 0;
        double* evalIndexValues = evalIndexValuesAll + (dim + 1) * threadId;

        // for faster index flattening, last element is for padding
        uint32_t* intermediates = intermediatesAll + (dim + 1) * threadId;

        for (size_t j = 0; j < dim; j++) {
          evalIndexValues[j] = 1.0;
          intermediates[j] = 0.0;
        }
      }

      uncachedMultTransposeInner(dim, datasetPtr, dataIndexBase, end_index_data, subspace,
                                 levelArrayContinuous, validIndicesCount, validIndices,
                                 levelIndices,  // nextIterationToRecalcReferences,
                                 componentResults, evalIndexValuesAll, intermediatesAll);

      if (subspace.type == SubspaceNodeCombined::SubspaceType::LIST) {
        for (std::pair<uint32_t, double>& tuple : subspace.indexFlatSurplusPairs) {
          listSubspace[tuple.first] = std::numeric_limits<double>::quiet_NaN();
        }
      }
    }  // end iterate grid

    for (size_t parallelIndex = 0; parallelIndex < X86COMBINED_PARALLEL_DATA_POINTS;
         parallelIndex++) {
      size_t dataIndex = dataIndexBase + parallelIndex;
      result.set(dataIndex, componentResults[parallelIndex]);
    }
  }  // end iterate data chunks

  delete[] evalIndexValuesAll;
  delete[] intermediatesAll;
  delete[] listSubspace;
}
}  // namespace datadriven
}  // namespace sgpp
