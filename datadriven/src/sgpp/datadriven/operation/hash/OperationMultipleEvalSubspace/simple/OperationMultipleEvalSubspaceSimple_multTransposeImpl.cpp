// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/datadriven/operation/hash/OperationMultipleEvalSubspace/simple/OperationMultipleEvalSubspaceSimple.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace datadriven {

    void OperationMultipleEvalSubspaceSimple::multTransposeImpl(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result,
        //const size_t start_index_grid,
        //const size_t end_index_grid,
        const size_t start_index_data, const size_t end_index_data) {

      size_t tid = omp_get_thread_num();

      if (tid == 0) {
        this->setCoefficients(result);
      }

      #pragma omp barrier

      size_t dim = this->dataset.getNcols();

      base::DataVector dataTuple(dim);
      float_t* dataTuplePtr = dataTuple.getPointer();

      size_t* indexPtr = new size_t[dim];
      float_t* evalIndexValues = new float_t[dim + 1];

      //float_t **flatLevels = this->flatLevels;

      //for faster index flattening
      size_t* intermediates = new size_t[dim + 1];

      //float_t maxIndex = sizeof(allSubspaces) / (subspaceSize * sizeof(size_t));
      float_t maxIndex = static_cast<float_t>(subspaceCount * subspaceSize);

      for (size_t dataIndex = start_index_data; dataIndex < end_index_data; dataIndex++) {

        evalIndexValues[0] = 1.0;
        intermediates[0] = 0;
        dataset.getRow(dataIndex, dataTuple);
        size_t levelIndex = 0;
        size_t nextIterationToRecalc = 0; //all index components are recalculated

        while (levelIndex < maxIndex) {

          size_t* hInversePtr = allSubspaces + levelIndex + dim;
          //size_t *levelPtr = allSubspaces + levelIndex;
          size_t linearSurplusIndex = *(allSubspaces + levelIndex + (2 * dim) + 1);
          //float_t *levelArray = flatLevels[levelFlat];
          float_t* levelArray = &(this->allSurplusses[linearSurplusIndex]);

          // calculate index
#if X86SIMPLE_ENABLE_PARTIAL_RESULT_REUSAGE == 1

          for (size_t i = nextIterationToRecalc; i < dim; i++) {
#else

          for (size_t i = 0; i < dim; i++) {
#endif
            float_t unadjusted = dataTuplePtr[i] * static_cast<float_t>(hInversePtr[i]);
            indexPtr[i] = calculateIndexComponent(unadjusted);
          }

          size_t indexFlat = this->flattenIndex(intermediates, dim, hInversePtr, indexPtr, nextIterationToRecalc);
          float_t surplus = levelArray[indexFlat];

          if (!std::isnan(surplus)) {

            //prepare the values for the individual components
#if X86SIMPLE_ENABLE_PARTIAL_RESULT_REUSAGE == 1
            float_t phiEval = evalIndexValues[nextIterationToRecalc];

            for (size_t i = nextIterationToRecalc; i < dim; i++) {
#else
            float_t phiEval = 1.0;

            for (size_t i = 0; i < dim; i++) {
#endif
              float_t phi1DEval = static_cast<float_t>(hInversePtr[i]) * dataTuplePtr[i]
                                  - static_cast<float_t>(indexPtr[i]);
              phi1DEval = std::max(0.0, 1.0 - fabs(phi1DEval));
              phiEval *= phi1DEval;
              evalIndexValues[i + 1] = phiEval;
            }

            float_t partialSurplus = phiEval * alpha[dataIndex];

            #pragma omp atomic
            levelArray[indexFlat] += partialSurplus;

            nextIterationToRecalc = allSubspaces[levelIndex + (2 * dim) + 2];
            levelIndex += subspaceSize;
          }
          else {
#if X86SIMPLE_ENABLE_SUBSPACE_SKIPPING == 1
            //skip to next relevant subspace
            nextIterationToRecalc = allSubspaces[levelIndex + (2 * dim) + 3];
            levelIndex = allSubspaces[levelIndex + 2 * dim];
#else
            nextIterationToRecalc = allSubspaces[levelIndex + (2 * dim) + 2];
            levelIndex += subspaceSize;
#endif
          }

        } // end iterate subspaces

      } // end iterate data points

      delete[] indexPtr;
      delete[] evalIndexValues;
      delete[] intermediates;

      #pragma omp barrier

      if (tid == 0) {
        this->unflatten(result);
      }

    }

  }
}
