// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/datadriven/operation/hash/OperationMultipleEvalSubspace/combined/OperationMultipleEvalSubspaceCombined.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace datadriven {

void OperationMultipleEvalSubspaceCombined::uncachedMultTransposeInner(
    size_t dim, const double* const datasetPtr, size_t dataIndexBase, size_t end_index_data,
    SubspaceNodeCombined& subspace, double* levelArrayContinuous, size_t validIndicesCount,
    size_t* validIndices,
    size_t* levelIndices,  // size_t *nextIterationToRecalcReferences,
    double* componentResults, double* evalIndexValuesAll, uint32_t* intermediatesAll) {
  for (size_t validIndex = 0; validIndex < validIndicesCount;
       validIndex += X86COMBINED_VEC_PADDING) {
    // for (size_t validIndex = 0; validIndex < validIndicesCount; validIndex += 4) {
    size_t parallelIndices[4];
    parallelIndices[0] = validIndices[validIndex];
    parallelIndices[1] = validIndices[validIndex + 1];
    parallelIndices[2] = validIndices[validIndex + 2];
    parallelIndices[3] = validIndices[validIndex + 3];

// cout << "firstRound: " << firstRound << endl;

//        size_t nextIterationToRecalc = nextIterationToRecalcReferences[parallelIndices[0]];

#if X86COMBINED_ENABLE_PARTIAL_RESULT_REUSAGE == 1
    size_t nextIterationToRecalc = subspace.arriveDiff;
#else
    size_t nextIterationToRecalc = 0;
#endif
    //        if (nextIterationToRecalc != subspace.arriveDiff) {
    //            cout << "your are a bit wrong today, arrive:" << subspace.arriveDiff << " ref: "
    //            << nextIterationToRecalc << endl;
    //
    //        }

    const double* const dataTuplePtr[4] = {datasetPtr + (dataIndexBase + parallelIndices[0]) * dim,
                                           datasetPtr + (dataIndexBase + parallelIndices[1]) * dim,
                                           datasetPtr + (dataIndexBase + parallelIndices[2]) * dim,
                                           datasetPtr + (dataIndexBase + parallelIndices[3]) * dim};

    double* evalIndexValues[4];
    evalIndexValues[0] = evalIndexValuesAll + (dim + 1) * parallelIndices[0];
    evalIndexValues[1] = evalIndexValuesAll + (dim + 1) * parallelIndices[1];
    evalIndexValues[2] = evalIndexValuesAll + (dim + 1) * parallelIndices[2];
    evalIndexValues[3] = evalIndexValuesAll + (dim + 1) * parallelIndices[3];

    // for faster index flattening, last element is for padding
    uint32_t* intermediates[4];
    intermediates[0] = intermediatesAll + (dim + 1) * parallelIndices[0];
    intermediates[1] = intermediatesAll + (dim + 1) * parallelIndices[1];
    intermediates[2] = intermediatesAll + (dim + 1) * parallelIndices[2];
    intermediates[3] = intermediatesAll + (dim + 1) * parallelIndices[3];

    uint32_t indexFlat[4];
    double phiEval[4];

#if X86COMBINED_UNROLL == 1
    size_t parallelIndices2[4];
    parallelIndices2[0] = validIndices[validIndex + 4];
    parallelIndices2[1] = validIndices[validIndex + 5];
    parallelIndices2[2] = validIndices[validIndex + 6];
    parallelIndices2[3] = validIndices[validIndex + 7];

    const double* const dataTuplePtr2[4] = {
        datasetPtr + (dataIndexBase + parallelIndices2[0]) * dim,
        datasetPtr + (dataIndexBase + parallelIndices2[1]) * dim,
        datasetPtr + (dataIndexBase + parallelIndices2[2]) * dim,
        datasetPtr + (dataIndexBase + parallelIndices2[3]) * dim};

    double* evalIndexValues2[4];
    evalIndexValues2[0] = evalIndexValuesAll + (dim + 1) * parallelIndices2[0];
    evalIndexValues2[1] = evalIndexValuesAll + (dim + 1) * parallelIndices2[1];
    evalIndexValues2[2] = evalIndexValuesAll + (dim + 1) * parallelIndices2[2];
    evalIndexValues2[3] = evalIndexValuesAll + (dim + 1) * parallelIndices2[3];

    uint32_t* intermediates2[4];
    intermediates2[0] = intermediatesAll + (dim + 1) * parallelIndices2[0];
    intermediates2[1] = intermediatesAll + (dim + 1) * parallelIndices2[1];
    intermediates2[2] = intermediatesAll + (dim + 1) * parallelIndices2[2];
    intermediates2[3] = intermediatesAll + (dim + 1) * parallelIndices2[3];

    uint32_t indexFlat2[4];
    double phiEval2[4];

    OperationMultipleEvalSubspaceCombined::calculateIndexCombined2(
        dim, nextIterationToRecalc, dataTuplePtr, dataTuplePtr2, subspace.hInverse, intermediates,
        intermediates2, evalIndexValues, evalIndexValues2, indexFlat, indexFlat2, phiEval,
        phiEval2);
#else
    OperationMultipleEvalSubspaceCombined::calculateIndexCombined(
        dim, nextIterationToRecalc, dataTuplePtr, subspace.hInverse, intermediates, evalIndexValues,
        indexFlat, phiEval);
#endif

    double surplus[4];
    surplus[0] = levelArrayContinuous[indexFlat[0]];
    surplus[1] = levelArrayContinuous[indexFlat[1]];
    surplus[2] = levelArrayContinuous[indexFlat[2]];
    surplus[3] = levelArrayContinuous[indexFlat[3]];

#if X86COMBINED_UNROLL == 1
    double surplus2[4];
    surplus2[0] = levelArrayContinuous[indexFlat2[0]];
    surplus2[1] = levelArrayContinuous[indexFlat2[1]];
    surplus2[2] = levelArrayContinuous[indexFlat2[2]];
    surplus2[3] = levelArrayContinuous[indexFlat2[3]];
#endif

    for (size_t innerIndex = 0; innerIndex < 4; innerIndex++) {
      if (!std::isnan(surplus[innerIndex])) {
        // cout << parallelIndices[innerIndex] << " -> phi: " << phiEval[innerIndex] << " surplus: "
        // << surplus[innerIndex] << " = " << phiEval[innerIndex] * surplus[innerIndex] << endl;
        /*if (nextRoundResults[innerIndex] != phiEval[innerIndex] * surplus[innerIndex] &&
        firstRound == false) {
            cout << "differ validindex: " << parallelIndices[innerIndex] << " => ref: " <<
        phiEval[innerIndex] * surplus[innerIndex] << " new: " << nextRoundResults[innerIndex] <<
        endl;

        }*/
        componentResults[parallelIndices[innerIndex]] += phiEval[innerIndex] * surplus[innerIndex];
        // nextIterationToRecalcReferences[parallelIndices[innerIndex]] = subspace.nextDiff;
        levelIndices[parallelIndices[innerIndex]] += 1;
      } else {
#if X86COMBINED_ENABLE_SUBSPACE_SKIPPING == 1
        // skip to next relevant subspace
        // nextIterationToRecalcReferences[parallelIndices[innerIndex]] = subspace.jumpDiff;
        levelIndices[parallelIndices[innerIndex]] = subspace.jumpTargetIndex;
#else
        // nextIterationToRecalcReferences[parallelIndices[innerIndex]] = subspace.nextDiff;
        levelIndices[parallelIndices[innerIndex]] += 1;
#endif
      }
    }

#if X86COMBINED_UNROLL == 1

    for (size_t innerIndex = 0; innerIndex < 4; innerIndex++) {
      if (!std::isnan(surplus2[innerIndex])) {
        // cout << parallelIndices2[innerIndex] << " -> phi: " << phiEval2[innerIndex] << " surplus:
        // " << surplus2[innerIndex] << " = " << phiEval2[innerIndex] * surplus2[innerIndex] <<
        // endl;

        // cout << parallelIndices2[innerIndex] << " ->second phi: " << phiEval2[innerIndex] << "
        // surplus: " << surplus2[innerIndex] << " = " << phiEval2[innerIndex] *
        // surplus2[innerIndex] << endl;
        // nextRoundResults[innerIndex] = phiEval2[innerIndex] * surplus2[innerIndex];

        componentResults[parallelIndices2[innerIndex]] +=
            phiEval2[innerIndex] * surplus2[innerIndex];
        // nextIterationToRecalcReferences[parallelIndices2[innerIndex]] = subspace.nextDiff;
        levelIndices[parallelIndices2[innerIndex]] += 1;
      } else {
#if X86COMBINED_ENABLE_SUBSPACE_SKIPPING == 1
        // skip to next relevant subspace
        // nextIterationToRecalcReferences[parallelIndices2[innerIndex]] = subspace.jumpDiff;
        levelIndices[parallelIndices2[innerIndex]] = subspace.jumpTargetIndex;
#else
        // nextIterationToRecalcReferences[parallelIndices2[innerIndex]] = subspace.nextDiff;
        levelIndices[parallelIndices2[innerIndex]] += 1;
#endif
      }
    }

#endif
  }  // end X86COMBINED_PARALLEL_DATA_POINTS
}
}
}
