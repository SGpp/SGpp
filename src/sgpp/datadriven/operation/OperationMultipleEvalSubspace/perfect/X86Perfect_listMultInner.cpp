#include "../../OperationMultipleEvalSubspace/perfect/X86Perfect.hpp"

void X86Perfect::listMultInner(size_t dim, const double * const datasetPtr, sg::base::DataVector &alpha,
        size_t dataIndexBase, size_t end_index_data, X86PerfectSubspaceNode &subspace, double *levelArrayContinuous,
        size_t validIndicesCount, size_t *validIndices, size_t *levelIndices, size_t *nextIterationToRecalcReferences,
        size_t nextIterationToRecalc, double *evalIndexValuesAll, uint32_t *intermediatesAll) {

    for (size_t validIndex = 0; validIndex < validIndicesCount; validIndex += X86PERFECT_VEC_PADDING) {

        size_t parallelIndices[4];
        parallelIndices[0] = validIndices[validIndex];
        parallelIndices[1] = validIndices[validIndex + 1];
        parallelIndices[2] = validIndices[validIndex + 2];
        parallelIndices[3] = validIndices[validIndex + 3];

        const double * const dataTuplePtr[4] = { datasetPtr + (dataIndexBase + parallelIndices[0]) * dim, datasetPtr
                + (dataIndexBase + parallelIndices[1]) * dim, datasetPtr + (dataIndexBase + parallelIndices[2]) * dim,
                datasetPtr + (dataIndexBase + parallelIndices[3]) * dim };

        double *evalIndexValues[4];
        evalIndexValues[0] = evalIndexValuesAll + (dim + 1) * parallelIndices[0];
        evalIndexValues[1] = evalIndexValuesAll + (dim + 1) * parallelIndices[1];
        evalIndexValues[2] = evalIndexValuesAll + (dim + 1) * parallelIndices[2];
        evalIndexValues[3] = evalIndexValuesAll + (dim + 1) * parallelIndices[3];

        //for faster index flattening, last element is for padding
        uint32_t *intermediates[4];
        intermediates[0] = intermediatesAll + (dim + 1) * parallelIndices[0];
        intermediates[1] = intermediatesAll + (dim + 1) * parallelIndices[1];
        intermediates[2] = intermediatesAll + (dim + 1) * parallelIndices[2];
        intermediates[3] = intermediatesAll + (dim + 1) * parallelIndices[3];

        uint32_t indexFlat[4];
        double phiEval[4];

        X86Perfect::calculateIndexPerfect(dim, nextIterationToRecalc, dataTuplePtr, subspace.hInverse, intermediates,
                evalIndexValues, indexFlat, phiEval);

        //surplus is only used to find out whether the value actually exists
        //will not have an update hazard here
        double surplus[4];
        bool isLeaf[4];
        surplus[0] = levelArrayContinuous[indexFlat[0]];
        isLeaf[0] = levelArrayContinuous[indexFlat[0] + 1] == 1.0 ? true : false;
        surplus[1] = levelArrayContinuous[indexFlat[1]];
        isLeaf[1] = levelArrayContinuous[indexFlat[1] + 1] == 1.0 ? true : false;
        surplus[2] = levelArrayContinuous[indexFlat[2]];
        isLeaf[2] = levelArrayContinuous[indexFlat[2] + 1] == 1.0 ? true : false;
        surplus[3] = levelArrayContinuous[indexFlat[3]];
        isLeaf[3] = levelArrayContinuous[indexFlat[3] + 1] == 1.0 ? true : false;

        for (size_t innerIndex = 0; innerIndex < 4; innerIndex++) {
            size_t parallelIndex = parallelIndices[innerIndex];
            if (!std::isnan(surplus[innerIndex])) {
                double partialSurplus = 0.0;
                if (dataIndexBase + parallelIndex < end_index_data && parallelIndex < X86PERFECT_PARALLEL_DATA_POINTS) {
                    partialSurplus = phiEval[innerIndex] * alpha[dataIndexBase + parallelIndex];

                    size_t localIndexFlat = indexFlat[innerIndex];

                    //no atomics required, subspace is locked before processing
                    //#pragma omp atomic
                    levelArrayContinuous[localIndexFlat] += partialSurplus;
                }
            }
            // leaf -> can jump, after jump I can continue at another non-existing point -> can jump
            if (isLeaf[innerIndex] || std::isnan(surplus[innerIndex])) {
                nextIterationToRecalcReferences[parallelIndices[innerIndex]] = subspace.jumpDiff;
                levelIndices[parallelIndices[innerIndex]] = subspace.jumpTargetIndex;
            } else {
                nextIterationToRecalcReferences[parallelIndices[innerIndex]] = subspace.nextDiff;
                levelIndices[parallelIndices[innerIndex]] += 1;
            }
        }
/*
        for (size_t innerIndex = 0; innerIndex < 4; innerIndex++) {
            size_t parallelIndex = parallelIndices[innerIndex];
            if (!std::isnan(surplus[innerIndex])) {
                double partialSurplus = 0.0;
                if (dataIndexBase + parallelIndex < end_index_data && parallelIndex < X86PERFECT_PARALLEL_DATA_POINTS) {
                    partialSurplus = phiEval[innerIndex] * alpha[dataIndexBase + parallelIndex];

                    size_t localIndexFlat = indexFlat[innerIndex];

                    //no atomics required, working on temporary arrays
                    //#pragma omp atomic
                    levelArrayContinuous[localIndexFlat] += partialSurplus;
                }
                nextIterationToRecalcReferences[parallelIndex] = subspace.nextDiff;
                levelIndices[parallelIndex] += 1;
            } else {
#if X86PERFECT_ENABLE_SUBSPACE_SKIPPING == 1
                //skip to next relevant subspace
                nextIterationToRecalcReferences[parallelIndex] = subspace.jumpDiff;
                levelIndices[parallelIndex] = subspace.jumpTargetIndex;
#else
                nextIterationToRecalcReferences[parallelIndex] = subspace.nextDiff;
                levelIndices[parallelIndex] += 1;
#endif
            }
        } // end innerIndex
        */
    } // end parallel
}

/*
 void listMultInnerRead(size_t dim,
 const double * const datasetPtr,
 sg::base::DataVector &alpha,
 size_t dataIndexBase,
 size_t end_index_data,
 //size_t levelIndex,
 X86PerfectSubspaceNode &subspace,
 double *levelArrayContinuous,
 size_t validIndicesCount,
 size_t *validIndices,
 size_t *levelIndices,
 size_t *nextIterationToRecalcReferences,
 size_t nextIterationToRecalc,
 double *evalIndexValuesAll,
 uint32_t *intermediatesAll
 ) {

 for (size_t validIndex = 0; validIndex < validIndicesCount; validIndex += X86PERFECT_VEC_PADDING) {

 size_t parallelIndices[4];
 parallelIndices[0] = validIndices[validIndex];
 parallelIndices[1] = validIndices[validIndex + 1];
 parallelIndices[2] = validIndices[validIndex + 2];
 parallelIndices[3] = validIndices[validIndex + 3];

 const double * const dataTuplePtr[4] =
 {datasetPtr + (dataIndexBase + parallelIndices[0]) * dim,
 datasetPtr + (dataIndexBase + parallelIndices[1]) * dim,
 datasetPtr + (dataIndexBase + parallelIndices[2]) * dim,
 datasetPtr + (dataIndexBase + parallelIndices[3]) * dim};

 double *evalIndexValues[4];
 evalIndexValues[0] = evalIndexValuesAll + (dim + 1) * parallelIndices[0];
 evalIndexValues[1] = evalIndexValuesAll + (dim + 1) * parallelIndices[1];
 evalIndexValues[2] = evalIndexValuesAll + (dim + 1) * parallelIndices[2];
 evalIndexValues[3] = evalIndexValuesAll + (dim + 1) * parallelIndices[3];

 //for faster index flattening, last element is for padding
 uint32_t *intermediates[4];
 intermediates[0] = intermediatesAll + (dim + 1) * parallelIndices[0];
 intermediates[1] = intermediatesAll + (dim + 1) * parallelIndices[1];
 intermediates[2] = intermediatesAll + (dim + 1) * parallelIndices[2];
 intermediates[3] = intermediatesAll + (dim + 1) * parallelIndices[3];

 uint32_t indexFlat[4];
 double phiEval[4];

 X86Perfect::calculateIndexPerfect(dim, nextIterationToRecalc,
 dataTuplePtr, subspace.hInverse,
 intermediates,
 evalIndexValues,
 indexFlat, phiEval
 );

 //surplus is only used to find out whether the value actually exists
 //will not have an update hazard here
 double surplus[4];
 // cout << "test: " << indexFlat[0] << endl;
 surplus[0] = levelArrayContinuous[indexFlat[0]];
 surplus[1] = levelArrayContinuous[indexFlat[1]];
 surplus[2] = levelArrayContinuous[indexFlat[2]];
 surplus[3] = levelArrayContinuous[indexFlat[3]];

 for (size_t innerIndex = 0; innerIndex < 4; innerIndex++) {

 //	cout << "read surplus: " << surplus[innerIndex] << endl;
 //TODO: rework ifs (at least remove inner)
 size_t parallelIndex = parallelIndices[innerIndex];
 if (!std::isnan(surplus[innerIndex])) {
 double partialSurplus = 0.0;
 if (dataIndexBase + parallelIndex < end_index_data
 && parallelIndex < X86PERFECT_PARALLEL_DATA_POINTS) {
 partialSurplus = phiEval[innerIndex] * alpha[dataIndexBase + parallelIndex];

 size_t localIndexFlat = indexFlat[innerIndex];

 //TODO: remove this
 // if (localIndexFlat == 0) {
 //   cout << "partial: " << partialSurplus << endl;
 // }
 // cout << "parital: " << partialSurplus << endl;
 // cout << "arr: " << levelArrayContinuous[localIndexFlat] << endl;
 #pragma omp atomic
 levelArrayContinuous[localIndexFlat] += partialSurplus;
 }
 nextIterationToRecalcReferences[parallelIndex] = subspace.nextDiff;
 // levelIndices[parallelIndex] += this->subspaceSize;
 } else {
 #if X86PERFECT_ENABLE_SUBSPACE_SKIPPING == 1
 //skip to next relevant subspace
 nextIterationToRecalcReferences[parallelIndex] = subspace.jumpDiff;
 levelIndices[parallelIndex] = subspace.nextSubspaceIndex;
 #else
 nextIterationToRecalcReferences[parallelIndex] = subspace.nextDiff;
 // levelIndices[parallelIndex] += this->subspaceSize;
 #endif
 }
 } // end innerIndex
 } // end parallel
 }
 */
