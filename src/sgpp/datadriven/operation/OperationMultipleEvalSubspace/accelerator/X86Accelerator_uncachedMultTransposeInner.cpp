#include "../../OperationMultipleEvalSubspace/accelerator/X86Accelerator.hpp"

void X86Accelerator::uncachedMultTransposeInner(size_t dim, const double * const datasetPtr, size_t dataIndexBase,
        size_t end_index_data, X86AcceleratorSubspaceNode &subspace, double *levelArrayContinuous,
        size_t validIndicesCount, size_t *validIndices, size_t *levelIndices, // size_t *nextIterationToRecalcReferences,
        double *componentResults, double *evalIndexValuesAll, uint32_t *intermediatesAll) {

    for (size_t validIndex = 0; validIndex < validIndicesCount; validIndex += X86ACCELERATOR_VEC_PADDING) {
    //for (size_t validIndex = 0; validIndex < validIndicesCount; validIndex += 4) {
        size_t parallelIndices[8];
        parallelIndices[0] = validIndices[validIndex];
        parallelIndices[1] = validIndices[validIndex + 1];
        parallelIndices[2] = validIndices[validIndex + 2];
        parallelIndices[3] = validIndices[validIndex + 3];

        parallelIndices[4] = validIndices[validIndex + 4];
        parallelIndices[5] = validIndices[validIndex + 5];
        parallelIndices[6] = validIndices[validIndex + 6];
        parallelIndices[7] = validIndices[validIndex + 7];

        //cout << "firstRound: " << firstRound << endl;

//        size_t nextIterationToRecalc = nextIterationToRecalcReferences[parallelIndices[0]];

#if X86ACCELERATOR_ENABLE_PARTIAL_RESULT_REUSAGE == 1
        size_t nextIterationToRecalc = subspace.arriveDiff;
#else
        size_t nextIterationToRecalc = 0;
#endif
//        if (nextIterationToRecalc != subspace.arriveDiff) {
//            cout << "your are a bit wrong today, arrive:" << subspace.arriveDiff << " ref: " << nextIterationToRecalc << endl;
//
//        }

        const double * const dataTuplePtr[8] = {
                datasetPtr + (dataIndexBase + parallelIndices[0]) * dim,
                datasetPtr + (dataIndexBase + parallelIndices[1]) * dim,
                datasetPtr + (dataIndexBase + parallelIndices[2]) * dim,
                datasetPtr + (dataIndexBase + parallelIndices[3]) * dim,

                datasetPtr + (dataIndexBase + parallelIndices[4]) * dim,
                datasetPtr + (dataIndexBase + parallelIndices[5]) * dim,
                datasetPtr + (dataIndexBase + parallelIndices[6]) * dim,
                datasetPtr + (dataIndexBase + parallelIndices[7]) * dim
 };

        double *evalIndexValues[8];
        evalIndexValues[0] = evalIndexValuesAll + (dim + 1) * parallelIndices[0];
        evalIndexValues[1] = evalIndexValuesAll + (dim + 1) * parallelIndices[1];
        evalIndexValues[2] = evalIndexValuesAll + (dim + 1) * parallelIndices[2];
        evalIndexValues[3] = evalIndexValuesAll + (dim + 1) * parallelIndices[3];

        evalIndexValues[4] = evalIndexValuesAll + (dim + 1) * parallelIndices[4];
        evalIndexValues[5] = evalIndexValuesAll + (dim + 1) * parallelIndices[5];
        evalIndexValues[6] = evalIndexValuesAll + (dim + 1) * parallelIndices[6];
        evalIndexValues[7] = evalIndexValuesAll + (dim + 1) * parallelIndices[7];

        //for faster index flattening, last element is for padding
        uint32_t *intermediates[8];
        intermediates[0] = intermediatesAll + (dim + 1) * parallelIndices[0];
        intermediates[1] = intermediatesAll + (dim + 1) * parallelIndices[1];
        intermediates[2] = intermediatesAll + (dim + 1) * parallelIndices[2];
        intermediates[3] = intermediatesAll + (dim + 1) * parallelIndices[3];

        intermediates[4] = intermediatesAll + (dim + 1) * parallelIndices[4];
        intermediates[5] = intermediatesAll + (dim + 1) * parallelIndices[5];
        intermediates[6] = intermediatesAll + (dim + 1) * parallelIndices[6];
        intermediates[7] = intermediatesAll + (dim + 1) * parallelIndices[7];

        uint32_t indexFlat[8];
        double phiEval[8];

        X86Accelerator::calculateIndexAccelerator(dim, nextIterationToRecalc, dataTuplePtr, subspace.hInverse, intermediates,
                evalIndexValues, indexFlat, phiEval);

        double surplus[8];
        surplus[0] = levelArrayContinuous[indexFlat[0]];
        surplus[1] = levelArrayContinuous[indexFlat[1]];
        surplus[2] = levelArrayContinuous[indexFlat[2]];
        surplus[3] = levelArrayContinuous[indexFlat[3]];

        surplus[4] = levelArrayContinuous[indexFlat[4]];
        surplus[5] = levelArrayContinuous[indexFlat[5]];
        surplus[6] = levelArrayContinuous[indexFlat[6]];
        surplus[7] = levelArrayContinuous[indexFlat[7]];

        for (size_t innerIndex = 0; innerIndex < 8; innerIndex++) {
            if (!std::isnan(surplus[innerIndex])) {
                //cout << parallelIndices[innerIndex] << " -> phi: " << phiEval[innerIndex] << " surplus: " << surplus[innerIndex] << " = " << phiEval[innerIndex] * surplus[innerIndex] << endl;
                /*if (nextRoundResults[innerIndex] != phiEval[innerIndex] * surplus[innerIndex] && firstRound == false) {
                    cout << "differ validindex: " << parallelIndices[innerIndex] << " => ref: " << phiEval[innerIndex] * surplus[innerIndex] << " new: " << nextRoundResults[innerIndex] << endl;

                }*/
                componentResults[parallelIndices[innerIndex]] += phiEval[innerIndex] * surplus[innerIndex];
                //nextIterationToRecalcReferences[parallelIndices[innerIndex]] = subspace.nextDiff;
                levelIndices[parallelIndices[innerIndex]] += 1;
            } else {

#if X86ACCELERATOR_ENABLE_SUBSPACE_SKIPPING == 1	      
                //skip to next relevant subspace
                //nextIterationToRecalcReferences[parallelIndices[innerIndex]] = subspace.jumpDiff;
                levelIndices[parallelIndices[innerIndex]] = subspace.jumpTargetIndex;
#else
                //nextIterationToRecalcReferences[parallelIndices[innerIndex]] = subspace.nextDiff;
                levelIndices[parallelIndices[innerIndex]] += 1;
#endif
            }
        }
    } //end X86ACCELERATOR_PARALLEL_DATA_POINTS
}
