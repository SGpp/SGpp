#include "../../OperationMultipleEvalSubspace/perfect/X86Perfect.hpp"

void X86Perfect::uncachedMultTransposeInner(size_t dim, const double * const datasetPtr, size_t dataIndexBase,
        size_t end_index_data, X86PerfectSubspaceNode &subspace, double *levelArrayContinuous, size_t validIndicesCount,
        size_t *validIndices, size_t *levelIndices, size_t *nextIterationToRecalcReferences, double *componentResults,
        double *evalIndexValuesAll, uint32_t *intermediatesAll) {
    for (size_t validIndex = 0; validIndex < validIndicesCount; validIndex += 4) {

        size_t parallelIndices[4];
        parallelIndices[0] = validIndices[validIndex];
        parallelIndices[1] = validIndices[validIndex + 1];
        parallelIndices[2] = validIndices[validIndex + 2];
        parallelIndices[3] = validIndices[validIndex + 3];

        size_t nextIterationToRecalc = nextIterationToRecalcReferences[parallelIndices[0]];

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

        /*cout << "flatlevel: " << subspace.flatLevel << endl;
         cout << "grid points on level:" << subspace.gridPointsOnLevel << endl;
         cout << "if0: " << indexFlat[0] << endl;
         cout << "if1: " << indexFlat[1] << endl;
         cout << "if2: " << indexFlat[2] << endl;
         cout << "if3: " << indexFlat[3] << endl;*/

        for (size_t innerIndex = 0; innerIndex < 4; innerIndex++) {
            if (!std::isnan(surplus[innerIndex])) {
            componentResults[parallelIndices[innerIndex]] += phiEval[innerIndex] * surplus[innerIndex];
            }
            // leaf -> can jump, after jump I can continue at another non-existing point -> can jump
            if (isLeaf[innerIndex] || std::isnan(surplus[innerIndex])) {
                nextIterationToRecalcReferences[parallelIndices[innerIndex]] = subspace.jumpDiff;
                levelIndices[parallelIndices[innerIndex]] = subspace.jumpTargetIndex;
            } else {
                nextIterationToRecalcReferences[parallelIndices[innerIndex]] = subspace.nextDiff;
                levelIndices[parallelIndices[innerIndex]] += 1;
            }

            /*
             if (!std::isnan(surplus[innerIndex])) {
             componentResults[parallelIndices[innerIndex]] += phiEval[innerIndex] * surplus[innerIndex];
             nextIterationToRecalcReferences[parallelIndices[innerIndex]] = subspace.nextDiff;
             levelIndices[parallelIndices[innerIndex]] += 1;
             } else {
             #if X86PERFECT_ENABLE_SUBSPACE_SKIPPING == 1
             //skip to next relevant subspace
             nextIterationToRecalcReferences[parallelIndices[innerIndex]] = subspace.jumpDiff;
             levelIndices[parallelIndices[innerIndex]] = subspace.jumpTargetIndex;
             #else
             nextIterationToRecalcReferences[parallelIndices[innerIndex]] = subspace.nextDiff;
             levelIndices[parallelIndices[innerIndex]] += 1;
             #endif
             }
             */

        }
    } //end X86PERFECT_PARALLEL_DATA_POINTS
}
