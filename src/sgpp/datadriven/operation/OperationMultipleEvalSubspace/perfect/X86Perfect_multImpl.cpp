#include "../../OperationMultipleEvalSubspace/perfect/X86Perfect.hpp"

void X86Perfect::multImpl(sg::base::DataVector &alpha, sg::base::DataVector &result, const size_t start_index_data,
        const size_t end_index_data) {

    size_t tid = omp_get_thread_num();
    if (tid == 0) {
        this->setCoefficients(result);
    }

#pragma omp barrier

    size_t dim = this->dataset->getNcols();
    const double * const datasetPtr = this->dataset->getPointer();

    size_t totalThreadNumber = X86PERFECT_PARALLEL_DATA_POINTS + X86PERFECT_VEC_PADDING;

    double *evalIndexValuesAll = new double[(dim + 1) * totalThreadNumber];
    for (size_t i = 0; i < (dim + 1) * totalThreadNumber; i++) {
        evalIndexValuesAll[i] = 1.0;
    }

    //for faster index flattening
    uint32_t *intermediatesAll = new uint32_t[(dim + 1) * totalThreadNumber];
    for (size_t i = 0; i < (dim + 1) * totalThreadNumber; i++) {
        intermediatesAll[i] = 0.0;
    }

    size_t validIndices[X86PERFECT_PARALLEL_DATA_POINTS + X86PERFECT_VEC_PADDING];
    size_t validIndicesCount;

    size_t levelIndices[X86PERFECT_PARALLEL_DATA_POINTS + X86PERFECT_VEC_PADDING];
    size_t nextIterationToRecalcReferences[X86PERFECT_PARALLEL_DATA_POINTS + X86PERFECT_VEC_PADDING];

    double *listSubspace = new double[this->maxGridPointsOnLevel * 2];
    for (size_t i = 0; i < this->maxGridPointsOnLevel * 2; i++) {
        listSubspace[i] = std::numeric_limits<double>::quiet_NaN();
    }

//    uint64_t jumpCount = 0;
//    uint64_t jumpDistance = 0;
//    uint64_t evaluationCounter = 0;

    for (size_t dataIndexBase = start_index_data; dataIndexBase < end_index_data; dataIndexBase +=
    X86PERFECT_PARALLEL_DATA_POINTS) {

        for (size_t i = 0; i < totalThreadNumber; i++) {
            levelIndices[i] = 0.0;
            nextIterationToRecalcReferences[i] = 0;
        }

        for (size_t subspaceIndex = 0; subspaceIndex < subspaceCount; subspaceIndex++) {
            X86PerfectSubspaceNode &subspace = this->allSubspaceNodes[subspaceIndex];

            //prepare the subspace array for a list type subspace
            if (subspace.type == X86PerfectSubspaceNode::SubspaceType::LIST) {

                //fill with surplusses
                for (ListSubspaceTuple &tuple : subspace.indexFlatSurplusPairs) {
                    //accumulator that are later added to the global surplusses
                    //listSubspace[tuple.first] = 0.0;
                    listSubspace[tuple.indexFlat] = 0.0;
                    listSubspace[tuple.indexFlat + 1] = tuple.isLeaf;
                }
            }

            validIndicesCount = 0;
            for (size_t parallelIndex = 0; parallelIndex < X86PERFECT_PARALLEL_DATA_POINTS; parallelIndex++) {
                size_t parallelLevelIndex = levelIndices[parallelIndex];

                if (parallelLevelIndex == subspaceIndex) {
                    validIndices[validIndicesCount] = parallelIndex;
                    validIndicesCount += 1;
                }
            }

            //padding for up to vector size, no padding required if all data tuples participate as
            //the number of data points is a multiple of the vector width
            size_t paddingSize = min((int) (validIndicesCount + X86PERFECT_VEC_PADDING),
            X86PERFECT_PARALLEL_DATA_POINTS + X86PERFECT_VEC_PADDING);
            for (size_t i = validIndicesCount; i < paddingSize; i++) {
                size_t threadId = X86PERFECT_PARALLEL_DATA_POINTS + (i - validIndicesCount);
                validIndices[i] = threadId;
                levelIndices[threadId] = 0;
                nextIterationToRecalcReferences[threadId] = 0;
                double *evalIndexValues = evalIndexValuesAll + (dim + 1) * threadId;

                //for faster index flattening, last element is for padding
                uint32_t *intermediates = intermediatesAll + (dim + 1) * threadId;
                for (size_t j = 0; j < dim; j++) {
                    evalIndexValues[j] = 1.0;
                    intermediates[j] = 0;
                }
            }

            if (subspace.type == X86PerfectSubspaceNode::SubspaceType::ARRAY) {

                //lock the current subspace, so that no atomic writes are necessary
                subspace.lockSubspace();

                uncachedMultInner(dim, datasetPtr, alpha, dataIndexBase, end_index_data, subspace, validIndicesCount,
                        validIndices, levelIndices, nextIterationToRecalcReferences, evalIndexValuesAll,
                        intermediatesAll);
                //unlocks the subspace lock for ARRAY and BLUEPRINT type subspaces
                subspace.unlockSubspace();

            } else if (subspace.type == X86PerfectSubspaceNode::SubspaceType::LIST) {

                size_t nextIterationToRecalc = nextIterationToRecalcReferences[validIndices[0]];

                listMultInner(dim, datasetPtr, alpha, dataIndexBase, end_index_data, subspace, listSubspace,
                        validIndicesCount, validIndices, levelIndices, nextIterationToRecalcReferences,
                        nextIterationToRecalc, evalIndexValuesAll, intermediatesAll);

                //write results into the global surplus array
                if (subspace.type == X86PerfectSubspaceNode::SubspaceType::LIST) {
                    //for (pair<uint32_t, double> &tuple : subspace.indexFlatSurplusPairs) {
                    for (ListSubspaceTuple &tuple : subspace.indexFlatSurplusPairs) {
                        //if (listSubspace[tuple.first] != 0.0) {
                        if (listSubspace[tuple.indexFlat] != 0.0) {
#pragma omp atomic
                            //tuple.second += listSubspace[tuple.first];
                            tuple.surplus += listSubspace[tuple.indexFlat];
                        }
                        //listSubspace[tuple.first] = std::numeric_limits<double>::quiet_NaN();
                        listSubspace[tuple.indexFlat] = std::numeric_limits<double>::quiet_NaN();
                        listSubspace[tuple.indexFlat + 1] = std::numeric_limits<double>::quiet_NaN();
                    }
                }

//                for (size_t validIndex = 0; validIndex < validIndicesCount; validIndex += X86PERFECT_VEC_PADDING) {
//                    for (size_t innerIndex = 0; innerIndex < 4; innerIndex++) {
//                        if (validIndex + innerIndex >= X86PERFECT_PARALLEL_DATA_POINTS) {
//                            break;
//                        }
//                        evaluationCounter += 1;
//                        if (levelIndices[validIndices[innerIndex]] != subspaceIndex + 1) {
//                            jumpCount += 1;
//                            jumpDistance += levelIndices[validIndices[innerIndex]] - subspaceIndex;
//                        }
//                    }
//                }
            }
        } // end iterate subspaces
    } // end iterate chunks

    delete evalIndexValuesAll;
    delete intermediatesAll;
    delete listSubspace;

//    if (jumpCount != 0) {
//        uint64_t totalEvaluations = jumpDistance + evaluationCounter;
//        cout << "avr. jump distance: " << (jumpDistance / jumpCount) << endl;
//        cout << "actual evaluations: " << evaluationCounter << endl;
//        cout << "skipped " << (jumpDistance) << " out of " << totalEvaluations << " total evaluations" << endl;
//        cout << "skipped " << ((double) (jumpDistance) / (double) totalEvaluations) * 100 << "% evaluations" << endl;
//    } else {
//        cout << "no jumps" << endl;
//    }
#pragma omp barrier

    if (tid == 0) {
        this->unflatten(result);
    }
}
