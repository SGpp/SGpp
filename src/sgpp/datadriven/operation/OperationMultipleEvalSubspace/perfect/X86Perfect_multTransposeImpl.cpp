#include <iomanip>
#include <algorithm>
#include "../../OperationMultipleEvalSubspace/perfect/X86Perfect.hpp"

//#include "X86PerfectTypes.hpp"

void X86Perfect::multTransposeImpl(sg::base::DataVector &source, sg::base::DataVector &result,
        const size_t start_index_data, const size_t end_index_data) {
    size_t tid = omp_get_thread_num();
    if (tid == 0) {
        this->setCoefficients(source);
    }

#pragma omp barrier

    size_t dim = dataset->getNcols();
    double * datasetPtr = dataset->getPointer();

    size_t totalThreadNumber = X86PERFECT_PARALLEL_DATA_POINTS + X86PERFECT_VEC_PADDING;

    double *evalIndexValuesAll = new double[(dim + 1) * totalThreadNumber];
    for (size_t i = 0; i < (dim + 1) * totalThreadNumber; i++) {
        evalIndexValuesAll[i] = 1.0;
    }

    //for faster index flattening, last element is for padding
    uint32_t *intermediatesAll = new uint32_t[(dim + 1) * totalThreadNumber];
    for (size_t i = 0; i < (dim + 1) * totalThreadNumber; i++) {
        intermediatesAll[i] = 0.0;
    }

    size_t validIndices[X86PERFECT_PARALLEL_DATA_POINTS + X86PERFECT_VEC_PADDING];
    size_t validIndicesCount;

    double componentResults[X86PERFECT_PARALLEL_DATA_POINTS + X86PERFECT_VEC_PADDING];
    size_t levelIndices[X86PERFECT_PARALLEL_DATA_POINTS + X86PERFECT_VEC_PADDING];
    size_t nextIterationToRecalcReferences[X86PERFECT_PARALLEL_DATA_POINTS + X86PERFECT_VEC_PADDING];

    double *listSubspace = new double[this->maxGridPointsOnLevel * 2];
    for (size_t i = 0; i < this->maxGridPointsOnLevel * 2; i++) {
        listSubspace[i] = std::numeric_limits<double>::quiet_NaN();
    }

    //process the next chunk of data tuples in parallel
    for (size_t dataIndexBase = start_index_data; dataIndexBase < end_index_data; dataIndexBase +=
    X86PERFECT_PARALLEL_DATA_POINTS) {

        for (size_t i = 0; i < totalThreadNumber; i++) {
            levelIndices[i] = 0.0;
            componentResults[i] = 0.0;
            nextIterationToRecalcReferences[i] = 0;
        }

        for (size_t subspaceIndex = 0; subspaceIndex < subspaceCount; subspaceIndex++) {
            X86PerfectSubspaceNode &subspace = this->allSubspaceNodes[subspaceIndex];

            double *levelArrayContinuous = nullptr;
            //prepare the subspace array for a list type subspace
            if (subspace.type == X86PerfectSubspaceNode::SubspaceType::LIST) {
                //cout << "list subspace" << endl;
                //fill with surplusses
                //for (pair<uint32_t, double> tuple : subspace.indexFlatSurplusPairs) {
                for (ListSubspaceTuple tuple : subspace.indexFlatSurplusPairs) {
                    //actual values are utilized, but only read
                    //listSubspace[tuple.first] = tuple.second;
                    listSubspace[tuple.indexFlat] = tuple.surplus;
                    listSubspace[tuple.indexFlat + 1] = tuple.isLeaf ? 1.0 : 0.0;
                }
                levelArrayContinuous = listSubspace;
            } else {
                //cout << "array subspace" << endl;
                levelArrayContinuous = subspace.subspaceArray;
            }

            validIndicesCount = 0;
            for (size_t parallelIndex = 0; parallelIndex < X86PERFECT_PARALLEL_DATA_POINTS; parallelIndex++) {
                size_t parallelLevelIndex = levelIndices[parallelIndex];
                if (parallelLevelIndex == subspaceIndex) {
                    validIndices[validIndicesCount] = parallelIndex;
                    validIndicesCount += 1;
                }
            }

            size_t paddingSize = min((int) (validIndicesCount + 4),
            X86PERFECT_PARALLEL_DATA_POINTS + X86PERFECT_VEC_PADDING);
            for (size_t i = validIndicesCount; i < paddingSize; i++) {
                size_t threadId = X86PERFECT_PARALLEL_DATA_POINTS + (i - validIndicesCount);
                validIndices[i] = threadId;
                componentResults[threadId] = 0.0;
                levelIndices[threadId] = 0;
                nextIterationToRecalcReferences[threadId] = 0;
                double *evalIndexValues = evalIndexValuesAll + (dim + 1) * threadId;

                //for faster index flattening, last element is for padding
                uint32_t *intermediates = intermediatesAll + (dim + 1) * threadId;
                for (size_t j = 0; j < dim; j++) {
                    evalIndexValues[j] = 1.0;
                    intermediates[j] = 0.0;
                }
            }

            uncachedMultTransposeInner(dim, datasetPtr, dataIndexBase, end_index_data, subspace, levelArrayContinuous,
                    validIndicesCount, validIndices, levelIndices, nextIterationToRecalcReferences, componentResults,
                    evalIndexValuesAll, intermediatesAll);

            if (subspace.type == X86PerfectSubspaceNode::SubspaceType::LIST) {
                //for (pair<uint32_t, double> &tuple : subspace.indexFlatSurplusPairs) {
                for (ListSubspaceTuple &tuple : subspace.indexFlatSurplusPairs) {
                    //listSubspace[tuple.first] = std::numeric_limits<double>::quiet_NaN();
                    listSubspace[tuple.indexFlat] = std::numeric_limits<double>::quiet_NaN();
                    listSubspace[tuple.indexFlat + 1] = std::numeric_limits<double>::quiet_NaN();
                }
            }
        } // end iterate grid

        for (size_t parallelIndex = 0; parallelIndex < X86PERFECT_PARALLEL_DATA_POINTS; parallelIndex++) {
            size_t dataIndex = dataIndexBase + parallelIndex;
            result.set(dataIndex, componentResults[parallelIndex]);
        }
    } // end iterate data chunks

    delete evalIndexValuesAll;
    delete intermediatesAll;
    delete listSubspace;
}
