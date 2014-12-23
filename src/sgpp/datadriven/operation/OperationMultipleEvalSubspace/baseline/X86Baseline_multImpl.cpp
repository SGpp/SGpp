#include "../../OperationMultipleEvalSubspace/baseline/X86Baseline.hpp"

void X86Baseline::multImpl(sg::base::DataVector &alpha, sg::base::DataVector &result,
//const size_t start_index_grid,
//const size_t end_index_grid,
        const size_t start_index_data, const size_t end_index_data) {

    size_t tid = omp_get_thread_num();
    if (tid == 0) {
        this->setCoefficients(result);
    }

#pragma omp barrier

    //cout << "multImpl" << endl;

    /*int totalActualEvaluated = 0;
     int totalNotNaNs = 0;
     int totalJumpDistance = 0;
     int totalJumps = 0;*/

    size_t dim = this->dataset->getNcols();
    double *datasetPtr = this->dataset->getPointer();

    //DataVector dataTuple(dim);
    //double *dataTuplePtr = dataTuple.getPointer();

    //TODO: replace by "new"
    //size_t *indexPtr = (size_t *) malloc(sizeof(size_t) * dim);
    size_t *indexPtrAll = new size_t[(dim + 1) * X86BASELINE_PARALLEL_DATA_POINTS];
    for (size_t i = 0; i < (dim + 1) * X86BASELINE_PARALLEL_DATA_POINTS; i++) {
        indexPtrAll[i] = 0;
    }
    //double *evalIndexValues = (double *) malloc(sizeof(double) * (dim + 1));
    double *evalIndexValuesAll = new double[(dim + 1) * X86BASELINE_PARALLEL_DATA_POINTS];
    for (size_t i = 0; i < (dim + 1) * X86BASELINE_PARALLEL_DATA_POINTS; i++) {
        evalIndexValuesAll[i] = 1.0;
    }

//  double **flatLevels = this->flatLevels;

    //for faster index flattening
    //double *intermediates = (double *) malloc(sizeof(double) * (dim + 1));
    double *intermediatesAll = new double[(dim + 1) * X86BASELINE_PARALLEL_DATA_POINTS];
    for (size_t i = 0; i < (dim + 1) * X86BASELINE_PARALLEL_DATA_POINTS; i++) {
        intermediatesAll[i] = 0.0;
    }

    //double maxIndex = sizeof(allSubspaces) / (subspaceSize * sizeof(size_t));
    double maxIndex = subspaceCount * subspaceSize;
    //for (size_t dataIndex = start_index_data; dataIndex < end_index_data; dataIndex++) {
    for (size_t dataIndexBase = start_index_data; dataIndexBase < end_index_data; dataIndexBase +=
            X86BASELINE_PARALLEL_DATA_POINTS) {

        //evalIndexValues[0] = 1.0;
        //intermediates[0] = 0.0;
        //int actualEvaluated = 0;
        //dataset->getRow(dataIndex, dataTuple);
        //size_t levelIndex = 0;
        //size_t nextIterationToRecalc = 0; //all index components are recalculated

        //size_t levelIndex = 0;
        size_t levelIndices[X86BASELINE_PARALLEL_DATA_POINTS];
        //size_t nextIterationToRecalc = 0; //all index components are recalculated
        size_t nextIterationToRecalcReferences[X86BASELINE_PARALLEL_DATA_POINTS]; //all index components are recalculated
        for (size_t i = 0; i < X86BASELINE_PARALLEL_DATA_POINTS; i++) {
            levelIndices[i] = 0.0;
            //componentResults[i] = 0.0;
            nextIterationToRecalcReferences[i] = 0;
        }

        //while (levelIndex < maxIndex) {
        for (size_t allLevelIndex = 0; allLevelIndex < maxIndex; allLevelIndex += subspaceSize) {

            size_t *hInversePtr = allSubspaces + allLevelIndex + dim;
            //size_t *levelPtr = allSubspaces + levelIndex;
            //size_t levelFlat = *(allSubspaces + allLevelIndex + (2 * dim) + 1);
            size_t linearLevelArray = *(allSubspaces + allLevelIndex + (2 * dim) + 1);
            double *levelArray = &(this->allSurplusses[linearLevelArray]);

            for (size_t parallelIndex = 0; parallelIndex < X86BASELINE_PARALLEL_DATA_POINTS; parallelIndex++) {
                //cout << "data points index: " << parallelIndex << endl;
                size_t levelIndex = levelIndices[parallelIndex];
                if (levelIndex != allLevelIndex) {
                    continue; // current element skips subspace
                }
                //actualEvaluated += 1;
                size_t nextIterationToRecalc = nextIterationToRecalcReferences[parallelIndex];
                double *dataTuplePtr = datasetPtr + (dataIndexBase + parallelIndex) * dim;
                size_t *indexPtr = indexPtrAll + (dim + 1) * parallelIndex;
                double *evalIndexValues = evalIndexValuesAll + (dim + 1) * parallelIndex;
                //for faster index flattening, last element is for padding
                double *intermediates = intermediatesAll + (dim + 1) * parallelIndex;

                // calculate index
                for (size_t i = nextIterationToRecalc; i < dim; i++) {
                    double unadjusted = dataTuplePtr[i] * hInversePtr[i];
                    indexPtr[i] = calculateIndexComponent(dim, unadjusted);
                }

                int indexFlat = X86Baseline::flattenIndex(intermediates, dim, hInversePtr, indexPtr,
                        nextIterationToRecalc);
                //double *levelArray = flatLevels[levelFlat];
                double surplus = levelArray[indexFlat];

                if (!std::isnan(surplus)) {
                    //totalNotNaNs += 1;
                    double phiEval = evalIndexValues[nextIterationToRecalc];
                    //prepare the values for the individual components
                    for (size_t i = nextIterationToRecalc; i < dim; i++) {
                        double phi1DEval = hInversePtr[i] * dataTuplePtr[i] - indexPtr[i];
                        phi1DEval = max(0.0, 1.0 - abs(phi1DEval));
                        phiEval *= phi1DEval;
                        evalIndexValues[i + 1] = phiEval;
                    }

                    double partialSurplus = phiEval * alpha[dataIndexBase + parallelIndex];
                    //componentResults[parallelIndex] += partialSurplus;

#pragma omp atomic
                    levelArray[indexFlat] += partialSurplus;

                    //nextIterationToRecalc = allSubspaces[levelIndex + (2 * dim) + 2];
                    nextIterationToRecalcReferences[parallelIndex] = allSubspaces[levelIndex + (2 * dim) + 2];
                    //levelIndex += subspaceSize;
                    levelIndices[parallelIndex] += subspaceSize;
                } else {
                    //size_t old = levelIndex;

#if X86BASELINE_ENABLE_SUBSPACE_SKIPPING == 1
                    //skip to next relevant subspace
                    //nextIterationToRecalc = allSubspaces[levelIndex + (2 * dim) + 3];
                    nextIterationToRecalcReferences[parallelIndex] = allSubspaces[levelIndex + (2 * dim) + 3];
                    //levelIndex = allSubspaces[levelIndex + 2 * dim];
                    levelIndices[parallelIndex] = allSubspaces[levelIndex + 2 * dim];
#else
                    nextIterationToRecalcReferences[parallelIndex] = allSubspaces[levelIndex + (2 * dim) + 2];
                    levelIndices[parallelIndex] += subspaceSize;
#endif

                    //nextIterationToRecalc = allSubspaces[levelIndex + (2 * dim) + 2];
                    //levelIndex += subspaceSize;
                    //totalJumpDistance += (levelIndex - old) / subspaceSize;
                    //totalJumps += 1;
                }
            } // end parallel

        } // end iterate subspaces
          //totalActualEvaluated += actualEvaluated;

    } // end iterate chunks

    /*if (omp_get_thread_num() == 0) {
     if (totalActualEvaluated > 0) {
     cout << "totalActualEvaluated: " << totalActualEvaluated << endl;
     cout << "avr. actual evaluated: " << ((double) totalActualEvaluated / ((double) end_index_data - start_index_data)) << endl;
     } else {
     cout << "avr. actual evaluated: empty range!" << endl;
     }
     if (totalNotNaNs > 0) {
     cout << "totalNotNaNs: " << totalNotNaNs << endl;
     cout << "avr. totalNotNaNs: " << ((double) totalNotNaNs / ((double) end_index_data - start_index_data)) << endl;
     } else {
     cout << "avr. totalNotNaNs: empty range!" << endl;
     }
     if (totalJumps > 0) {
     cout << "avr. jump distance: " << (totalJumpDistance / totalJumps) << endl;
     } else {
     cout << "no jumps!" << endl;
     }
     }
     */

    delete indexPtrAll;
    delete evalIndexValuesAll;
    delete intermediatesAll;

#pragma omp barrier

    if (tid == 0) {
        this->unflatten(result);
    }

}
