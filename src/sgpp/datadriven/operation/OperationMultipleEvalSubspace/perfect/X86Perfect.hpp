#pragma once

#include <iostream>
#include <vector>
#include <map>
#include "omp.h"
#include <immintrin.h>

#include "../../OperationMultipleEvalSubspace/perfect/X86PerfectParameters.hpp"
#include "../../OperationMultipleEvalSubspace/perfect/X86PerfectSubspaceNode.hpp"
#include "FastKernelBase.hpp"

using namespace std;
using namespace sg::base;
using namespace sg::parallel;

namespace sg {
namespace parallel {

class X86Perfect: public FastKernelBase {
private:

    size_t subspaceSize = -1;

    size_t maxGridPointsOnLevel;

    map<uint32_t, uint32_t> allLevelsIndexMap;

    DataMatrix *dataset = nullptr;
    size_t dim = -1;
    size_t maxLevel = 0;

    vector<X86PerfectSubspaceNode> allSubspaceNodes;
    size_t subspaceCount = -1;

    /// Pointer to the grid's gridstorage object
    sg::base::GridStorage* storage = nullptr;
    uint32_t totalRegularGridPoints = -1;

#ifdef X86PERFECT_WRITE_STATS
    size_t refinementStep = 0;
    ofstream statsFile;
    string csvSep = "& ";
#endif

    void prepareSubspaceIterator();

    void uncachedMultInner(size_t dim, const double * const datasetPtr, sg::base::DataVector &alpha,
            size_t dataIndexBase, size_t end_index_data, X86PerfectSubspaceNode &subspace, size_t validIndicesCount,
            size_t *validIndices, size_t *levelIndices, size_t *nextIterationToRecalcReferences,
            double *evalIndexValuesAll, uint32_t *intermediatesAll);

    void listMultInner(size_t dim, const double * const datasetPtr, sg::base::DataVector &alpha, size_t dataIndexBase,
            size_t end_index_data, X86PerfectSubspaceNode &subspace, double *levelArrayContinuous,
            size_t validIndicesCount, size_t *validIndices, size_t *levelIndices,
            size_t *nextIterationToRecalcReferences, size_t nextIterationToRecalc, double *evalIndexValuesAll,
            uint32_t *intermediatesAll);

    void uncachedMultTransposeInner(size_t dim, const double * const datasetPtr, size_t dataIndexBase,
            size_t end_index_data, X86PerfectSubspaceNode &subspace, double *levelArrayContinuous,
            size_t validIndicesCount, size_t *validIndices, size_t *levelIndices,
            size_t *nextIterationToRecalcReferences, double *componentResults, double *evalIndexValuesAll,
            uint32_t *intermediatesAll);

public:

#include "../../OperationMultipleEvalSubspace/perfect/X86Perfect_calculateIndexPerfect.hpp"

    static const KernelType kernelType = Standard;

    X86Perfect(GridStorage* storage, DataMatrix* dataset, size_t dim);

    ~X86Perfect();

    void prepare();

    void setCoefficients(DataVector &surplusVector);

    void unflatten(DataVector &result);

    static uint32_t flattenIndex(size_t dim, DataVector &maxIndices, DataVector &index);

    void setSurplus(DataVector &level, DataVector &maxIndices, DataVector &index, double value, bool isLeaf);

    void getSurplus(DataVector &level, DataVector &maxIndices, DataVector &index, double &value, bool &isVirtual);

    uint32_t flattenLevel(size_t dim, size_t maxLevel, DataVector &level);

    void multImpl(sg::base::DataVector &alpha, sg::base::DataVector &result, const size_t start_index_data,
            const size_t end_index_data);

    void multTransposeImpl(sg::base::DataVector &source, sg::base::DataVector &result, const size_t start_index_data,
            const size_t end_index_data);

    void padDataset();

    size_t getAlignment() {
        return X86PERFECT_PARALLEL_DATA_POINTS;
    }
};

}
}
