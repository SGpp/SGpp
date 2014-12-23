#include "../../OperationMultipleEvalSubspace/accelerator/X86Accelerator.hpp"

X86Accelerator::X86Accelerator(GridStorage* storage, DataMatrix* dataset, size_t dim) {
    this->storage = storage;
    this->dataset = dataset;
    this->dim = dim;

    this->subspaceSize = (2 * this->dim) + 8;
    this->maxGridPointsOnLevel = 0;

#ifdef X86ACCELERATOR_WRITE_STATS
    string prefix("results/data/stats_");
    string fileName(X86ACCELERATOR_WRITE_STATS);
    this->statsFile.open(prefix + fileName, ios::out);;

    this->statsFile << "# name: " << X86ACCELERATOR_WRITE_STATS_NAME << endl;
    this->statsFile << "refinementStep & ";
    this->statsFile << "nonVirtualGridPoints & ";
    this->statsFile << "totalRegularGridPoints & ";
    this->statsFile << "actualGridPoints & ";
    this->statsFile << "largestArraySubspace & ";
    this->statsFile << "largestListSubspace & ";
    this->statsFile << "numberOfListSubspaces & ";
    this->statsFile << "subspaceCount & ";
    this->statsFile << "avrPointsPerSubspace & ";
    this->statsFile << "memoryEstimate & ";
    this->statsFile << "memoryEfficiency";
    this->statsFile << endl;

#endif
}

X86Accelerator::~X86Accelerator() {
#ifdef X86ACCELERATOR_WRITE_STATS
    this->statsFile.close();
#endif
}

void X86Accelerator::prepare() {
    this->allLevelsIndexMap.clear();
    this->allSubspaceNodes.clear();
    this->prepareSubspaceIterator();
}

void X86Accelerator::setCoefficients(DataVector &surplusVector) {
    DataVector level(dim);
    DataVector maxIndex(dim);
    DataVector index(dim);

    //TODO: use appropriate types here
    unsigned int curLevel;
    unsigned int curIndex;
    for (size_t gridIndex = 0; gridIndex < this->storage->size(); gridIndex++) {
        sg::base::GridIndex *point = this->storage->get(gridIndex);
        for (size_t d = 0; d < this->dim; d++) {
            point->get(d, curLevel, curIndex);
            level.set(d, curLevel);
            index.set(d, curIndex);
            maxIndex.set(d, 1 << curLevel);
        }

        this->setSurplus(level, maxIndex, index, surplusVector.get(gridIndex));
    }
}

//writes a result vector in the order of the points in the grid storage
void X86Accelerator::unflatten(DataVector &result) {
    DataVector level(dim);
    DataVector maxIndex(dim);
    DataVector index(dim);

    unsigned int curLevel;
    unsigned int curIndex;
    for (size_t gridIndex = 0; gridIndex < this->storage->size(); gridIndex++) {
        sg::base::GridIndex *point = this->storage->get(gridIndex);
        for (size_t d = 0; d < this->dim; d++) {
            point->get(d, curLevel, curIndex);
            level.set(d, curLevel);
            index.set(d, curIndex);
            maxIndex.set(d, 1 << curLevel);
        }
        double surplus;
        bool isVirtual;
        this->getSurplus(level, maxIndex, index, surplus, isVirtual);

        result.set(gridIndex, surplus);
    }
}

void X86Accelerator::setSurplus(DataVector &level, DataVector &maxIndices, DataVector &index, double value) {
    uint32_t levelFlat = this->flattenLevel(this->dim, this->maxLevel, level);
    uint32_t indexFlat = this->flattenIndex(this->dim, maxIndices, index);
    uint32_t subspaceIndex = this->allLevelsIndexMap.find(levelFlat)->second;
    X86AcceleratorSubspaceNode &subspace = this->allSubspaceNodes[subspaceIndex];
    subspace.setSurplus(indexFlat, value);
}

void X86Accelerator::getSurplus(DataVector &level, DataVector &maxIndices, DataVector &index, double &value,
        bool &isVirtual) {
    uint32_t levelFlat = this->flattenLevel(this->dim, this->maxLevel, level);
    uint32_t indexFlat = this->flattenIndex(this->dim, maxIndices, index);
    uint32_t subspaceIndex = this->allLevelsIndexMap.find(levelFlat)->second;
    X86AcceleratorSubspaceNode &subspace = this->allSubspaceNodes[subspaceIndex];
    value = subspace.getSurplus(indexFlat);
    if (std::isnan(value)) {
        isVirtual = true;
    } else {
        isVirtual = false;
    }
}

uint32_t X86Accelerator::flattenIndex(size_t dim, DataVector &maxIndices, DataVector &index) {
    uint32_t indexFlat = index.get(0);
    indexFlat >>= 1;
    for (size_t i = 1; i < dim; i++) {
        int actualDirectionGridPoints = maxIndices.get(i);
        actualDirectionGridPoints >>= 1;
        indexFlat *= actualDirectionGridPoints;
        uint32_t actualIndex = index.get(i);
        actualIndex >>= 1; //divide index by 2, skip even indices
        indexFlat += actualIndex;
    }
    return indexFlat;
}

uint32_t X86Accelerator::flattenLevel(size_t dim, size_t maxLevel, DataVector &level) {
    uint32_t levelFlat = 0;
    levelFlat += level.get(dim - 1);
    // loop terminates at -1
    for (int i = dim - 2; i >= 0; i--) {
        levelFlat *= maxLevel;
        levelFlat += level.get(i);
    }
    return levelFlat;
}

void X86Accelerator::padDataset() {
    size_t chunkSize = X86ACCELERATOR_PARALLEL_DATA_POINTS;

    // Assure that data has a even number of instances -> padding might be needed
    size_t remainder = this->dataset->getNrows() % chunkSize;
    size_t loopCount = chunkSize - remainder;

    if (loopCount != chunkSize) {
        sg::base::DataVector lastRow(this->dataset->getNcols());
        size_t oldSize = this->dataset->getNrows();
        this->dataset->getRow(this->dataset->getNrows() - 1, lastRow);
        this->dataset->resize(this->dataset->getNrows() + loopCount);

        for (size_t i = 0; i < loopCount; i++) {
            this->dataset->setRow(oldSize + i, lastRow);
        }
    }
    //initialize padding area - hacky
    this->dataset->addSize(X86ACCELERATOR_VEC_PADDING * 2);
    for (size_t i = this->dataset->getNrows(); i < this->dataset->getNrows() + this->dataset->getUnused(); i++) {
        for (size_t j = 0; j < this->dataset->getNcols(); j++) {
            this->dataset->set(i, j, 0.0);
        }
    }
}
