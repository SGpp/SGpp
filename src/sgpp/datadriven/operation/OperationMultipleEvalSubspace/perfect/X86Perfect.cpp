#include "../../OperationMultipleEvalSubspace/perfect/X86Perfect.hpp"

X86Perfect::X86Perfect(GridStorage* storage, DataMatrix* dataset, size_t dim) {
    this->storage = storage;
    this->dataset = dataset;
    this->dim = dim;

    this->subspaceSize = (2 * this->dim) + 8;
    this->maxGridPointsOnLevel = 0;

#ifdef X86PERFECT_WRITE_STATS
    string prefix("results/data/stats_");
    string fileName(X86PERFECT_WRITE_STATS);
    this->statsFile.open(prefix + fileName, ios::out);;

    this->statsFile << "# name: " << X86PERFECT_WRITE_STATS_NAME << endl;
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

X86Perfect::~X86Perfect() {
#ifdef X86PERFECT_WRITE_STATS
    this->statsFile.close();
#endif
}

void X86Perfect::prepare() {
    this->allLevelsIndexMap.clear();
    this->allSubspaceNodes.clear();
    this->prepareSubspaceIterator();
}

void X86Perfect::setCoefficients(DataVector &surplusVector) {
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

        /*static uint32_t counter = 0;
         if (counter < 1000) {
         cout << "----------" << endl;
         cout << "is leaf: " << point->isLeaf() << endl;
         cout << "index: " << index.toString() << endl;
         counter++;
         }*/
        this->setSurplus(level, maxIndex, index, surplusVector.get(gridIndex), point->isLeaf());
    }
}

//writes a result vector in the order of the points in the grid storage
void X86Perfect::unflatten(DataVector &result) {
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

void X86Perfect::setSurplus(DataVector &level, DataVector &maxIndices, DataVector &index, double value, bool isLeaf) {
    uint32_t levelFlat = this->flattenLevel(this->dim, this->maxLevel, level);
    uint32_t indexFlat = this->flattenIndex(this->dim, maxIndices, index);
    uint32_t subspaceIndex = this->allLevelsIndexMap.find(levelFlat)->second;
    X86PerfectSubspaceNode &subspace = this->allSubspaceNodes[subspaceIndex];
    subspace.setSurplus(indexFlat, value, isLeaf);
}

void X86Perfect::getSurplus(DataVector &level, DataVector &maxIndices, DataVector &index, double &value,
        bool &isVirtual) {
    uint32_t levelFlat = this->flattenLevel(this->dim, this->maxLevel, level);
    uint32_t indexFlat = this->flattenIndex(this->dim, maxIndices, index);
    uint32_t subspaceIndex = this->allLevelsIndexMap.find(levelFlat)->second;
    X86PerfectSubspaceNode &subspace = this->allSubspaceNodes[subspaceIndex];
    value = subspace.getSurplus(indexFlat);
    if (std::isnan(value)) {
        isVirtual = true;
    } else {
        isVirtual = false;
    }
}

uint32_t X86Perfect::flattenIndex(size_t dim, DataVector &maxIndices, DataVector &index) {
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
    indexFlat <<= 1; //times 2 to represent surplus-isLeaf tuples
    return indexFlat;
}

uint32_t X86Perfect::flattenLevel(size_t dim, size_t maxLevel, DataVector &level) {
    uint32_t levelFlat = 0;
    levelFlat += level.get(dim - 1);
// loop terminates at -1
    for (int i = dim - 2; i >= 0; i--) {
        levelFlat *= maxLevel;
        levelFlat += level.get(i);
    }
    return levelFlat;
}

void X86Perfect::padDataset() {
    size_t chunkSize = X86PERFECT_PARALLEL_DATA_POINTS;

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
    this->dataset->addSize(X86PERFECT_VEC_PADDING * 2);
    for (size_t i = this->dataset->getNrows(); i < this->dataset->getNrows() + this->dataset->getUnused(); i++) {
        for (size_t j = 0; j < this->dataset->getNcols(); j++) {
            this->dataset->set(i, j, 0.0);
        }
    }
}
