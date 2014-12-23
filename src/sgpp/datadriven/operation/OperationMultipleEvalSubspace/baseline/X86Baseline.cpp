#include "../../OperationMultipleEvalSubspace/baseline/X86Baseline.hpp"
#include <limits>

namespace sg {
namespace parallel {
namespace x86baseline {

size_t problemDimForSorting;

int subspaceComparator(const void *first, const void *second) {
    for (size_t i = 0; i < problemDimForSorting; i++) {
        if (((size_t *) first)[i] >= ((size_t *) second)[i]) {
            if (((size_t *) first)[i] > ((size_t *) second)[i]) {
                return 1;
            }
        } else {
            return 0;
        }
    }
    return 0;
}

size_t getHighestDifferentIndex(size_t *base, size_t cur, size_t last) {
    //for (int i = problemDimForSorting - 1; i >= 0; i--) {
    for (size_t i = 0; i < problemDimForSorting; i++) {
        if (base[cur + (size_t) i] != base[last + (size_t) i]) {
            return i;
        }
    }
    return -1; //should never be returned
}

void printSubspace(size_t *base, size_t dim, size_t cur, int subspaceSize) {
    cout << "l: ";
    for (size_t j = 0; j < dim; j++) {
        if (j > 0) {
            cout << ", ";
        }
        cout << base[cur + j];
    }
    if (subspaceSize != -1) {
        cout << " next: ";
        cout << (base[cur + dim] / subspaceSize);
    }
}

}
}
}

using namespace sg::parallel::x86baseline;

X86Baseline::X86Baseline(GridStorage* storage, DataMatrix* dataset, size_t dim) {
    this->storage = storage;
    this->dataset = dataset;
    this->dim = dim;
}

X86Baseline::~X86Baseline() {
//  if (this->flatLevels != nullptr) {
//    DataVector level(this->dim);
//    for (size_t i = 0; i < this->subspaceCount; i++) {
//      //allLevels.getRow(i, level);
//      int flatLevel = this->allSubspaces[i * this->subspaceSize + (2 * dim) + 1];
//      //int flatLevel = this->flattenLevel(this->dim, this->maxLevel, level);
//      delete this->flatLevels[flatLevel];
//    }
//    delete this->flatLevels;
//    delete this->levelGridPoints;
//  }
//
//  if (this->allSubspaces != nullptr) {
//    delete this->allSubspaces;
//  }

    if (this->allSurplusses != nullptr) {
        delete[] this->allSurplusses;
        //delete[] this->levelGridPoints;
    }
    this->allSurplussesIndexMap.clear();

    if (this->allSubspaces != nullptr) {
        delete[] this->allSubspaces;
    }
}

void X86Baseline::prepare() {

//  //this requires old data, so free it before overriding everything else
//  if (this->flatLevels != nullptr) {
//    DataVector level(this->dim);
//    for (size_t i = 0; i < this->subspaceCount; i++) {
//      //allLevels.getRow(i, level);
//      int flatLevel = this->allSubspaces[i * this->subspaceSize + (2 * dim) + 1];
//      //int flatLevel = this->flattenLevel(this->dim, this->maxLevel, level);
//      delete this->flatLevels[flatLevel];
//    }
//    delete this->flatLevels;
//    delete this->levelGridPoints;
//  }

    if (this->allSurplusses != nullptr) {
        delete[] this->allSurplusses;
        //delete[] this->levelGridPoints;
    }
    this->allSurplussesIndexMap.clear();

    this->prepareSubspaceIterator();

    this->createFlatStorage();
}

/*
 levels needs to have size of the number of subspaces, stores a level tuple for each subspace
 level level part of a level index tuple, needs to have size of grid
 index index part of a level index tuple, needs to have size of grid

 Important: In contrast to other eval preparation method, this creates the actual level, not 2**level
 */
void X86Baseline::prepareSubspaceIterator() {

    /////////////////////////////////////////////////////
    // extract the subspace and grid points
    // and put them in a map with flatLevel -> subspace
    /////////////////////////////////////////////////////

    map<size_t, SubspaceNode> allLevelsMap;

    size_t maxLevel = 0;

    DataVector level(this->dim);
    DataVector index(this->dim);
    DataVector maxIndex(this->dim);

    unsigned int curLevel;
    unsigned int curIndex;

    //calculate the maxLevel first
    //TODO: this should be removed, maxLevel shouldn't be needed
    for (size_t gridIndex = 0; gridIndex < this->storage->size(); gridIndex++) {
        sg::base::GridIndex *point = this->storage->get(gridIndex);
        for (size_t d = 0; d < this->dim; d++) {
            point->get(d, curLevel, curIndex);
            if (curLevel > maxLevel) {
                maxLevel = curLevel;
            }
        }
    }

    //create map of subspaces -> now we know which subspaces actually exist in the grid
    for (size_t gridIndex = 0; gridIndex < this->storage->size(); gridIndex++) {
        sg::base::GridIndex *point = this->storage->get(gridIndex);

        for (size_t d = 0; d < this->dim; d++) {
            point->get(d, curLevel, curIndex);
            level.set(d, curLevel);
            index.set(d, curIndex);
            maxIndex.set(d, 1 << curLevel);
        }

        size_t flatLevel = X86Baseline::flattenLevel(this->dim, maxLevel, level);

        if (allLevelsMap.find(flatLevel) == allLevelsMap.end()) {
            SubspaceNode newNode(level, maxIndex, index);
            allLevelsMap.insert(std::pair<size_t, SubspaceNode>(flatLevel, newNode));
            //allLevelsMap.emplace(flatLevel, SubspaceNode(level, maxIndex, index));
        } else {
            SubspaceNode &subspace = allLevelsMap.find(flatLevel)->second;
            subspace.addGridPoint(index);
        }
    }

    size_t subspaceCount = allLevelsMap.size();
    this->subspaceCount = subspaceCount;
    //cout << "no. of subspaces: " << subspaceCount << endl;
    //cout << "maxLevel: " << maxLevel << endl;
    this->maxLevel = maxLevel;

    //////////////////////////////////////////////////////////////////
    // iterate though the subspaces for information (can be removed)
    //////////////////////////////////////////////////////////////////

    /*	double averageSubspaceUtilizationPercent = 0;
     int averagePointsOnSubspace = 0;

     for (typename std::map<size_t, SubspaceNode>::iterator it = allLevelsMap.begin();
     it != allLevelsMap.end(); ++it) {
     SubspaceNode &subspace = it->second;

     size_t gridPointsOnLevel = 1;
     for (size_t j = 0; j < this->dim; j++) {
     int dimTemp = subspace.hInverse[j];
     dimTemp >>= 1; //skip even indices
     gridPointsOnLevel *= dimTemp;
     }

     double subspaceUtilizationPercent = (subspace.actualGridPointsOnLevel / (double) gridPointsOnLevel) * 100.0;
     cout << level.toString() << " populated: " << subspace.actualGridPointsOnLevel;
     cout << " (max: " << gridPointsOnLevel << ")";
     cout << " -> " << subspaceUtilizationPercent << "%" << endl;
     averageSubspaceUtilizationPercent += subspaceUtilizationPercent;
     averagePointsOnSubspace += subspace.actualGridPointsOnLevel;
     }

     double overallUtilization = averageSubspaceUtilizationPercent / (double) this->allLevels->getNrows();
     cout << "averageSubspaceUtilizationPercent: " << overallUtilization << "%" << endl;
     double overallPointsOnSubspace = averagePointsOnSubspace / (double) this->allLevels->getNrows();
     cout << "averagePointsOnSubspace: " << overallPointsOnSubspace << endl;
     */

    ///////////////////////////////////////////////////////////////
    // now build the complex and fast iterable structure
    ///////////////////////////////////////////////////////////////
    // layout [#level #hInverse #next #flatLevel #nextDiff #jumpDiff #level2 ...]
    size_t subspaceSize = (2 * this->dim) + 4;
    this->subspaceSize = subspaceSize;

    // one special subspace at the end (basically for padding)
    //size_t *allSubspaces = (size_t *) malloc((subspaceCount + 1) * subspaceSize * sizeof(size_t));
    if (allSubspaces != nullptr) {
        delete allSubspaces;
    }
    size_t *allSubspaces = new size_t[(subspaceCount + 1) * subspaceSize];

    // create linear structure by iterate though the previously create array
    size_t i = 0;
    size_t linearLevelIndexCounter = 0;
    for (typename std::map<size_t, SubspaceNode>::iterator it = allLevelsMap.begin(); it != allLevelsMap.end(); ++it) {
        SubspaceNode &subspace = it->second;

        for (size_t j = 0; j < this->dim; j++) {
            level.set(j, subspace.level[j]);
            allSubspaces[i * subspaceSize + j] = subspace.level[j];
            allSubspaces[i * subspaceSize + this->dim + j] = subspace.hInverse[j];
        }
        //next dummy
        allSubspaces[i * subspaceSize + (2 * this->dim)] = 9999;
        //flatLevel
        size_t flatLevel = this->flattenLevel(this->dim, maxLevel, level);
        size_t linearLevelIndex = linearLevelIndexCounter;
        allSubspaces[i * subspaceSize + (2 * this->dim) + 1] = linearLevelIndex;
//        allSubspaces[i * subspaceSize + (2 * this->dim) + 1] = flatLevel;
        this->allSurplussesIndexMap[flatLevel] = linearLevelIndex;
        linearLevelIndexCounter += subspace.gridPointsOnLevel;

        //these diffs tell, which components (diff to <dim) of the index have to be recomputed
        //if a jump is performed or if no jump is performed
        //the diff information has to be available at the node before the jump
        //(at the jump destination it is unknown from which index the jump occured)
        //nextDiff dummy
        allSubspaces[i * subspaceSize + (2 * this->dim) + 2] = 9999;
        //jumpDiff dummy
        allSubspaces[i * subspaceSize + (2 * this->dim) + 3] = 9999;

        i += 1;
    }

    problemDimForSorting = this->dim;
    qsort(allSubspaces, subspaceCount, subspaceSize * sizeof(size_t), subspaceComparator);

    size_t computationFinishedMarker = subspaceCount * subspaceSize;

    allSubspaces[(subspaceCount - 1) * subspaceSize + (2 * this->dim)] = computationFinishedMarker;
    size_t *lastChangeIndex = new size_t[this->dim]; //(size_t *) malloc(sizeof(size_t) * this->dim);
    for (size_t i = 0; i < this->dim; i++) {
        lastChangeIndex[i] = computationFinishedMarker;
    }

    for (size_t i = subspaceCount - 1; i > 0; i--) {
        size_t currentItemIndex = ((size_t) i) * subspaceSize;
        size_t lastItemIndex = ((size_t) i - 1) * subspaceSize;
        size_t nextItemIndex = ((size_t) i + 1) * subspaceSize;
        size_t lastDiff = getHighestDifferentIndex(allSubspaces, currentItemIndex, lastItemIndex);

        // diff to the next element, if no jump is taken, required to get (near) O(1) index calculation
        size_t nextDiff = this->dim; //will result in doing nothing for the padding
        // diff to the jump destination (if the jump should be taken), required to get (near) O(1) index calculation
        size_t jumpDiff = this->dim;

        if (i != subspaceCount - 1) {
            nextDiff = getHighestDifferentIndex(allSubspaces, currentItemIndex, nextItemIndex);
            size_t jumpDestination = lastChangeIndex[lastDiff];
            if (jumpDestination != computationFinishedMarker) {
                jumpDiff = getHighestDifferentIndex(allSubspaces, currentItemIndex, jumpDestination);
            }
        }

        //next diff
        allSubspaces[currentItemIndex + (2 * this->dim) + 2] = nextDiff;

        //jump diff
        allSubspaces[currentItemIndex + (2 * this->dim) + 3] = jumpDiff;
        allSubspaces[currentItemIndex + (2 * this->dim)] = lastChangeIndex[lastDiff];

        for (size_t j = lastDiff + 1; j < dim; j++) {
            lastChangeIndex[j] = currentItemIndex;
        }

    }
    delete lastChangeIndex;

    // first subspace always leads to finished
    allSubspaces[2 * this->dim] = computationFinishedMarker;

    allSubspaces[2 * this->dim + 2] = 0;
    allSubspaces[2 * this->dim + 3] = 0;
    if (subspaceCount > 1) {
        allSubspaces[2 * this->dim + 2] = getHighestDifferentIndex(allSubspaces, 0, subspaceSize);
    }

    // setup padding subspace - last (pseudo) subspace
    size_t paddingIndex = subspaceCount * subspaceSize;
    for (size_t i = 0; i < this->dim; i++) { //level dummy = 1, ..., 1
        allSubspaces[paddingIndex + i] = 1;
    }
    for (size_t i = 0; i < this->dim; i++) { //hInverse dummy = 2, ..., 2
        allSubspaces[paddingIndex + this->dim + i] = 2;
    }
    allSubspaces[paddingIndex + 2 * this->dim] = computationFinishedMarker;
    allSubspaces[paddingIndex + 2 * this->dim + 1] = 0;
    //next diff
    allSubspaces[paddingIndex + (2 * this->dim) + 2] = this->dim;
    //jump diff
    allSubspaces[paddingIndex + (2 * this->dim) + 3] = this->dim;

    this->allSubspaces = allSubspaces;
}

void X86Baseline::createFlatStorage() {
    //this->levelGridPoints = new size_t[subspaceCount];

    //create data structure

    //TODO: calculation flawed! does only count sparse grid levels, index calculation uses full grid subspaces!
    //this->flatLevels = new double*[MAX_SUBSPACES];

    totalGridPoints = 0;
    for (size_t i = 0; i < this->subspaceCount; i++) {

        //DataVector level(dim);
        //allLevels.getRow(i, level);

        size_t *hInversePtr = this->allSubspaces + i * this->subspaceSize + dim;
        //DataVector hInverse(dim);
        //allHInverse.getRow(i, hInverse);
        //int flatLevel = this->flattenLevel(this->dim, this->maxLevel, level);
        //int flatLevel = this->allSubspaces[i * this->subspaceSize + (2 * dim) + 1];

        size_t gridPointsOnLevel = 1;
        for (size_t j = 0; j < this->dim; j++) {
            //cout << hInversePtr[j] << endl;
            //int dimTemp = hInverse[j];
            int dimTemp = hInversePtr[j];
            dimTemp >>= 1; //skip even indices
            gridPointsOnLevel *= dimTemp;
        }
        //this->levelGridPoints[i] = gridPointsOnLevel;

        totalGridPoints += gridPointsOnLevel;

        //cout << "creating fL: "  << flatLevel << " points: " << gridPointsOnLevel << endl;
//        double *levelArray = new double[gridPointsOnLevel];
//        this->flatLevels[flatLevel] = levelArray;
    }

    // the linearLevelIndex in allSubspaces encodes the start of the individual level array
    this->allSurplusses = new double[totalGridPoints];
}

void X86Baseline::setCoefficients(DataVector &surplusVector) {
    DataVector level(dim);
    DataVector maxIndex(dim);
    DataVector index(dim);

    for (size_t i = 0; i < this->totalGridPoints; i++) {
        this->allSurplusses[i] = std::numeric_limits<double>::quiet_NaN();
    }

    //cout << "surplus size: " << surplusVector.getSize() << endl;

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
        /*cout << "l: " << level.toString() << endl;
         cout << "i: " << index.toString() << endl;
         cout << "maxIndex: " << maxIndex.toString() << endl;*/

        this->setSurplus(level, maxIndex, index, surplusVector.get(gridIndex));
    }
}

//writes a result vector in the order of the points in the grid storage
void X86Baseline::unflatten(DataVector &result) {
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

        //cout << "result i: " << gridIndex << " => " << surplus << endl;
        result.set(gridIndex, surplus);
        //cout << "partialFin[" << gridIndex << "]: " << surplus << endl;
    }
    //cout << "partialFin: " << result.get(0) << endl;
}

void X86Baseline::setSurplus(DataVector &level, DataVector &maxIndices, DataVector &index, double value) {
    int levelFlat = this->flattenLevel(this->dim, this->maxLevel, level);
    size_t linearLevelIndex = this->allSurplussesIndexMap[levelFlat];
    int indexFlat = this->flattenIndex(this->dim, maxIndices, index);

//    double *levelArray = this->flatLevels[levelFlat];
    double *levelArray = &(this->allSurplusses[linearLevelIndex]);
    levelArray[indexFlat] = value;
}

void X86Baseline::getSurplus(DataVector &level, DataVector &maxIndices, DataVector &index, double &value,
        bool &isVirtual) {
    int levelFlat = this->flattenLevel(this->dim, this->maxLevel, level);
    size_t linearLevelIndex = this->allSurplussesIndexMap[levelFlat];
    int indexFlat = this->flattenIndex(this->dim, maxIndices, index);
//    double *levelArray = this->flatLevels[levelFlat];
    double *levelArray = &(this->allSurplusses[linearLevelIndex]);
    value = levelArray[indexFlat];
    if (std::isnan(value)) {
        isVirtual = true;
    } else {
        isVirtual = false;
    }
}

size_t X86Baseline::flattenIndex(size_t dim, DataVector &maxIndices, DataVector &index) {
    size_t indexFlat = index.get(0);
    indexFlat >>= 1;
    for (size_t i = 1; i < dim; i++) {
        int actualDirectionGridPoints = maxIndices.get(i);
        actualDirectionGridPoints >>= 1;
        indexFlat *= actualDirectionGridPoints;
        size_t actualIndex = index.get(i);
        actualIndex >>= 1; //divide index by 2, skip even indices
        indexFlat += actualIndex;
    }
    return indexFlat;
}

size_t X86Baseline::flattenLevel(size_t dim, size_t maxLevel, DataVector &level) {
    size_t levelFlat = 0;
    levelFlat += level.get(dim - 1);
    // loop terminates at -1
    for (int i = dim - 2; i >= 0; i--) {
        levelFlat *= maxLevel;
        levelFlat += level.get(i);
    }
    return levelFlat;
}

void X86Baseline::padDataset() {
    //size_t vecWidth = getVecWidth(vecType);
    size_t chunkSize = X86BASELINE_PARALLEL_DATA_POINTS;

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
}
