#include "../../OperationMultipleEvalSubspace/vectorized/X86Vectorized.hpp"
#include <limits>

#define MAX_SUBSPACES 100000000

X86Vectorized::X86Vectorized(GridStorage* storage, 
				     DataMatrix* dataset,
				     size_t dim) {      
  this->storage = storage;
  this->dataset = dataset;
  this->dim = dim;


  //this->subspaceSize = (2 * this->dim) + 4;
  this->subspaceSize = (2 * this->dim) + 6;
  this->levelOffset = 0;
  this->hInverseOffset = this->dim;
  this->nextOffset = 2 * this->dim;
  this->flatLevelOffset = (2 * this->dim) + 1;
  this->nextDiffOffset = (2 * this->dim) + 2;
  this->jumpDiffOffset = (2 * this->dim) + 3;
  this->flatLevelPointerOffset = (2 * this->dim) + 4; //double * => doubled weight in subspace size
}

X86Vectorized::~X86Vectorized() {
  // if (this->flatLevels != nullptr) {
  //   DataVector level(this->dim);
  //   for (size_t i = 0; i < this->subspaceCount; i++) {
  //     //allLevels.getRow(i, level);
  //     size_t flatLevel = this->allSubspaces[GET_FLATLEVEL(i)];
  //     //int flatLevel = this->flattenLevel(this->dim, this->maxLevel, level);
  //     delete this->flatLevels[flatLevel];
  //   }
  //   delete this->flatLevels;
  //   delete this->levelGridPoints;
  // }
  if (this->allSurplusses != nullptr) {
    delete this->allSurplusses;
    delete this->levelGridPoints;
  }

  if (this->allSubspaces != nullptr) {
    delete this->allSubspaces;
  }
}

void X86Vectorized::prepare() {

  //this requires old data, so free it before overriding everything else
  // if (this->flatLevels != nullptr) {
  //   DataVector level(this->dim);
  //   for (size_t i = 0; i < this->subspaceCount; i++) {
  //     //allLevels.getRow(i, level);
  //     size_t flatLevel = this->allSubspaces[GET_FLATLEVEL(i)];
  //     //int flatLevel = this->flattenLevel(this->dim, this->maxLevel, level);
  //     delete this->flatLevels[flatLevel];
  //   }
  //   delete this->flatLevels;
  //   delete this->levelGridPoints;
  // }

  if (this->allSurplusses != nullptr) {
    delete this->allSurplusses;
    delete this->levelGridPoints;
  }
  
  this->prepareSubspaceIterator();

  this->createFlatStorage();
}

void X86Vectorized::createFlatStorage() {
  this->levelGridPoints = new uint32_t[subspaceCount];

  //create data structure

  //TODO: calculation flawed! does only count sparse grid levels, index calculation uses full grid subspaces!
  //this->flatLevels = new double*[MAX_SUBSPACES];

  // for (size_t i = 0; i < this->subspaceCount; i++) {

  //   uint32_t *hInversePtr = GET_HINVERSE_PTR(i);
  //   int flatLevel = this->allSubspaces[GET_FLATLEVEL(i)];

  //   uint32_t gridPointsOnLevel = 1;
  //   for (size_t j = 0; j < this->dim; j++) {
  //     int dimTemp = hInversePtr[j];
  //     dimTemp >>= 1; //skip even indices
  //     gridPointsOnLevel *= dimTemp;
  //   }
  //   this->levelGridPoints[i] = gridPointsOnLevel;

  //   //cout << "creating fL: "  << flatLevel << " points: " << gridPointsOnLevel << endl;
  //   double *levelArray = new double[gridPointsOnLevel];
  //   this->flatLevels[flatLevel] = levelArray;
  // }

  this->totalRegularGridPoints = 0;
  for (size_t i = 0; i < this->subspaceCount; i++) {

    uint32_t *hInversePtr = GET_HINVERSE_PTR(i);

    uint32_t gridPointsOnLevel = 1;
    for (size_t j = 0; j < this->dim; j++) {
      int dimTemp = hInversePtr[j];
      dimTemp >>= 1; //skip even indices
      gridPointsOnLevel *= dimTemp;
    }
    this->totalRegularGridPoints += gridPointsOnLevel;
  }

  flatLevelArrayMap.clear();
  this->allSurplusses = new double[this->totalRegularGridPoints];
  cout << "totalRegularGridPoints: " << this->totalRegularGridPoints << endl;
  double *nextSubspaceStart = allSurplusses;
  for (size_t i = 0; i < this->subspaceCount; i++) {

    uint32_t *hInversePtr = GET_HINVERSE_PTR(i);
    int flatLevel = this->allSubspaces[GET_FLATLEVEL(i)];
    this->flatLevelArrayMap[flatLevel] = nextSubspaceStart;

    uint32_t gridPointsOnLevel = 1;
    for (size_t j = 0; j < this->dim; j++) {
      int dimTemp = hInversePtr[j];
      dimTemp >>= 1; //skip even indices
      gridPointsOnLevel *= dimTemp;
    }

    //set the pointer to the current subspace in the allSubspaces structure
    *GET_FLATLEVELPTR_ON_PTR(i) = nextSubspaceStart;
    // move the pointer to the start of the next subspace
    nextSubspaceStart += gridPointsOnLevel;
  }
  
}

void X86Vectorized::setCoefficients(DataVector &surplusVector) {
  DataVector level(dim);
  DataVector maxIndex(dim);
  DataVector index(dim);

  // //zero all coefficients (including virtual coefficients)
  // for (size_t i = 0; i < this->subspaceCount; i++) {
  //   int flatLevel = allSubspaces[GET_FLATLEVEL(i)];
  //   //initialize with NaN, so that virtual grid points can be recognized
  //   double *levelArray = this->flatLevels[flatLevel];
  //   for (size_t j = 0; j < this->levelGridPoints[i]; j++) {
  //     levelArray[j] = std::numeric_limits<double>::quiet_NaN();
  //   }

  // }

  //zero all coefficients (including virtual coefficients)
  for (size_t i = 0; i < this->totalRegularGridPoints; i++) {
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
void X86Vectorized::unflatten(DataVector &result) {
  DataVector level(dim);
  DataVector maxIndex(dim);
  DataVector index(dim);  

  unsigned int curLevel;
  unsigned int curIndex;
  for (size_t gridIndex = 0; gridIndex < this->storage->size(); gridIndex ++) {
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

void X86Vectorized::setSurplus(DataVector &level, DataVector &maxIndices, DataVector &index, double value) {
  int levelFlat = this->flattenLevel(this->dim, this->maxLevel, level);
  int indexFlat = this->flattenIndex(this->dim, maxIndices, index);
  // double *levelArray = this->flatLevels[levelFlat];
  // levelArray[indexFlat] = value;
  this->flatLevelArrayMap[levelFlat][indexFlat] = value;
}

void X86Vectorized::getSurplus(DataVector &level, DataVector &maxIndices, DataVector &index, double &value, bool &isVirtual) {
  int levelFlat = this->flattenLevel(this->dim, this->maxLevel, level);
  int indexFlat = this->flattenIndex(this->dim, maxIndices, index);  
  // double *levelArray = this->flatLevels[levelFlat];
  //value = levelArray[indexFlat];
  //double value2 = this->flatLevelArrayMap[levelFlat][indexFlat];
  value = this->flatLevelArrayMap[levelFlat][indexFlat];
  if (std::isnan(value)) {
    isVirtual = true;
  } else {
    isVirtual = false;
  }
}

uint32_t X86Vectorized::flattenIndex(size_t dim, DataVector &maxIndices, DataVector &index) {
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

uint32_t X86Vectorized::flattenLevel(size_t dim, size_t maxLevel, DataVector &level) {
  uint32_t levelFlat = 0;
  levelFlat += level.get(dim - 1);
  // loop terminates at -1
  for (int i = dim - 2; i >= 0; i--) {
    levelFlat *= maxLevel;
    levelFlat += level.get(i);
  }
  return levelFlat;
}  

void X86Vectorized::padDataset() {
  //size_t vecWidth = getVecWidth(vecType);
  size_t chunkSize = X86VECTORIZED_PARALLEL_DATA_POINTS;

  // Assure that data has a even number of instances -> padding might be needed
  size_t remainder = this->dataset->getNrows() % chunkSize;
  size_t loopCount = chunkSize - remainder;
  //cout << "loopCount: " << loopCount << endl;
  //cout << "vecOff: " << X86VECTORIZED_PARALLEL_DATA_POINTS << endl;

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
  this->dataset->addSize(X86VECTORIZED_VEC_PADDING * 2);
  for (size_t i = this->dataset->getNrows(); i < this->dataset->getNrows() + this->dataset->getUnused(); i++) {
    for (size_t j = 0; j < this->dataset->getNcols(); j++) {
      this->dataset->set(i, j, 0.0);
    }
  }
}
