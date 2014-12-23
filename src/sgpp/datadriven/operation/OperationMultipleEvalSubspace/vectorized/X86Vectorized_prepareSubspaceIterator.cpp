#include "../../OperationMultipleEvalSubspace/vectorized/X86Vectorized.hpp"

namespace sg {
  namespace parallel {
    namespace x86vectorized {

      size_t problemDimForSorting;

      int subspaceComparator(const void *first, const void *second) {
	for (size_t i = 0; i < problemDimForSorting; i++) {
	  if (((uint32_t *) first)[i] >= ((uint32_t *) second)[i]) {
	    if (((uint32_t *) first)[i] > ((uint32_t *) second)[i]) {
	      return 1;
	    }
	  } else {
	    return 0;
	  }
	}
	return 0;
      }

      uint32_t getHighestDifferentIndex(uint32_t *base, size_t cur, size_t last) {
	//for (int i = problemDimForSorting - 1; i >= 0; i--) {
	for (size_t i = 0; i < problemDimForSorting; i++) {
	  if (base[cur + (size_t) i] != base[last + (size_t) i]) {
	    return i;
	  }
	}
	return -1; //should never be returned
      }

      void printSubspace(uint32_t *base, size_t dim, size_t cur, int subspaceSize) {
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

using namespace sg::parallel::x86vectorized;

/*
  levels needs to have size of the number of subspaces, stores a level tuple for each subspace 
  level level part of a level index tuple, needs to have size of grid
  index index part of a level index tuple, needs to have size of grid

  Important: In contrast to other eval preparation method, this creates the actual level, not 2**level
*/
void X86Vectorized::prepareSubspaceIterator() {

  /////////////////////////////////////////////////////
  // extract the subspace and grid points
  // and put them in a map with flatLevel -> subspace
  /////////////////////////////////////////////////////
	
  map<uint32_t, X86VectorizedSubspaceNode> allLevelsMap;

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

    uint32_t flatLevel = X86Vectorized::flattenLevel(this->dim, maxLevel, level);

    if (allLevelsMap.find(flatLevel) == allLevelsMap.end()) {
      X86VectorizedSubspaceNode newNode(level, maxIndex, index);
      allLevelsMap.insert(std::pair<size_t, X86VectorizedSubspaceNode>(flatLevel, newNode));
      //allLevelsMap.emplace(flatLevel, X86VectorizedSubspaceNode(level, maxIndex, index));
    } else {
      X86VectorizedSubspaceNode &subspace = allLevelsMap.find(flatLevel)->second;
      subspace.addGridPoint(index);
    }
  }

  this->subspaceCount = allLevelsMap.size();
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
	
  // one special subspace at the end (basically for padding)
  //size_t *allSubspaces = (size_t *) malloc((subspaceCount + 1) * subspaceSize * sizeof(size_t));
  if (allSubspaces != nullptr) {
    delete allSubspaces;
  }
  uint32_t *allSubspaces = new uint32_t[(this->subspaceCount + 1) * this->subspaceSize];

  // create linear structure by iterate though the previously create array
  size_t i = 0;
  for (typename std::map<uint32_t, X86VectorizedSubspaceNode>::iterator it = allLevelsMap.begin(); 
       it != allLevelsMap.end(); ++it) {	  
    X86VectorizedSubspaceNode &subspace = it->second;

    for (size_t j = 0; j < this->dim; j++) {
      level.set(j, subspace.level[j]);
      allSubspaces[GET_SUBS_ARRAY(i, this->levelOffset, j)] = subspace.level[j];
      allSubspaces[GET_SUBS_ARRAY(i, this->hInverseOffset, j)] = subspace.hInverse[j];
    }
    //next dummy
    allSubspaces[GET_NEXT(i)] = 9999;
    //flatLevel
    uint32_t flatLevel = this->flattenLevel(this->dim, maxLevel, level);
    allSubspaces[GET_FLATLEVEL(i)] = flatLevel;

    //these diffs tell, which components (diff to <dim) of the index have to be recomputed
    //if a jump is performed or if no jump is performed
    //the diff information has to be available at the node before the jump
    //(at the jump destination it is unknown from which index the jump occured)
    //nextDiff dummy 
    allSubspaces[GET_NEXTDIFF(i)] = 9999;
    //jumpDiff dummy
    allSubspaces[GET_JUMPDIFF(i)] = 9999;

    i += 1;
  }

  problemDimForSorting = this->dim;
  qsort(allSubspaces, subspaceCount, subspaceSize * sizeof(uint32_t), subspaceComparator);

  uint32_t computationFinishedMarker = subspaceCount * subspaceSize;
	  
  //allSubspaces[(subspaceCount - 1) * subspaceSize + (2 * this->dim)] = computationFinishedMarker;
  allSubspaces[GET_NEXT(subspaceSize - 1)] = computationFinishedMarker;
  uint32_t *lastChangeIndex = new uint32_t[this->dim];
  for (size_t i = 0; i < this->dim; i++) {
    lastChangeIndex[i] = computationFinishedMarker;
  }

  for (size_t i = subspaceCount - 1; i > 0; i--) {
    size_t currentItemIndex = ((size_t) i) * subspaceSize;
    size_t lastItemIndex = ((size_t) i - 1) * subspaceSize;
    size_t nextItemIndex = ((size_t) i + 1) * subspaceSize;
    uint32_t lastDiff = getHighestDifferentIndex(allSubspaces, currentItemIndex, lastItemIndex);

    // diff to the next element, if no jump is taken, required to get (near) O(1) index calculation
    uint32_t nextDiff = this->dim; //will result in doing nothing for the padding
    // diff to the jump destination (if the jump should be taken), required to get (near) O(1) index calculation
    uint32_t jumpDiff = this->dim;

    if (i != subspaceCount - 1) {
      nextDiff  = getHighestDifferentIndex(allSubspaces, currentItemIndex, nextItemIndex);
      uint32_t jumpDestination = lastChangeIndex[lastDiff];
      if (jumpDestination != computationFinishedMarker) {
	jumpDiff  = getHighestDifferentIndex(allSubspaces, currentItemIndex, jumpDestination);
      }
    }

    //next diff
    allSubspaces[GET_NEXTDIFF_LOOP(currentItemIndex)] = nextDiff;
	  
    //jump diff
    allSubspaces[GET_JUMPDIFF_LOOP(currentItemIndex)] = jumpDiff;    
    allSubspaces[GET_NEXT_LOOP(currentItemIndex)] = lastChangeIndex[lastDiff];

    for (size_t j = lastDiff + 1; j < dim; j++) {
      lastChangeIndex[j] = currentItemIndex;
    }

  }
  delete lastChangeIndex;

  // first subspace always leads to finished
  allSubspaces[GET_NEXT(0)] = computationFinishedMarker;

  allSubspaces[GET_NEXTDIFF(0)] = 0;
  allSubspaces[GET_JUMPDIFF(0)] = 0;
  if (subspaceCount > 1) {
    allSubspaces[GET_NEXTDIFF(0)] = getHighestDifferentIndex(allSubspaces, 0, subspaceSize);
  }

  // setup padding subspace - last (pseudo) subspace
  size_t paddingIndex = this->subspaceCount;
  for (size_t i = 0; i < this->dim; i++) { //level dummy = 1, ..., 1
    allSubspaces[GET_SUBS_ARRAY(paddingIndex, this->levelOffset, i)] = 1; 
  }
  for (size_t i = 0; i < this->dim; i++) { //hInverse dummy = 2, ..., 2
    allSubspaces[GET_SUBS_ARRAY(paddingIndex, this->hInverseOffset, i)] = 2;
  }
  allSubspaces[GET_NEXT(paddingIndex)] = computationFinishedMarker;
  allSubspaces[GET_FLATLEVEL(paddingIndex)] = 0;
  //next diff
  allSubspaces[GET_NEXTDIFF(paddingIndex)] = this->dim;	  
  //jump diff
  allSubspaces[GET_JUMPDIFF(paddingIndex)] = this->dim;

  this->allSubspaces = allSubspaces;  
}
