#include "../../OperationMultipleEvalSubspace/baseline/X86Baseline.hpp"

void X86Baseline::multTransposeImpl(sg::base::DataVector &source,
				    sg::base::DataVector &result,
				    //as only the non-zero grid points are evaluated,
				    //blocking over grid make no sense
				    //const size_t start_index_grid,
				    //const size_t end_index_grid,
				    const size_t start_index_data,
				    const size_t end_index_data) {

  //cout << "======================multTransposeImpl========================" << endl;

  size_t tid = omp_get_thread_num();
  if (tid == 0) {
    this->setCoefficients(source);
  }

#pragma omp barrier

  size_t dim = dataset->getNcols();
  double *datasetPtr = dataset->getPointer();

  //DataVector dataTuple(dim);
  //double *dataTuplePtr = dataTuple.getPointer();
  //size_t *indexPtr = (size_t *) malloc(sizeof(size_t) * dim);
  //indexPtr[0] = 0;
  size_t *indexPtrAll = new size_t[(dim+1) * X86BASELINE_PARALLEL_DATA_POINTS];
  for (size_t i = 0; i < (dim + 1) * X86BASELINE_PARALLEL_DATA_POINTS; i++) {
    indexPtrAll[i] = 0;
  }
  //double *evalIndexValues = (double *) malloc(sizeof(double) * (dim + 1));
  double *evalIndexValuesAll = new double[(dim + 1) * X86BASELINE_PARALLEL_DATA_POINTS];
  //evalIndexValues[0] = 1.0;
  for (size_t i = 0; i < (dim + 1) * X86BASELINE_PARALLEL_DATA_POINTS; i++) {
    evalIndexValuesAll[i] = 1.0;
  }

  //double **flatLevels = this->flatLevels;

  //for faster index flattening, last element is for padding
  //double *intermediates = (double *) malloc(sizeof(double) * (dim + 1));
  double *intermediatesAll = new double[(dim + 1) * X86BASELINE_PARALLEL_DATA_POINTS];
  //intermediates[0] = 0.0;
  for (size_t i = 0; i < (dim + 1) * X86BASELINE_PARALLEL_DATA_POINTS; i++) {
    intermediatesAll[i] = 0.0;
  }

  double maxIndex = subspaceCount * subspaceSize;

  //TODO: padding?
	
  //process the next chunk of data tuples in parallel
  //for (size_t dataIndex = start_index_data; dataIndex < end_index_data; dataIndex++) {
  for (size_t dataIndexBase = start_index_data; dataIndexBase < end_index_data; 
       dataIndexBase += X86BASELINE_PARALLEL_DATA_POINTS) {
    //cout << "chunk start: " << dataIndexBase << endl;

    //double componentResult = 0.0;
    double componentResults[X86BASELINE_PARALLEL_DATA_POINTS];	 
    //size_t levelIndex = 0;
    size_t levelIndices[X86BASELINE_PARALLEL_DATA_POINTS];
    //size_t nextIterationToRecalc = 0; //all index components are recalculated
    size_t nextIterationToRecalcReferences[X86BASELINE_PARALLEL_DATA_POINTS]; //all index components are recalculated
    for (size_t i = 0; i < X86BASELINE_PARALLEL_DATA_POINTS; i++) {
      levelIndices[i] = 0.0;
      componentResults[i] = 0.0;
      nextIterationToRecalcReferences[i] = 0;
    }

    //while (levelIndex < maxIndex) {
    for (size_t allLevelIndex = 0; allLevelIndex < maxIndex; allLevelIndex += subspaceSize) {
      //cout << "subspace index: " << (allLevelIndex / subspaceSize) << endl;

      size_t *hInversePtr = allSubspaces + allLevelIndex + dim;
      size_t linearLevelIndex = *(allSubspaces + allLevelIndex + (2 * dim) + 1);
      double *levelArray = &(this->allSurplusses[linearLevelIndex]);

      for (size_t parallelIndex = 0; parallelIndex < X86BASELINE_PARALLEL_DATA_POINTS; parallelIndex++) {
	//cout << "data points index: " << parallelIndex << endl;
	size_t levelIndex = levelIndices[parallelIndex];
	if (levelIndex != allLevelIndex) {
	  continue; // current element skips subspace
	}
	//size_t dataIndex = dataIndexBase + parallelIndex;

	//size_t *levelPtr = allSubspaces + levelIndex;
	size_t nextIterationToRecalc = nextIterationToRecalcReferences[parallelIndex];

	//double *dataTuplePtr = dataPoints + parallelIndex * dim;
	double *dataTuplePtr = datasetPtr + (dataIndexBase + parallelIndex) * dim;

	size_t *indexPtr = indexPtrAll + (dim + 1) * parallelIndex;
	double *evalIndexValues = evalIndexValuesAll + (dim + 1) * parallelIndex;

	//for faster index flattening, last element is for padding
	double *intermediates = intermediatesAll + (dim + 1) * parallelIndex;

	//TODO: calculate local data points at the beginning

	//TODO: could cache unadjusted for phi1D calculation
	for (size_t i = nextIterationToRecalc; i < dim; i++) {
	  //double unadjusted = dataTuplePtr[i] * hInversePtr[i];
	  double unadjusted = dataTuplePtr[i] * hInversePtr[i];
	  //cout << "dataPoints[" << (i) << "] = " << dataTuplePtr[i] << endl;
	  indexPtr[i] = calculateIndexComponent(dim, unadjusted);
	}

	/*cout << "l: ";
	  for (size_t i = 0; i < dim; i++) {
	  if (i > 0) {
	  cout << ", ";
	  }
	  cout << levelPtr[i];
	  }
	  cout << " ";

	  cout << "i: ";
	  for (size_t i = 0; i < dim; i++) {
	  if (i > 0) {
	  cout << ", ";
	  }
	  cout << indexPtr[i];
	  }
	  cout << endl;*/
	    
	int indexFlat = X86Baseline::flattenIndex(intermediates, dim, hInversePtr, indexPtr, nextIterationToRecalc);  
	//double *levelArray = flatLevels[levelFlat];
	double surplus = levelArray[indexFlat];

	if (!std::isnan(surplus)) {

	  //double phiEval = 1.0;
	  double phiEval = evalIndexValues[nextIterationToRecalc];
	  //prepare the values for the individual components
	  for (size_t i = nextIterationToRecalc; i < dim; i++) {
	    double phi1DEval = hInversePtr[i] * dataTuplePtr[i] - indexPtr[i];
	    phi1DEval = max(0.0, 1.0 - abs(phi1DEval));
	    phiEval *= phi1DEval;
	    evalIndexValues[i + 1] = phiEval;
	  }

	  //componentResult += phiEval * surplus;
	  componentResults[parallelIndex] += phiEval * surplus;
	  //nextIterationToRecalc = allSubspaces[levelIndex + (2 * dim) + 2];
	  nextIterationToRecalcReferences[parallelIndex] = allSubspaces[levelIndex + (2 * dim) + 2];
	  //nextIterationToRecalcReferences[parallelIndex] = 0;

	  //levelIndex += subspaceSize;
	  levelIndices[parallelIndex] += subspaceSize;

	} else {

#if X86BASELINE_ENABLE_SUBSPACE_SKIPPING == 1
	  //skip to next relevant subspace
	  //nextIterationToRecalc = allSubspaces[levelIndex + (2 * dim) + 3];
	  nextIterationToRecalcReferences[parallelIndex]  = allSubspaces[levelIndex + (2 * dim) + 3];
	  //nextIterationToRecalc = 0;
	  //nextIterationToRecalcReferences[parallelIndex]  = 0;
	  //levelIndex = allSubspaces[levelIndex + 2 * dim];
	  levelIndices[parallelIndex] = allSubspaces[levelIndex + 2 * dim];
	  //levelIndices[parallelIndex] += subspaceSize;

#else
	  nextIterationToRecalcReferences[parallelIndex] = allSubspaces[levelIndex + (2 * dim) + 2];
	  levelIndices[parallelIndex] += subspaceSize;
#endif

	}

      } //end X86BASELINE_PARALLEL_DATA_POINTS

    } // end iterate grid


    //result.set(dataIndex, componentResult);	    
    for (size_t parallelIndex = 0; parallelIndex < X86BASELINE_PARALLEL_DATA_POINTS; parallelIndex++) {
      size_t dataIndex = dataIndexBase + parallelIndex;
      result.set(dataIndex, componentResults[parallelIndex]);	    
    }
  } // end iterate data chunks
  
  delete indexPtrAll;
  delete evalIndexValuesAll;
  delete intermediatesAll;

}
