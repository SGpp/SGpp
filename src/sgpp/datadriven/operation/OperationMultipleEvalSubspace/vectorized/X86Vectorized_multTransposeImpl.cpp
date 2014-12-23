#include <iomanip>
#include <algorithm>

//#include <immintrin.h>

#include "../../OperationMultipleEvalSubspace/vectorized/X86Vectorized.hpp"
#include "../../OperationMultipleEvalSubspace/vectorized/X86VectorizedTypes.hpp"

void X86Vectorized::multTransposeImpl(sg::base::DataVector &source,
				    sg::base::DataVector &result,
				    const size_t start_index_data,
				    const size_t end_index_data) {
  size_t tid = omp_get_thread_num();
  if (tid == 0) {
    this->setCoefficients(source);
  }

#pragma omp barrier

  size_t dim = dataset->getNcols();
  double * datasetPtr = dataset->getPointer();

  size_t totalThreadNumber = X86VECTORIZED_PARALLEL_DATA_POINTS + X86VECTORIZED_VEC_PADDING;

  double *evalIndexValuesAll = new double[(dim + 1) * totalThreadNumber];
  for (size_t i = 0; i < (dim + 1) * totalThreadNumber; i++) {
    evalIndexValuesAll[i] = 1.0;
  }

  //double **flatLevels = this->flatLevels;

  //for faster index flattening, last element is for padding
  uint32_t *intermediatesAll = new uint32_t[(dim + 1) * totalThreadNumber];
  for (size_t i = 0; i < (dim + 1) * totalThreadNumber; i++) {
    intermediatesAll[i] = 0.0;
  }

  size_t maxIndex = subspaceCount * this->subspaceSize;

  size_t validIndices[X86VECTORIZED_PARALLEL_DATA_POINTS + X86VECTORIZED_VEC_PADDING];
  size_t validIndicesCount;

  double componentResults[X86VECTORIZED_PARALLEL_DATA_POINTS + X86VECTORIZED_VEC_PADDING];	 
  size_t levelIndices[X86VECTORIZED_PARALLEL_DATA_POINTS + X86VECTORIZED_VEC_PADDING];
  size_t nextIterationToRecalcReferences[X86VECTORIZED_PARALLEL_DATA_POINTS + X86VECTORIZED_VEC_PADDING]; 

  // double chunkData[(X86VECTORIZED_PARALLEL_DATA_POINTS + X86VECTORIZED_VEC_PADDING) * dim];
	
  //process the next chunk of data tuples in parallel
  for (size_t dataIndexBase = start_index_data; dataIndexBase < end_index_data; 
       dataIndexBase += X86VECTORIZED_PARALLEL_DATA_POINTS) {

    for (size_t i = 0; i < totalThreadNumber; i++) {
      levelIndices[i] = 0.0;
      componentResults[i] = 0.0;
      nextIterationToRecalcReferences[i] = 0;
      // for (size_t j = 0; j < dim; j++) {
      // 	chunkData[i * dim + j] = 
      // 	  datasetPtr[(dataIndexBase + i) * dim + j];
      // }
    }

    for (size_t levelIndex = 0; levelIndex < maxIndex; levelIndex += this->subspaceSize) {

      uint32_t *hInversePtr = allSubspaces + GET_HINVERSE_LOOP(levelIndex);
      //uint32_t levelFlat = allSubspaces[GET_FLATLEVEL_LOOP(levelIndex)];
      //double *levelArray = flatLevels[levelFlat];
      double *levelArrayContinuous = *GET_FLATLEVELPTR_ON_PTR((levelIndex / this->subspaceSize));

      validIndicesCount = 0;
      for (size_t parallelIndex = 0; parallelIndex < X86VECTORIZED_PARALLEL_DATA_POINTS; 
	   parallelIndex++) {
	size_t parallelLevelIndex = levelIndices[parallelIndex];	
	if (parallelLevelIndex == levelIndex) {
	  validIndices[validIndicesCount] = parallelIndex;
	  validIndicesCount += 1;
	}
      }

      size_t paddingSize = min((int) (validIndicesCount + 4), X86VECTORIZED_PARALLEL_DATA_POINTS);
      for (size_t i = validIndicesCount; i < paddingSize; i++) {
	size_t threadId = X86VECTORIZED_PARALLEL_DATA_POINTS + (i - validIndicesCount);
	validIndices[i] = threadId;
	componentResults[threadId] = 0.0;
	levelIndices[threadId] = 0;
	nextIterationToRecalcReferences[threadId] = 0;
	double *evalIndexValues = evalIndexValuesAll + (dim + 1) * threadId;

	//for faster index flattening, last element is for padding
	uint32_t *intermediates = intermediatesAll + (dim + 1) * threadId;	  
	for (size_t j = 0; j < dim; j++) {
	  //indexPtr[j] = 1;
	  evalIndexValues[j] = 1.0;
	  intermediates[j] = 0.0;
	}
      }
            
      for (size_t validIndex = 0; validIndex < validIndicesCount; validIndex += 4) {

	size_t parallelIndices[4];
	parallelIndices[0] = validIndices[validIndex];
	parallelIndices[1] = validIndices[validIndex + 1];
	parallelIndices[2] = validIndices[validIndex + 2];
	parallelIndices[3] = validIndices[validIndex + 3];

	size_t nextIterationToRecalc = nextIterationToRecalcReferences[parallelIndices[0]];

	// double *dataTuplePtr[4];
	// dataTuplePtr[0] = datasetPtr + (dataIndexBase + parallelIndices[0]) * dim;
	// dataTuplePtr[1] = datasetPtr + (dataIndexBase + parallelIndices[1]) * dim;
	// dataTuplePtr[2] = datasetPtr + (dataIndexBase + parallelIndices[2]) * dim;
	// dataTuplePtr[3] = datasetPtr + (dataIndexBase + parallelIndices[3]) * dim;

	// double * dataTuplePtr[4];
	// dataTuplePtr[0] = chunkData + parallelIndices[0] * dim;
	// dataTuplePtr[1] = chunkData + parallelIndices[1] * dim;
	// dataTuplePtr[2] = chunkData + parallelIndices[2] * dim;
	// dataTuplePtr[3] = chunkData + parallelIndices[3] * dim;

	// const double * const dataTuplePtr[4] = 
	//   {chunkData + parallelIndices[0] * dim,
	//    chunkData + parallelIndices[1] * dim,
	//    chunkData + parallelIndices[2] * dim,
	//    chunkData + parallelIndices[3] * dim};


	const double * const dataTuplePtr[4] = 
	  {datasetPtr + (dataIndexBase + parallelIndices[0]) * dim,
	   datasetPtr + (dataIndexBase + parallelIndices[1]) * dim,
	   datasetPtr + (dataIndexBase + parallelIndices[2]) * dim,
	   datasetPtr + (dataIndexBase + parallelIndices[3]) * dim};

	double *evalIndexValues[4];
	evalIndexValues[0] = evalIndexValuesAll + (dim + 1) * parallelIndices[0];
	evalIndexValues[1] = evalIndexValuesAll + (dim + 1) * parallelIndices[1];
	evalIndexValues[2] = evalIndexValuesAll + (dim + 1) * parallelIndices[2];
	evalIndexValues[3] = evalIndexValuesAll + (dim + 1) * parallelIndices[3];

	//for faster index flattening, last element is for padding
	uint32_t *intermediates[4];
	intermediates[0] = intermediatesAll + (dim + 1) * parallelIndices[0];
	intermediates[1] = intermediatesAll + (dim + 1) * parallelIndices[1];
	intermediates[2] = intermediatesAll + (dim + 1) * parallelIndices[2];
	intermediates[3] = intermediatesAll + (dim + 1) * parallelIndices[3];

	uint32_t indexFlat[4];
	double phiEval[4];
	X86Vectorized::calculateIndexVectorized(dim, nextIterationToRecalc, 
						dataTuplePtr, hInversePtr, 
						intermediates, 
						evalIndexValues,
						indexFlat, phiEval
						);

	// X86Vectorized::flattenIndexVectorized(intermediates, dim, hInversePtr, indexPtr, 
	//   				      nextIterationToRecalc, indexFlat);  

	// X86Vectorized::evaluateBaseVectorized(dim, nextIterationToRecalc, dataTuplePtr, 
	// 				      evalIndexValues, hInversePtr, indexPtr, phiEval);

	/*double surplus[4];
	surplus[0] = levelArray[indexFlat[0]];
	surplus[1] = levelArray[indexFlat[1]];
	surplus[2] = levelArray[indexFlat[2]];
	surplus[3] = levelArray[indexFlat[3]];*/

	double surplus[4];
	surplus[0] = levelArrayContinuous[indexFlat[0]];
	surplus[1] = levelArrayContinuous[indexFlat[1]];
	surplus[2] = levelArrayContinuous[indexFlat[2]];
	surplus[3] = levelArrayContinuous[indexFlat[3]];

	//TODO: either transform to speedup or make loop

	for (size_t innerIndex = 0; innerIndex < 4; innerIndex++) {
	  if (!std::isnan(surplus[innerIndex])) {
	    componentResults[parallelIndices[innerIndex]] += phiEval[innerIndex] * surplus[innerIndex];
	    nextIterationToRecalcReferences[parallelIndices[innerIndex]] = allSubspaces[levelIndex + (2 * dim) + 2];
	    levelIndices[parallelIndices[innerIndex]] += this->subspaceSize;
	  } else {
#if X86VECTORIZED_ENABLE_SUBSPACE_SKIPPING == 1	      
	    //skip to next relevant subspace
	    nextIterationToRecalcReferences[parallelIndices[innerIndex]]  = 
	      allSubspaces[levelIndex + (2 * dim) + 3];
	    levelIndices[parallelIndices[innerIndex]] = allSubspaces[levelIndex + 2 * dim];
#else
	    nextIterationToRecalcReferences[parallelIndices[innerIndex]] = 
	      allSubspaces[levelIndex + (2 * dim) + 2];
	    levelIndices[parallelIndices[innerIndex]] += this->subspaceSize;
#endif
	  }
	}

// 	if (!std::isnan(surplus[0])) {
// 	  componentResults[parallelIndices[0]] += phiEval[0] * surplus[0];
// 	  nextIterationToRecalcReferences[parallelIndices[0]] = allSubspaces[levelIndex + (2 * dim) + 2];
// 	  levelIndices[parallelIndices[0]] += this->subspaceSize;
// 	}
// 	if (!std::isnan(surplus[1])) {
// 	  componentResults[parallelIndices[1]] += phiEval[1] * surplus[1];
// 	  nextIterationToRecalcReferences[parallelIndices[1]] = allSubspaces[levelIndex + (2 * dim) + 2];
// 	  levelIndices[parallelIndices[1]] += this->subspaceSize;
// 	}
// 	if (!std::isnan(surplus[2])) {
// 	  componentResults[parallelIndices[2]] += phiEval[2] * surplus[2];
// 	  nextIterationToRecalcReferences[parallelIndices[2]] = allSubspaces[levelIndex + (2 * dim) + 2];
// 	  levelIndices[parallelIndices[2]] += this->subspaceSize;
// 	}
// 	if (!std::isnan(surplus[3])) {
// 	  componentResults[parallelIndices[3]] += phiEval[3] * surplus[3];
// 	  nextIterationToRecalcReferences[parallelIndices[3]] = allSubspaces[levelIndex + (2 * dim) + 2];
// 	  levelIndices[parallelIndices[3]] += this->subspaceSize;
// 	}

// 	if (std::isnan(surplus[0])) {
// #if X86VECTORIZED_ENABLE_SUBSPACE_SKIPPING == 1	      
// 	  //skip to next relevant subspace
// 	  nextIterationToRecalcReferences[parallelIndices[0]]  = 
// 	    allSubspaces[levelIndex + (2 * dim) + 3];
// 	  levelIndices[parallelIndices[0]] = allSubspaces[levelIndex + 2 * dim];
// #else
// 	  nextIterationToRecalcReferences[parallelIndices[0]] = 
// 	    allSubspaces[levelIndex + (2 * dim) + 2];
// 	  levelIndices[parallelIndices[0]] += this->subspaceSize;
// #endif
// 	}
// 	if (std::isnan(surplus[1])) {
// #if X86VECTORIZED_ENABLE_SUBSPACE_SKIPPING == 1	      
// 	  //skip to next relevant subspace
// 	  nextIterationToRecalcReferences[parallelIndices[1]]  = 
// 	    allSubspaces[levelIndex + (2 * dim) + 3];
// 	  levelIndices[parallelIndices[1]] = allSubspaces[levelIndex + 2 * dim];
// #else
// 	  nextIterationToRecalcReferences[parallelIndices[1]] = 
// 	    allSubspaces[levelIndex + (2 * dim) + 2];
// 	  levelIndices[parallelIndices[1]] += this->subspaceSize;
// #endif
// 	}
// 	if (std::isnan(surplus[2])) {
// #if X86VECTORIZED_ENABLE_SUBSPACE_SKIPPING == 1	      
// 	  //skip to next relevant subspace
// 	  nextIterationToRecalcReferences[parallelIndices[2]]  = 
// 	    allSubspaces[levelIndex + (2 * dim) + 3];
// 	  levelIndices[parallelIndices[2]] = allSubspaces[levelIndex + 2 * dim];
// #else
// 	  nextIterationToRecalcReferences[parallelIndices[2]] = 
// 	    allSubspaces[levelIndex + (2 * dim) + 2];
// 	  levelIndices[parallelIndices[2]] += this->subspaceSize;
// #endif
// 	}
// 	if (std::isnan(surplus[3])) {
// #if X86VECTORIZED_ENABLE_SUBSPACE_SKIPPING == 1	      
// 	  //skip to next relevant subspace
// 	  nextIterationToRecalcReferences[parallelIndices[3]]  = 
// 	    allSubspaces[levelIndex + (2 * dim) + 3];
// 	  levelIndices[parallelIndices[3]] = allSubspaces[levelIndex + 2 * dim];
// #else
// 	  nextIterationToRecalcReferences[parallelIndices[3]] = 
// 	    allSubspaces[levelIndex + (2 * dim) + 2];
// 	  levelIndices[parallelIndices[3]] += this->subspaceSize;
// #endif
// 	}

      } //end X86VECTORIZED_PARALLEL_DATA_POINTS
    } // end iterate grid
	  
    for (size_t parallelIndex = 0; parallelIndex < X86VECTORIZED_PARALLEL_DATA_POINTS; parallelIndex++) {
      size_t dataIndex = dataIndexBase + parallelIndex;
      result.set(dataIndex, componentResults[parallelIndex]);	    
    }
  } // end iterate data chunks
  
  //delete indexPtrAll;
  delete evalIndexValuesAll;
  delete intermediatesAll;
}
