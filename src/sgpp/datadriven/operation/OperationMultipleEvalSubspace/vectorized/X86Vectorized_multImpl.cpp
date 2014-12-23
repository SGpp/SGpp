#include "../../OperationMultipleEvalSubspace/vectorized/X86Vectorized.hpp"

void X86Vectorized::multImpl(sg::base::DataVector &alpha,
				 sg::base::DataVector &result,
				 const size_t start_index_data,
				 const size_t end_index_data) {	  

  size_t tid = omp_get_thread_num();
  if (tid == 0) {
    this->setCoefficients(result);
  }

#pragma omp barrier

  size_t dim = this->dataset->getNcols();
  const double * const datasetPtr = this->dataset->getPointer();

  size_t totalThreadNumber = X86VECTORIZED_PARALLEL_DATA_POINTS + X86VECTORIZED_VEC_PADDING;

  double *evalIndexValuesAll = new double[(dim + 1) * totalThreadNumber];
  for (size_t i = 0; i < (dim + 1) * totalThreadNumber; i++) {
    evalIndexValuesAll[i] = 1.0;
  }

  //for faster index flattening
  uint32_t *intermediatesAll = new uint32_t[(dim + 1) * totalThreadNumber];
  for (size_t i = 0; i < (dim + 1) * totalThreadNumber; i++) {
    intermediatesAll[i] = 0.0;
  }

  size_t validIndices[X86VECTORIZED_PARALLEL_DATA_POINTS + X86VECTORIZED_VEC_PADDING];
  size_t validIndicesCount;

  size_t levelIndices[X86VECTORIZED_PARALLEL_DATA_POINTS + X86VECTORIZED_VEC_PADDING];
  size_t nextIterationToRecalcReferences[X86VECTORIZED_PARALLEL_DATA_POINTS + X86VECTORIZED_VEC_PADDING]; 

  // double chunkData[(X86VECTORIZED_PARALLEL_DATA_POINTS + X86VECTORIZED_VEC_PADDING) * dim];

  size_t maxIndex = this->subspaceCount * this->subspaceSize;
  for (size_t dataIndexBase = start_index_data; dataIndexBase < end_index_data; 
       dataIndexBase += X86VECTORIZED_PARALLEL_DATA_POINTS) {

    for (size_t i = 0; i < totalThreadNumber; i++) {
      levelIndices[i] = 0.0;
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
      double *levelArrayContinuous = *GET_FLATLEVELPTR_ON_PTR(levelIndex / this->subspaceSize);
      
      validIndicesCount = 0;
      for (size_t parallelIndex = 0; parallelIndex < X86VECTORIZED_PARALLEL_DATA_POINTS; 
	   parallelIndex++) {
	size_t parallelLevelIndex = levelIndices[parallelIndex];	
	if (parallelLevelIndex == levelIndex) {
	  validIndices[validIndicesCount] = parallelIndex;
	  validIndicesCount += 1;
	}
      }

      //padding for up to vector size, no padding required if all data tuples participate as
      //the number of data points is a multiple of the vector width
      size_t paddingSize = min((int) (validIndicesCount + X86VECTORIZED_VEC_PADDING), X86VECTORIZED_PARALLEL_DATA_POINTS);
      for (size_t i = validIndicesCount; i < paddingSize; i++) {
	size_t threadId = X86VECTORIZED_PARALLEL_DATA_POINTS + (i - validIndicesCount);
	validIndices[i] = threadId;
	levelIndices[threadId] = 0;
	nextIterationToRecalcReferences[threadId] = 0;
	double *evalIndexValues = evalIndexValuesAll + (dim + 1) * threadId;

	//for faster index flattening, last element is for padding
	uint32_t *intermediates = intermediatesAll + (dim + 1) * threadId;	  
	for (size_t j = 0; j < dim; j++) {
	  evalIndexValues[j] = 1.0;
	  intermediates[j] = 0;
	}
      }

      for (size_t validIndex = 0; validIndex < validIndicesCount; validIndex += X86VECTORIZED_VEC_PADDING) {

	size_t parallelIndices[4];
	parallelIndices[0] = validIndices[validIndex];
	parallelIndices[1] = validIndices[validIndex + 1];
	parallelIndices[2] = validIndices[validIndex + 2];
	parallelIndices[3] = validIndices[validIndex + 3];

	//pipe
	// size_t parallelIndices2[4];
	// parallelIndices2[0] = validIndices[validIndex + 4];
	// parallelIndices2[1] = validIndices[validIndex + 5];
	// parallelIndices2[2] = validIndices[validIndex + 6];
	// parallelIndices2[3] = validIndices[validIndex + 7];

	// size_t parallelIndices3[4];
	// parallelIndices3[0] = validIndices[validIndex + 8];
	// parallelIndices3[1] = validIndices[validIndex + 9];
	// parallelIndices3[2] = validIndices[validIndex + 10];
	// parallelIndices3[3] = validIndices[validIndex + 11];

	// size_t parallelIndices4[4];
	// parallelIndices4[0] = validIndices[validIndex + 12];
	// parallelIndices4[1] = validIndices[validIndex + 13];
	// parallelIndices4[2] = validIndices[validIndex + 14];
	// parallelIndices4[3] = validIndices[validIndex + 15];

	// size_t parallelIndices5[4];
	// parallelIndices5[0] = validIndices[validIndex + 16];
	// parallelIndices5[1] = validIndices[validIndex + 17];
	// parallelIndices5[2] = validIndices[validIndex + 18];
	// parallelIndices5[3] = validIndices[validIndex + 19];

	// size_t parallelIndices6[4];
	// parallelIndices6[0] = validIndices[validIndex + 20];
	// parallelIndices6[1] = validIndices[validIndex + 21];
	// parallelIndices6[2] = validIndices[validIndex + 22];
	// parallelIndices6[3] = validIndices[validIndex + 23];

	size_t nextIterationToRecalc = nextIterationToRecalcReferences[parallelIndices[0]];

	const double * const dataTuplePtr[4] = 
	  {datasetPtr + (dataIndexBase + parallelIndices[0]) * dim,
	   datasetPtr + (dataIndexBase + parallelIndices[1]) * dim,
	   datasetPtr + (dataIndexBase + parallelIndices[2]) * dim,
	   datasetPtr + (dataIndexBase + parallelIndices[3]) * dim};

	//pipe
	// const double * const dataTuplePtr2[4] = 
	//   {datasetPtr + (dataIndexBase + parallelIndices2[0]) * dim,
	//    datasetPtr + (dataIndexBase + parallelIndices2[1]) * dim,
	//    datasetPtr + (dataIndexBase + parallelIndices2[2]) * dim,
	//    datasetPtr + (dataIndexBase + parallelIndices2[3]) * dim};

	// const double * const dataTuplePtr3[4] = 
	//   {datasetPtr + (dataIndexBase + parallelIndices3[0]) * dim,
	//    datasetPtr + (dataIndexBase + parallelIndices3[1]) * dim,
	//    datasetPtr + (dataIndexBase + parallelIndices3[2]) * dim,
	//    datasetPtr + (dataIndexBase + parallelIndices3[3]) * dim};

	// const double * const dataTuplePtr4[4] = 
	//   {datasetPtr + (dataIndexBase + parallelIndices4[0]) * dim,
	//    datasetPtr + (dataIndexBase + parallelIndices4[1]) * dim,
	//    datasetPtr + (dataIndexBase + parallelIndices4[2]) * dim,
	//    datasetPtr + (dataIndexBase + parallelIndices4[3]) * dim};

	// const double * const dataTuplePtr5[4] = 
	//   {datasetPtr + (dataIndexBase + parallelIndices5[0]) * dim,
	//    datasetPtr + (dataIndexBase + parallelIndices5[1]) * dim,
	//    datasetPtr + (dataIndexBase + parallelIndices5[2]) * dim,
	//    datasetPtr + (dataIndexBase + parallelIndices5[3]) * dim};

	// const double * const dataTuplePtr6[4] = 
	//   {datasetPtr + (dataIndexBase + parallelIndices6[0]) * dim,
	//    datasetPtr + (dataIndexBase + parallelIndices6[1]) * dim,
	//    datasetPtr + (dataIndexBase + parallelIndices6[2]) * dim,
	//    datasetPtr + (dataIndexBase + parallelIndices6[3]) * dim};


	double *evalIndexValues[4];
	evalIndexValues[0] = evalIndexValuesAll + (dim + 1) * parallelIndices[0];
	evalIndexValues[1] = evalIndexValuesAll + (dim + 1) * parallelIndices[1];
	evalIndexValues[2] = evalIndexValuesAll + (dim + 1) * parallelIndices[2];
	evalIndexValues[3] = evalIndexValuesAll + (dim + 1) * parallelIndices[3];

	//pipe
	// double *evalIndexValues2[4];
	// evalIndexValues2[0] = evalIndexValuesAll + (dim + 1) * parallelIndices2[0];
	// evalIndexValues2[1] = evalIndexValuesAll + (dim + 1) * parallelIndices2[1];
	// evalIndexValues2[2] = evalIndexValuesAll + (dim + 1) * parallelIndices2[2];
	// evalIndexValues2[3] = evalIndexValuesAll + (dim + 1) * parallelIndices2[3];

	// double *evalIndexValues3[4];
	// evalIndexValues3[0] = evalIndexValuesAll + (dim + 1) * parallelIndices3[0];
	// evalIndexValues3[1] = evalIndexValuesAll + (dim + 1) * parallelIndices3[1];
	// evalIndexValues3[2] = evalIndexValuesAll + (dim + 1) * parallelIndices3[2];
	// evalIndexValues3[3] = evalIndexValuesAll + (dim + 1) * parallelIndices3[3];

	// double *evalIndexValues4[4];
	// evalIndexValues4[0] = evalIndexValuesAll + (dim + 1) * parallelIndices4[0];
	// evalIndexValues4[1] = evalIndexValuesAll + (dim + 1) * parallelIndices4[1];
	// evalIndexValues4[2] = evalIndexValuesAll + (dim + 1) * parallelIndices4[2];
	// evalIndexValues4[3] = evalIndexValuesAll + (dim + 1) * parallelIndices4[3];

	// double *evalIndexValues5[4];
	// evalIndexValues5[0] = evalIndexValuesAll + (dim + 1) * parallelIndices5[0];
	// evalIndexValues5[1] = evalIndexValuesAll + (dim + 1) * parallelIndices5[1];
	// evalIndexValues5[2] = evalIndexValuesAll + (dim + 1) * parallelIndices5[2];
	// evalIndexValues5[3] = evalIndexValuesAll + (dim + 1) * parallelIndices5[3];

	// double *evalIndexValues6[4];
	// evalIndexValues6[0] = evalIndexValuesAll + (dim + 1) * parallelIndices6[0];
	// evalIndexValues6[1] = evalIndexValuesAll + (dim + 1) * parallelIndices6[1];
	// evalIndexValues6[2] = evalIndexValuesAll + (dim + 1) * parallelIndices6[2];
	// evalIndexValues6[3] = evalIndexValuesAll + (dim + 1) * parallelIndices6[3];


	//for faster index flattening, last element is for padding
	uint32_t *intermediates[4];
	intermediates[0] = intermediatesAll + (dim + 1) * parallelIndices[0];
	intermediates[1] = intermediatesAll + (dim + 1) * parallelIndices[1];
	intermediates[2] = intermediatesAll + (dim + 1) * parallelIndices[2];
	intermediates[3] = intermediatesAll + (dim + 1) * parallelIndices[3];

	//pipe
	// uint32_t *intermediates2[4];
	// intermediates2[0] = intermediatesAll + (dim + 1) * parallelIndices2[0];
	// intermediates2[1] = intermediatesAll + (dim + 1) * parallelIndices2[1];
	// intermediates2[2] = intermediatesAll + (dim + 1) * parallelIndices2[2];
	// intermediates2[3] = intermediatesAll + (dim + 1) * parallelIndices2[3];

	// uint32_t *intermediates3[4];
	// intermediates3[0] = intermediatesAll + (dim + 1) * parallelIndices3[0];
	// intermediates3[1] = intermediatesAll + (dim + 1) * parallelIndices3[1];
	// intermediates3[2] = intermediatesAll + (dim + 1) * parallelIndices3[2];
	// intermediates3[3] = intermediatesAll + (dim + 1) * parallelIndices3[3];

	// uint32_t *intermediates4[4];
	// intermediates4[0] = intermediatesAll + (dim + 1) * parallelIndices4[0];
	// intermediates4[1] = intermediatesAll + (dim + 1) * parallelIndices4[1];
	// intermediates4[2] = intermediatesAll + (dim + 1) * parallelIndices4[2];
	// intermediates4[3] = intermediatesAll + (dim + 1) * parallelIndices4[3];

	// uint32_t *intermediates5[4];
	// intermediates5[0] = intermediatesAll + (dim + 1) * parallelIndices5[0];
	// intermediates5[1] = intermediatesAll + (dim + 1) * parallelIndices5[1];
	// intermediates5[2] = intermediatesAll + (dim + 1) * parallelIndices5[2];
	// intermediates5[3] = intermediatesAll + (dim + 1) * parallelIndices5[3];

	// uint32_t *intermediates6[4];
	// intermediates6[0] = intermediatesAll + (dim + 1) * parallelIndices6[0];
	// intermediates6[1] = intermediatesAll + (dim + 1) * parallelIndices6[1];
	// intermediates6[2] = intermediatesAll + (dim + 1) * parallelIndices6[2];
	// intermediates6[3] = intermediatesAll + (dim + 1) * parallelIndices6[3];


	uint32_t indexFlat[4];
	double phiEval[4];

	//pipe
	// uint32_t indexFlat2[4];
	// double phiEval2[4];

	// uint32_t indexFlat3[4];
	// double phiEval3[4];

	// uint32_t indexFlat4[4];
	// double phiEval4[4];

	// uint32_t indexFlat5[4];
	// double phiEval5[4];

	// uint32_t indexFlat6[4];
	// double phiEval6[4];
	
	// X86Vectorized::calculateIndexVectorized6(dim, nextIterationToRecalc, 
	// 					 dataTuplePtr, 
	// 					 dataTuplePtr2,
	// 					 dataTuplePtr3,
	// 					 dataTuplePtr4,
	// 					 dataTuplePtr5,
	// 					 dataTuplePtr6,
	// 					 hInversePtr, 
	// 					 intermediates, 
	// 					 intermediates2, 
	// 					 intermediates3, 
	// 					 intermediates4, 
	// 					 intermediates5, 
	// 					 intermediates6, 
	// 					 evalIndexValues,
	// 					 evalIndexValues2,
	// 					 evalIndexValues3,
	// 					 evalIndexValues4,
	// 					 evalIndexValues5,
	// 					 evalIndexValues6,
	// 					 indexFlat, 
	// 					 indexFlat2, 
	// 					 indexFlat3, 
	// 					 indexFlat4, 
	// 					 indexFlat5, 
	// 					 indexFlat6, 
	// 					 phiEval,
	// 					 phiEval2,
	// 					 phiEval3,
	// 					 phiEval4,
	// 					 phiEval5,
	// 					 phiEval6
	// 					 );

	X86Vectorized::calculateIndexVectorized(dim, nextIterationToRecalc, 
						dataTuplePtr, hInversePtr, 
						intermediates, 
						evalIndexValues,
						indexFlat, phiEval
						);

	double surplus[4];
	surplus[0] = levelArrayContinuous[indexFlat[0]];
	surplus[1] = levelArrayContinuous[indexFlat[1]];
	surplus[2] = levelArrayContinuous[indexFlat[2]];
	surplus[3] = levelArrayContinuous[indexFlat[3]];

	// double surplus2[4];
	// surplus2[0] = levelArrayContinuous[indexFlat2[0]];
	// surplus2[1] = levelArrayContinuous[indexFlat2[1]];
	// surplus2[2] = levelArrayContinuous[indexFlat2[2]];
	// surplus2[3] = levelArrayContinuous[indexFlat2[3]];

	// double surplus3[4];
	// surplus3[0] = levelArrayContinuous[indexFlat3[0]];
	// surplus3[1] = levelArrayContinuous[indexFlat3[1]];
	// surplus3[2] = levelArrayContinuous[indexFlat3[2]];
	// surplus3[3] = levelArrayContinuous[indexFlat3[3]];

	// double surplus4[4];
	// surplus4[0] = levelArrayContinuous[indexFlat4[0]];
	// surplus4[1] = levelArrayContinuous[indexFlat4[1]];
	// surplus4[2] = levelArrayContinuous[indexFlat4[2]];
	// surplus4[3] = levelArrayContinuous[indexFlat4[3]];

	// double surplus5[4];
	// surplus5[0] = levelArrayContinuous[indexFlat5[0]];
	// surplus5[1] = levelArrayContinuous[indexFlat5[1]];
	// surplus5[2] = levelArrayContinuous[indexFlat5[2]];
	// surplus5[3] = levelArrayContinuous[indexFlat5[3]];

	// double surplus6[4];
	// surplus6[0] = levelArrayContinuous[indexFlat6[0]];
	// surplus6[1] = levelArrayContinuous[indexFlat6[1]];
	// surplus6[2] = levelArrayContinuous[indexFlat6[2]];
	// surplus6[3] = levelArrayContinuous[indexFlat6[3]];

	for (size_t innerIndex = 0; innerIndex < 4; innerIndex++) {

	  //TODO: rework ifs (at least remove inner)
	  size_t parallelIndex = parallelIndices[innerIndex];
	  if (!std::isnan(surplus[innerIndex])) {
	    double partialSurplus = 0.0;
	    if (dataIndexBase + parallelIndex < end_index_data 
		&& parallelIndex < X86VECTORIZED_PARALLEL_DATA_POINTS) {
	      partialSurplus = phiEval[innerIndex] * alpha[dataIndexBase + parallelIndex]; 
	      
	      size_t localIndexFlat = indexFlat[innerIndex];
// #pragma omp atomic
// 	      levelArray[localIndexFlat] += partialSurplus;
	      
#pragma omp atomic
	      levelArrayContinuous[localIndexFlat] += partialSurplus;
	      
	    }

	    nextIterationToRecalcReferences[parallelIndex] = allSubspaces[GET_NEXTDIFF_LOOP(levelIndex)];
	    levelIndices[parallelIndex] += this->subspaceSize;
	  } else {
#if X86VECTORIZED_ENABLE_SUBSPACE_SKIPPING == 1
	    //skip to next relevant subspace
	    nextIterationToRecalcReferences[parallelIndex] = allSubspaces[GET_JUMPDIFF_LOOP(levelIndex)];
	    levelIndices[parallelIndex] = allSubspaces[GET_NEXT_LOOP(levelIndex)];
#else
	    nextIterationToRecalcReferences[parallelIndex] = allSubspaces[GET_NEXTDIFF_LOOP(levelIndex)];
	    levelIndices[parallelIndex] += this->subspaceSize;
#endif
	  }
	} // end innerIndex

	//pipe
// 	for (size_t innerIndex = 0; innerIndex < 4; innerIndex++) {
// 	  size_t parallelIndex = parallelIndices2[innerIndex];
// 	  if (!std::isnan(surplus2[innerIndex])) {
// 	    double partialSurplus = 0.0;
// 	    if (dataIndexBase + parallelIndex < end_index_data 
// 		&& parallelIndex < X86VECTORIZED_PARALLEL_DATA_POINTS) {
// 	      partialSurplus = phiEval2[innerIndex] * alpha[dataIndexBase + parallelIndex]; 	      
// 	      size_t localIndexFlat = indexFlat2[innerIndex];
// #pragma omp atomic
// 	      levelArrayContinuous[localIndexFlat] += partialSurplus;	      
// 	    }
// 	    nextIterationToRecalcReferences[parallelIndex] = allSubspaces[GET_NEXTDIFF_LOOP(levelIndex)];
// 	    levelIndices[parallelIndex] += this->subspaceSize;
// 	  } else {
// #if X86VECTORIZED_ENABLE_SUBSPACE_SKIPPING == 1
// 	    //skip to next relevant subspace
// 	    nextIterationToRecalcReferences[parallelIndex] = allSubspaces[GET_JUMPDIFF_LOOP(levelIndex)];
// 	    levelIndices[parallelIndex] = allSubspaces[GET_NEXT_LOOP(levelIndex)];
// #else
// 	    nextIterationToRecalcReferences[parallelIndex] = allSubspaces[GET_NEXTDIFF_LOOP(levelIndex)];
// 	    levelIndices[parallelIndex] += this->subspaceSize;
// #endif
// 	  }
// 	} // end innerIndex

// 	//pipe
// 	for (size_t innerIndex = 0; innerIndex < 4; innerIndex++) {
// 	  size_t parallelIndex = parallelIndices3[innerIndex];
// 	  if (!std::isnan(surplus3[innerIndex])) {
// 	    double partialSurplus = 0.0;
// 	    if (dataIndexBase + parallelIndex < end_index_data 
// 		&& parallelIndex < X86VECTORIZED_PARALLEL_DATA_POINTS) {
// 	      partialSurplus = phiEval3[innerIndex] * alpha[dataIndexBase + parallelIndex]; 	      
// 	      size_t localIndexFlat = indexFlat3[innerIndex];
// #pragma omp atomic
// 	      levelArrayContinuous[localIndexFlat] += partialSurplus;	      
// 	    }
// 	    nextIterationToRecalcReferences[parallelIndex] = allSubspaces[GET_NEXTDIFF_LOOP(levelIndex)];
// 	    levelIndices[parallelIndex] += this->subspaceSize;
// 	  } else {
// #if X86VECTORIZED_ENABLE_SUBSPACE_SKIPPING == 1
// 	    //skip to next relevant subspace
// 	    nextIterationToRecalcReferences[parallelIndex] = allSubspaces[GET_JUMPDIFF_LOOP(levelIndex)];
// 	    levelIndices[parallelIndex] = allSubspaces[GET_NEXT_LOOP(levelIndex)];
// #else
// 	    nextIterationToRecalcReferences[parallelIndex] = allSubspaces[GET_NEXTDIFF_LOOP(levelIndex)];
// 	    levelIndices[parallelIndex] += this->subspaceSize;
// #endif
// 	  }
// 	} // end innerIndex

// 	//pipe
// 	for (size_t innerIndex = 0; innerIndex < 4; innerIndex++) {
// 	  size_t parallelIndex = parallelIndices4[innerIndex];
// 	  if (!std::isnan(surplus4[innerIndex])) {
// 	    double partialSurplus = 0.0;
// 	    if (dataIndexBase + parallelIndex < end_index_data 
// 		&& parallelIndex < X86VECTORIZED_PARALLEL_DATA_POINTS) {
// 	      partialSurplus = phiEval4[innerIndex] * alpha[dataIndexBase + parallelIndex]; 	      
// 	      size_t localIndexFlat = indexFlat4[innerIndex];
// #pragma omp atomic
// 	      levelArrayContinuous[localIndexFlat] += partialSurplus;	      
// 	    }
// 	    nextIterationToRecalcReferences[parallelIndex] = allSubspaces[GET_NEXTDIFF_LOOP(levelIndex)];
// 	    levelIndices[parallelIndex] += this->subspaceSize;
// 	  } else {
// #if X86VECTORIZED_ENABLE_SUBSPACE_SKIPPING == 1
// 	    //skip to next relevant subspace
// 	    nextIterationToRecalcReferences[parallelIndex] = allSubspaces[GET_JUMPDIFF_LOOP(levelIndex)];
// 	    levelIndices[parallelIndex] = allSubspaces[GET_NEXT_LOOP(levelIndex)];
// #else
// 	    nextIterationToRecalcReferences[parallelIndex] = allSubspaces[GET_NEXTDIFF_LOOP(levelIndex)];
// 	    levelIndices[parallelIndex] += this->subspaceSize;
// #endif
// 	  }
// 	} // end innerIndex
	
// 	//pipe
// 	for (size_t innerIndex = 0; innerIndex < 4; innerIndex++) {
// 	  size_t parallelIndex = parallelIndices5[innerIndex];
// 	  if (!std::isnan(surplus5[innerIndex])) {
// 	    double partialSurplus = 0.0;
// 	    if (dataIndexBase + parallelIndex < end_index_data 
// 		&& parallelIndex < X86VECTORIZED_PARALLEL_DATA_POINTS) {
// 	      partialSurplus = phiEval5[innerIndex] * alpha[dataIndexBase + parallelIndex]; 	      
// 	      size_t localIndexFlat = indexFlat5[innerIndex];
// #pragma omp atomic
// 	      levelArrayContinuous[localIndexFlat] += partialSurplus;	      
// 	    }
// 	    nextIterationToRecalcReferences[parallelIndex] = allSubspaces[GET_NEXTDIFF_LOOP(levelIndex)];
// 	    levelIndices[parallelIndex] += this->subspaceSize;
// 	  } else {
// #if X86VECTORIZED_ENABLE_SUBSPACE_SKIPPING == 1
// 	    //skip to next relevant subspace
// 	    nextIterationToRecalcReferences[parallelIndex] = allSubspaces[GET_JUMPDIFF_LOOP(levelIndex)];
// 	    levelIndices[parallelIndex] = allSubspaces[GET_NEXT_LOOP(levelIndex)];
// #else
// 	    nextIterationToRecalcReferences[parallelIndex] = allSubspaces[GET_NEXTDIFF_LOOP(levelIndex)];
// 	    levelIndices[parallelIndex] += this->subspaceSize;
// #endif
// 	  }
// 	} // end innerIndex

// 	//pipe
// 	for (size_t innerIndex = 0; innerIndex < 4; innerIndex++) {
// 	  size_t parallelIndex = parallelIndices6[innerIndex];
// 	  if (!std::isnan(surplus6[innerIndex])) {
// 	    double partialSurplus = 0.0;
// 	    if (dataIndexBase + parallelIndex < end_index_data 
// 		&& parallelIndex < X86VECTORIZED_PARALLEL_DATA_POINTS) {
// 	      partialSurplus = phiEval6[innerIndex] * alpha[dataIndexBase + parallelIndex]; 	      
// 	      size_t localIndexFlat = indexFlat6[innerIndex];
// #pragma omp atomic
// 	      levelArrayContinuous[localIndexFlat] += partialSurplus;	      
// 	    }
// 	    nextIterationToRecalcReferences[parallelIndex] = allSubspaces[GET_NEXTDIFF_LOOP(levelIndex)];
// 	    levelIndices[parallelIndex] += this->subspaceSize;
// 	  } else {
// #if X86VECTORIZED_ENABLE_SUBSPACE_SKIPPING == 1
// 	    //skip to next relevant subspace
// 	    nextIterationToRecalcReferences[parallelIndex] = allSubspaces[GET_JUMPDIFF_LOOP(levelIndex)];
// 	    levelIndices[parallelIndex] = allSubspaces[GET_NEXT_LOOP(levelIndex)];
// #else
// 	    nextIterationToRecalcReferences[parallelIndex] = allSubspaces[GET_NEXTDIFF_LOOP(levelIndex)];
// 	    levelIndices[parallelIndex] += this->subspaceSize;
// #endif
// 	  }
// 	} // end innerIndex

      } // end parallel
    } // end iterate subspaces
  } // end iterate chunks

  //delete indexPtrAll;
  delete evalIndexValuesAll;
  delete intermediatesAll;

#pragma omp barrier

  if (tid == 0) {
    this->unflatten(result);
  }

}
