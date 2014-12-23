#pragma once

#include <iostream>
#include <vector>
#include <map>

//#include "PseudoHashMap.hpp"

#include "omp.h"

#include "../../OperationMultipleEvalSubspace/vectorized/X86VectorizedParameters.hpp"
#include "../../OperationMultipleEvalSubspace/vectorized/X86VectorizedSubspaceNode.hpp"
#include "../../OperationMultipleEvalSubspace/vectorized/X86VectorizedTypes.hpp"
#include "FastKernelBase.hpp"
#include "X86VectorizedTypes.hpp"

using namespace std;

using namespace sg::base;
using namespace sg::parallel;

namespace sg {
  namespace parallel {

    class X86Vectorized : public FastKernelBase {
    private:

      size_t subspaceSize = -1;

      // list of offsets for subspace structure
      size_t levelOffset;
      size_t hInverseOffset;
      size_t nextOffset;
      size_t flatLevelOffset;
      size_t nextDiffOffset;
      size_t jumpDiffOffset;
      size_t flatLevelPointerOffset;

      map<uint32_t, X86VectorizedSubspaceNode> allLevelsMap;

      DataMatrix *dataset = nullptr;
      size_t dim = -1;
      size_t maxLevel = 0;

      uint32_t *allSubspaces = nullptr;
      size_t subspaceCount = -1;


      /// Pointer to the grid's gridstorage object
      sg::base::GridStorage* storage = nullptr;

      //double **flatLevels = nullptr;  
      uint32_t *levelGridPoints = nullptr;

      //all coefficients of all subspaces are stored here
      double *allSurplusses = nullptr;
      map<uint32_t, double *> flatLevelArrayMap;
      uint32_t totalRegularGridPoints = -1;

      void prepareSubspaceIterator();

      void createFlatStorage();

      static inline uint32_t calculateIndexComponent(size_t dim, double unadjusted) {
	//implies flooring
	uint32_t rounded = static_cast<uint32_t>(unadjusted);	

	uint32_t mask = 0x1;
	uint32_t sign = mask ^ (mask & rounded);
	
	uint32_t componentIndex = rounded + sign;
	return componentIndex;
      }

      static inline void calculateIndex(size_t dim, size_t nextIterationToRecalc, 
					double *dataTuplePtr, 
					uint32_t *hInversePtr, size_t *indexPtr) {
	for (size_t i = nextIterationToRecalc; i < dim; i++) {
	  double unadjusted = dataTuplePtr[i] * hInversePtr[i];

          //implies flooring
	  uint32_t rounded = static_cast<uint32_t>(unadjusted);	

          uint32_t mask = 0x1;
	  uint32_t sign = mask ^ (mask & rounded);
	
	  uint32_t componentIndex = rounded + sign;
          indexPtr[i] = componentIndex;
	}
      }

      #include "../../OperationMultipleEvalSubspace/vectorized/X86Vectorized_calculateIndexVectorized.hpp"

      static inline double evaluateBase(size_t dim, size_t nextIterationToRecalc,
					double *dataTuplePtr, double *evalIndexValues, 
					uint32_t *hInversePtr, size_t *indexPtr) {
	double phiEval = evalIndexValues[nextIterationToRecalc];
	//prepare the values for the individual components
	for (size_t i = nextIterationToRecalc; i < dim; i++) {
	  double phi1DEval = hInversePtr[i] * dataTuplePtr[i] - indexPtr[i];
	  phi1DEval = max(0.0, 1.0 - abs(phi1DEval));
	  phiEval *= phi1DEval;
	  evalIndexValues[i + 1] = phiEval;
	}
	return phiEval;
      }

      #include "../../OperationMultipleEvalSubspace/vectorized/X86Vectorized_evaluateBaseVectorized.hpp"
       
    public:
      static const KernelType kernelType = Standard;

      X86Vectorized(GridStorage* storage, 
		  DataMatrix* dataset,
		  size_t dim);

      ~X86Vectorized();

      void prepare();

      void setCoefficients(DataVector &surplusVector);

      void unflatten(DataVector &result);

      static uint32_t flattenIndex(size_t dim, DataVector &maxIndices, DataVector &index);

      void setSurplus(DataVector &level, DataVector &maxIndices, DataVector &index, double value);

      void getSurplus(DataVector &level, DataVector &maxIndices, DataVector &index, double &value, bool &isVirtual);

      uint32_t flattenLevel(size_t dim, size_t maxLevel, DataVector &level);

      void multImpl(sg::base::DataVector &alpha,
		    sg::base::DataVector &result,
		    //const size_t start_index_grid,
		    //const size_t end_index_grid,
		    const size_t start_index_data,
		    const size_t end_index_data);

      void multTransposeImpl(sg::base::DataVector &source,
			     sg::base::DataVector &result,
			     //as only the non-zero grid points are evaluated,
			     //blocking over grid make no sense
			     //const size_t start_index_grid,
			     //const size_t end_index_grid,
			     const size_t start_index_data,
			     const size_t end_index_data);

      void padDataset();

      static inline uint32_t flattenIndex(uint32_t *intermediates, size_t dim, 
					//size_t *maxIndicesPtr, size_t *indexPtr, size_t toRecalc) {
					  uint32_t *maxIndicesPtr, uint32_t *indexPtr, size_t toRecalc) {
	uint32_t indexFlat = intermediates[toRecalc]; // toRecalc 0 -> indexFlat 0
	for (size_t i = toRecalc; i < dim; i++) {
	  int actualDirectionGridPoints = maxIndicesPtr[i];
	  actualDirectionGridPoints >>= 1;
	  indexFlat *= actualDirectionGridPoints;
	  uint32_t actualIndex = indexPtr[i];
	  actualIndex >>= 1; //divide index by 2, skip even indices
	  indexFlat += actualIndex;
	  intermediates[i + 1] = indexFlat;
	}
	return indexFlat;
      }

      #include "../../OperationMultipleEvalSubspace/vectorized/X86Vectorized_flattenIndexVectorized.hpp"

      static inline uint32_t flattenIndex(double *intermediates, size_t dim, 
					uint32_t *maxIndicesPtr, size_t *indexPtr, size_t toRecalc) {
	uint32_t indexFlat = intermediates[toRecalc]; // toRecalc 0 -> indexFlat 0
	for (size_t i = toRecalc; i < dim; i++) {
	  int actualDirectionGridPoints = maxIndicesPtr[i];
	  actualDirectionGridPoints >>= 1;
	  indexFlat *= actualDirectionGridPoints;
	  size_t actualIndex = indexPtr[i];
	  actualIndex >>= 1; //divide index by 2, skip even indices
	  indexFlat += actualIndex;
	  intermediates[i + 1] = indexFlat;
	}    
	return indexFlat;
      }

      static inline uint32_t flattenLevel(size_t dim, size_t maxLevel, uint32_t *levelPtr) {
	uint32_t levelFlat = 0;
	//levelFlat += level.get(dim - 1);
	levelFlat += levelPtr[dim - 1];
	// loop terminates at -1
	for (int i = dim - 2; i >= 0; i--) {
	  levelFlat *= maxLevel;
	  //levelFlat += level.get(i);
	  levelFlat += levelPtr[i];
	}
	return levelFlat;
      }

      size_t getAlignment() {
	return X86VECTORIZED_PARALLEL_DATA_POINTS;
      }
    };
  }
}
