#ifndef X86BASELINE_HPP
#define X86BASELINE_HPP

#include <iostream>
#include <vector>
#include <map>

//#include "PseudoHashMap.hpp"

#include "omp.h"

#include "../../OperationMultipleEvalSubspace/baseline/X86BaselineParameters.hpp"
#include "../../OperationMultipleEvalSubspace/SubspaceNode.hpp"
#include "FastKernelBase.hpp"

using namespace std;

using namespace sg::base;
using namespace sg::parallel;

namespace sg {
  namespace parallel {

    class X86Baseline : public FastKernelBase {
    private:

      map<size_t, SubspaceNode> allLevelsMap;

      DataMatrix *dataset = nullptr;
      size_t dim = -1;
      size_t maxLevel = 0;

      size_t *allSubspaces = nullptr;
      size_t subspaceCount = -1;
      size_t subspaceSize = -1;

      /// Pointer to the grid's gridstorage object
      sg::base::GridStorage* storage = nullptr;

//      double **flatLevels = nullptr;
//      size_t *levelGridPoints = nullptr;

      size_t totalGridPoints = 0;
      double *allSurplusses = nullptr;
      map<uint32_t, uint32_t> allSurplussesIndexMap;

      void prepareSubspaceIterator();

      void createFlatStorage();

      static inline size_t calculateIndexComponent(size_t dim, double unadjusted) {
	//implies flooring
	size_t rounded = static_cast<size_t>(unadjusted);	

	size_t mask = 0x1;
	size_t sign = mask ^ (mask & rounded);
	
	size_t componentIndex = rounded + sign;
	return componentIndex;
      }
      
    public:
      static const KernelType kernelType = Standard;

      X86Baseline(GridStorage* storage, 
		  DataMatrix* dataset,
		  size_t dim);

      ~X86Baseline();

      void prepare();

      void setCoefficients(DataVector &surplusVector);

      void unflatten(DataVector &result);

      static size_t flattenIndex(size_t dim, DataVector &maxIndices, DataVector &index);

      void setSurplus(DataVector &level, DataVector &maxIndices, DataVector &index, double value);

      void getSurplus(DataVector &level, DataVector &maxIndices, DataVector &index, double &value, bool &isVirtual);

      size_t flattenLevel(size_t dim, size_t maxLevel, DataVector &level);

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

      static inline size_t flattenIndex(double *intermediates, size_t dim, 
					size_t *maxIndicesPtr, size_t *indexPtr, size_t toRecalc) {
	size_t indexFlat = intermediates[toRecalc]; // toRecalc 0 -> indexFlat 0
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

      static inline size_t flattenLevel(size_t dim, size_t maxLevel, size_t *levelPtr) {
	size_t levelFlat = 0;
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
	return X86BASELINE_PARALLEL_DATA_POINTS;
      }
    };
  }
}

#endif // X86Baseline_HPP
