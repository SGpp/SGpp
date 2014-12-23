#pragma once

#include <iostream>
#include <vector>
#include <map>

//////////////////////////////////////////////////////////////////////
// Caution: Subspace-skipping is disabled by default for this kernel
//////////////////////////////////////////////////////////////////////

#include "omp.h"

#include "datadriven/operation/OperationMultipleEvalSubspace/AbstractOperationMultipleEvalSubspace.hpp"
#include "datadriven/operation/OperationMultipleEvalSubspace/SubspaceNode.hpp"
#include "datadriven/operation/OperationMultipleEvalSubspace/simple/OperationMultipleEvalSubspaceSimpleParameters.hpp"

namespace sg {
namespace datadriven {

class OperationMultipleEvalSubspaceSimple: public AbstractOperationMultipleEvalSubspace {
private:

	size_t dim = -1;
	size_t maxLevel = 0;

	size_t *allSubspaces = NULL;
	size_t subspaceCount = -1;
	size_t subspaceSize = -1;

	size_t totalGridPoints = 0;
	double *allSurplusses = nullptr;
	map<uint32_t, uint32_t> allSurplussesIndexMap;

	void prepareSubspaceIterator();

	void createFlatStorage();

	void setSurplus(DataVector &level, DataVector &maxIndices, DataVector &index, double value);

	void getSurplus(DataVector &level, DataVector &maxIndices, DataVector &index, double &value, bool &isVirtual);

	void setCoefficients(DataVector &surplusVector);

	void unflatten(DataVector &result);

	size_t flattenIndex(size_t dim, DataVector &maxIndices, DataVector &index);

	size_t flattenIndex(double *intermediates, size_t dim, size_t *maxIndicesPtr, size_t *indexPtr,
			size_t toRecalc);

	size_t flattenLevel(size_t dim, size_t maxLevel, DataVector &level);

	static inline size_t calculateIndexComponent(size_t dim, double unadjusted) {
		//implies flooring
		size_t rounded = static_cast<size_t>(unadjusted);

		size_t mask = 0x1;
		size_t sign = mask ^ (mask & rounded);

		size_t componentIndex = rounded + sign;
		return componentIndex;
	}

public:

	OperationMultipleEvalSubspaceSimple(Grid &grid, DataMatrix &dataset);

	~OperationMultipleEvalSubspaceSimple();

	void prepare() override;

	void multImpl(sg::base::DataVector &alpha, sg::base::DataVector &result, const size_t start_index_data,
			const size_t end_index_data) override;

	void multTransposeImpl(sg::base::DataVector &source, sg::base::DataVector &result, const size_t start_index_data,
			const size_t end_index_data) override;

	void padDataset() override;

	size_t getAlignment() override;

	std::string getImplementationName() override;
};
}
}
