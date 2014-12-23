/* ****************************************************************************
 * Copyright (C) 2013 Technische Universitaet Muenchen                         *
 * This file is part of the SG++ project. For conditions of distribution and   *
 * use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
 **************************************************************************** */
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)
// @author Roman Karlstetter (karlstetter@mytum.de)
#ifndef ABSTRACTOPERATIONMULTIPLEEVALSUBSPACE_HPP
#define ABSTRACTOPERATIONMULTIPLEEVALSUBSPACE_HPP

#include "../OperationMultipleEvalSubspace/CommonParameters.hpp"
#include "base/grid/Grid.hpp"
#include "base/grid/GridStorage.hpp"
#include "base/tools/SGppStopwatch.hpp"
#include "datadriven/tools/PartitioningTool.hpp"
#include "base/operation/OperationMultipleEval.hpp"

namespace sg {
namespace datadriven {

class AbstractOperationMultipleEvalSubspace: public base::OperationMultipleEval {
protected:
	//should be set to true to avoid called "prepare()" with every "mult()" or "multTranspose()"
	bool isPrepared;
	base::GridStorage *storage;
private:
	base::SGppStopwatch timer;
	double duration;

public:
	AbstractOperationMultipleEvalSubspace(base::Grid &grid, base::DataMatrix &dataset) :
			base::OperationMultipleEval(grid, dataset), isPrepared(false), storage(grid.getStorage()), duration(-1.0) {
	}

	~AbstractOperationMultipleEvalSubspace() {
	}

	virtual void prepareExecution(size_t gridFrom = 0, size_t gridTo = std::numeric_limits<size_t>::max()) {

		this->prepare();

		if (gridTo == std::numeric_limits<size_t>::max()) {
			gridTo = this->storage->size();
		}

		this->resetKernel();
	}

	virtual void prepare() = 0;

	//TODO possible bug "=0;" results in undefined symbol when importing the python library
	//TODO does this have to be public?
	virtual void padDataset() {
	}

	//TODO where and when to do the padding? (currently in multTranspose done and reverted)

	virtual void multImpl(base::DataVector &alpha, base::DataVector &result, const size_t start_index_data,
			const size_t end_index_data) = 0;

	virtual void multTransposeImpl(sg::base::DataVector &source, sg::base::DataVector &result,
			const size_t start_index_data, const size_t end_index_data) = 0;

	void mult(sg::base::DataVector& alpha, sg::base::DataVector& result) override {
		if (!this->isPrepared) {
			this->padDataset();
			this->prepareExecution();
		}

		const size_t start_index_data = 0;
		//TODO handle transposed matrix? or use transposed matrix?
		const size_t end_index_data = this->dataset.getNrows();

		this->timer.start();
		result.setAll(0.0);

#pragma omp parallel
		{
			size_t start;
			size_t end;
			PartitioningTool::getOpenMPPartitionSegment(start_index_data, end_index_data, &start, &end,
					this->getAlignment());
			this->multImpl(alpha, result, start, end);
		}

		this->duration = this->timer.stop();
	}

	//TODO assumes padded "result" vector -> document somewhere or improve eval interface
	void multTranspose(sg::base::DataVector& source, sg::base::DataVector& result) override {
		bool resultResized = false;
		size_t originalResultSize = result.getSize();
		if (!this->isPrepared) {

			this->padDataset();
			this->prepareExecution();
			if (result.getSize() != this->dataset.getSize()) {
				result.resize(this->dataset.getSize());
				for (size_t i = originalResultSize; i < result.getSize(); i++) {
					result[i] = 0.0;
				}
				resultResized = true;
			}
		}

		const size_t start_index_data = 0;
		const size_t end_index_data = this->dataset.getNrows();

		this->timer.start();
		result.setAll(0.0);

#pragma omp parallel
		{
			size_t start;
			size_t end;
			PartitioningTool::getOpenMPPartitionSegment(start_index_data, end_index_data, &start, &end,
					this->getAlignment());
			this->multTransposeImpl(source, result, start, end);
		}

		if (resultResized) {
			result.resize(originalResultSize);
		}
		this->duration = this->timer.stop();
	}

	virtual size_t getAlignment() = 0;

	void resetKernel() {
	}

	static inline size_t getChunkGridPoints() {
		return 12;
	}
	static inline size_t getChunkDataPoints() {
		return 24; // must be divisible by 24
	}

};

}
}

#endif // ABSTRACTOPERATIONMULTIPLEEVALSUBSPACE_HPP
