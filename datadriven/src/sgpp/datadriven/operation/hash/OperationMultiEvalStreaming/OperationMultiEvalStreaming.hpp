// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <omp.h>

#include <sgpp/base/operation/hash/OperationMultipleEval.hpp>
#include <sgpp/base/tools/SGppStopwatch.hpp>
#include <sgpp/base/exception/operation_exception.hpp>

#include <sgpp/globaldef.hpp>

namespace SGPP {
namespace datadriven {

class OperationMultiEvalStreaming: public base::OperationMultipleEval {
protected:
	SGPP::base::DataMatrix preparedDataset;
	/// Member to store the sparse grid's levels for better vectorization
	SGPP::base::DataMatrix* level_ = nullptr;
	/// Member to store the sparse grid's indices for better vectorization
	SGPP::base::DataMatrix* index_ = nullptr;
	/// Timer object to handle time measurements
	SGPP::base::SGppStopwatch myTimer_;

	base::GridStorage* storage;

	float_t duration;
public:

	OperationMultiEvalStreaming(base::Grid& grid, base::DataMatrix& dataset);

	~OperationMultiEvalStreaming();

	size_t getChunkGridPoints();

	size_t getChunkDataPoints();

	void mult(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result)
			override;

	void multTranspose(SGPP::base::DataVector& source,
			SGPP::base::DataVector& result) override;

	void prepare() override;

	float_t getDuration();

private:
	void getPartitionSegment(size_t start, size_t end, size_t segmentCount,
			size_t segmentNumber, size_t* segmentStart, size_t* segmentEnd,
			size_t blockSize);

	size_t padDataset(SGPP::base::DataMatrix& dataset);

	void getOpenMPPartitionSegment(size_t start, size_t end,
			size_t* segmentStart, size_t* segmentEnd, size_t blocksize);

	void multImpl(SGPP::base::DataMatrix* level, SGPP::base::DataMatrix* index,
			SGPP::base::DataMatrix* dataset,
			SGPP::base::DataVector& alpha, SGPP::base::DataVector& result,
			const size_t start_index_grid, const size_t end_index_grid,
			const size_t start_index_data, const size_t end_index_data);

	void multTransposeImpl(SGPP::base::DataMatrix* level,
			SGPP::base::DataMatrix* index, SGPP::base::DataMatrix* dataset,
			SGPP::base::DataVector& source, SGPP::base::DataVector& result,
			const size_t start_index_grid, const size_t end_index_grid,
			const size_t start_index_data, const size_t end_index_data);

	void recalculateLevelAndIndex();
};

}
}
