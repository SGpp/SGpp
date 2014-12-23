#ifndef OPERATIONMULTIEVALLINEARNOBOUNDARYHASHGRIDCPU_HPP
#define OPERATIONMULTIEVALLINEARNOBOUNDARYHASHGRIDCPU_HPP

#include <omp.h>

#include "base/operation/OperationMultipleEval.hpp"
#include "base/tools/SGppStopwatch.hpp"
#include "base/exception/operation_exception.hpp"

namespace sg {
namespace datadriven {

class OperationMultiEvalLinearNoBoundaryHashGridCPU: public base::OperationMultipleEval {
protected:
	/// Member to store the sparse grid's levels for better vectorization
	sg::base::DataMatrix* level_ = nullptr;
	/// Member to store the sparse grid's indices for better vectorization
	sg::base::DataMatrix* index_ = nullptr;
	/// Timer object to handle time measurements
	sg::base::SGppStopwatch myTimer_;

	base::GridStorage *storage;

	double duration;
public:

	OperationMultiEvalLinearNoBoundaryHashGridCPU(base::Grid &grid, base::DataMatrix &dataset);

	~OperationMultiEvalLinearNoBoundaryHashGridCPU();

	size_t getChunkGridPoints();

	size_t getChunkDataPoints();

	void mult(sg::base::DataVector& alpha, sg::base::DataVector& result) override;

	void multTranspose(sg::base::DataVector& source, sg::base::DataVector& result) override;

	void recalculateLevelAndIndex();

	double getDuration();

private:
	void getPartitionSegment(size_t start, size_t end, size_t segmentCount, size_t segmentNumber, size_t* segmentStart,
			size_t* segmentEnd, size_t blockSize);

	size_t padDataset(sg::base::DataMatrix &dataset);

	void getOpenMPPartitionSegment(size_t start, size_t end, size_t* segmentStart, size_t* segmentEnd,
			size_t blocksize);

	void multImpl(sg::base::DataMatrix* level, sg::base::DataMatrix* index, sg::base::DataMatrix* dataset,
			sg::base::DataVector& alpha, sg::base::DataVector& result, const size_t start_index_grid,
			const size_t end_index_grid, const size_t start_index_data, const size_t end_index_data);

	void multTransposeImpl(sg::base::DataMatrix* level, sg::base::DataMatrix* index, sg::base::DataMatrix* dataset,
			sg::base::DataVector& source, sg::base::DataVector& result, const size_t start_index_grid,
			const size_t end_index_grid, const size_t start_index_data, const size_t end_index_data);
};

}
}
#endif // OPERATIONMULTIEVALLINEARNOBOUNDARYHASHGRIDCPU_HPP
