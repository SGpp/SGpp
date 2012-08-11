/******************************************************************************
* Copyright (C) 2011 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef OPERATIONMULTIPLEEVALVECTORIZED_HPP
#define OPERATIONMULTIPLEEVALVECTORIZED_HPP

#include "base/datatypes/DataVector.hpp"
#include "base/datatypes/DataMatrix.hpp"
#include "base/operation/OperationMatrix.hpp"

#define CHUNKDATAPOINTS_X86 24 // must be divide-able by 24
#define CHUNKGRIDPOINTS_X86 12

namespace sg
{
namespace parallel
{

/**
 * @brief Interface for multiplication with Matrices @f$B@f$ and @f$B^T@f$.
 *
 * If there are @f$N@f$ basis functions, @f$\{\varphi(\vec{x})\}_{i=1,\ldots,N}@f$ and @f$m@f$ data points
 *
 * This class defines an interface similar to OperationMultipleEval in order to support SIMD architectures
 * for datadriven task (multiple function evaluations, classification, regression). Target
 * architectures may be Intel SSE, Intel AVX, nVidia CUDA, OpenCL.
 */
class OperationMultipleEvalVectorized
{
protected:
	void adaptDatasetBoundaries();
	/**
	 * @brief calcOpenMPLoopDistribution
	 *
	 *@todo proper documentation
	 * call this from within a parallel section
	 *
	 * @param processStart
	 * @param processEnd
	 * @param chunkSize
	 * @param start
	 * @param end
	 */
	void calcOpenMPLoopDistribution(int processStart, int processEnd, int chunkSize, size_t *start, size_t *end);


	/// Pointer to the dataset that should be evaluated on the grid
	sg::base::DataMatrix* dataset_;
	/// Member to store the sparse grid's levels for better vectorization
	sg::base::DataMatrix* level_;
	/// Member to store the sparse grid's indices for better vectorization
	sg::base::DataMatrix* index_;

	int m_gridFrom;
	int m_gridTo;
	int m_datasetFrom;
	int m_datasetTo;

public:
	/**
	 * Constructor
	 *
	 * @param dataset data set that should be evaluated on the sparse grid
	 */
	OperationMultipleEvalVectorized(sg::base::DataMatrix* dataset) {
		this->dataset_ = dataset;
		this->level_ = NULL;
		this->index_ = NULL;
	}

	/**
	 * Destructor
	 *
	 * cleans up level_ and index_ members
	 */
	virtual ~OperationMultipleEvalVectorized() {
		if (this->level_ != NULL)
			delete this->level_;

		if (this->index_ != NULL)
			delete this->index_;
	}

	/**
	 * Determines the size of stored dataset, including event. padding.
	 *
	 * @return returns the size of the stored dataset
	 */
	size_t getDatasetSize()
	{
		return this->dataset_->getNrows();
	}

	/**
	 * Multiplication of @f$B^T@f$ with vector @f$\alpha@f$
	 *
	 * IMPORTANT REMARK:
	 * In order to use this routine you have to keep following points in mind (for multVectorized and multTransposeVectorized):
 	 * 	- data MUST a have even number of points AND it must be transposed
 	 *  - result MUST have the same size as data points that should be evaluated
	 *
	 * @param alpha vector, to which @f$B@f$ is applied. Typically the coefficient vector
	 * @param result the result vector of the matrix vector multiplication
	 */
	virtual double multVectorized(sg::base::DataVector& alpha, sg::base::DataVector& result) = 0;

	/**
	 * Multiplication of @f$B@f$ with vector @f$\alpha@f$
	 *
	 * IMPORTANT REMARK:
	 * In order to use this routine you have to keep following points in mind (for multVectorized and multTransposeVectorized):
 	 * 	- data MUST a have even number of points AND it must be transposed
 	 *  - result MUST have the same size as data points that should be evaluated
	 *
	 * @param source vector, to which @f$B^T@f$ is applied. Typically the coefficient vector
	 * @param result the result vector of the matrix vector multiplication
	 */
	virtual double multTransposeVectorized(sg::base::DataVector& source, sg::base::DataVector& result) = 0;

	/**
	 * rebuilds the DataMatrix for Level and Index in Derivatives
	 * needed for vectorization.
	 */
	virtual void rebuildLevelAndIndex() = 0;

	/**
	 * @brief updates the compute boundaries for the grid, after this has been resized
	 *
	 * @todo for now, the default implementation does nothing. perhaps remove default implementation.
	 * it would be an idea to integrate this with rebuildLevelAndIndex
	 *
	 * @param gridFrom
	 * @param gridTo
	 */
	virtual void updateGridComputeBoundaries(int gridFrom, int gridTo){}
};

}
}

#endif /* OPERATIONMULTIPLEEVALVECTORIZED_HPP */
