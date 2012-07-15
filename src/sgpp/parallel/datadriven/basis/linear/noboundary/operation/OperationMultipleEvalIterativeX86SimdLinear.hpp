/******************************************************************************
* Copyright (C) 2010 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef OPERATIONMULTIPLEEVALITERATIVEX86SIMDLINEAR_HPP
#define OPERATIONMULTIPLEEVALITERATIVEX86SIMDLINEAR_HPP

#include "parallel/datadriven/operation/OperationMultipleEvalVectorized.hpp"
#include "base/grid/GridStorage.hpp"
#include "base/tools/SGppStopwatch.hpp"

namespace sg
{
namespace parallel
{

/**
 * This class implements OperationMultipleEvalVectorizedSP for grids with linear basis ansatzfunctions with or
 * without boundaries with double precision coefficients.
 *
 * However, in this case highly efficient vector code (SSE or AVX instructions) are generated
 * to implement a iterative OperationB version. In addition cache blocking is used
 * in order to assure a most efficient cache usage.
 *
 * IMPORTANT REMARK:
 * In order to use this routine you have to keep following points in mind (for multVectorized and multTransposeVectorized):
 * @li data MUST a have even number of points AND it must be transposed
 * @li result MUST have the same size as data points that should be evaluated
 */
class OperationMultipleEvalIterativeX86SimdLinear : public sg::parallel::OperationMultipleEvalVectorized
{
public:
	/**
	 * Constructor of OperationMultipleEvalIterativeSPX86Simd
	 *
	 * Within the constructor sg::base::DataMatrixSP Level and sg::base::DataMatrixSP Index are set up.
	 * If the grid changes during your calculations and you don't want to create
	 * a new instance of this class, you have to call rebuildLevelAndIndex before
	 * doing any further mult or multTranspose calls.
	 *
	 * @param storage Pointer to the grid's gridstorage object
	 * @param dataset dataset that should be evaluated
	 */
    OperationMultipleEvalIterativeX86SimdLinear(sg::base::GridStorage* storage, sg::base::DataMatrix* dataset,
                                                int storageFrom, int storageTo, int datasetFrom, int datasetTo);

	/**
	 * Destructor
	 */
	virtual ~OperationMultipleEvalIterativeX86SimdLinear();

	virtual double multVectorized(sg::base::DataVector& alpha, sg::base::DataVector& result);

	virtual double multTransposeVectorized(sg::base::DataVector& source, sg::base::DataVector& result);

	virtual void rebuildLevelAndIndex();

protected:
	/// Pointer to the grid's GridStorage object
	sg::base::GridStorage* storage;
	/// Timer object to handle time measurements
	sg::base::SGppStopwatch* myTimer;

    int m_storageFrom;
    int m_storageTo;
    int m_datasetFrom;
    int m_datasetTo;
private:
    void adaptDatasetBoundaries();
};

}
}

#endif /* OPERATIONMULTIPLEEVALITERATIVEX86SIMDLINEAR_HPP */
