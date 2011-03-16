/******************************************************************************
* Copyright (C) 2010 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef DMSYSTEMMATRIXSPVECTORIZEDIDENTITY_HPP
#define DMSYSTEMMATRIXSPVECTORIZEDIDENTITY_HPP

#include "data/DataVectorSP.hpp"
#include "grid/Grid.hpp"
#include "operation/datadriven/OperationMultipleEvalVectorizedSP.hpp"
#include "operation/common/OperationMatrixSP.hpp"
#include "tools/common/SGppStopwatch.hpp"

#include <string>
using namespace sg::base;

namespace sg
{

/**
 * Class that implements the virtual class OperationMatrix for the
 * application of classification for the Systemmatrix
 *
 * The Identity matrix is used as regularization operator.
 *
 * For the Operation B's mult and mutlTransposed functions
 * vectorized formulations in SSE, AVX, OpenCL or Intel Array Building Blocks
 * are used.
 *
 * In this class single precision DataVectors are used.
 */
class DMSystemMatrixSPVectorizedIdentity : public OperationMatrixSP
{
private:
	/// the lambda, the regularisation parameter
	float lamb;
	/// OperationB for calculating the data matrix
	OperationMultipleEvalVectorizedSP* B;
	/// Pointer to the data matrix
	DataMatrixSP* data;
	/// Number of original training instances
	size_t numTrainingInstances;
	/// Number of patched and used training instances
	size_t numPatchedTrainingInstances;
	/// vectorization mode, possible values are SSE, AVX, OCL, ArBB
	std::string vecMode;
	/// vector width, class internal variable to enable padding and patching of vectors
	size_t vecWidth;
	// save some timings during computation
	/// time needed for Mult
	double completeTimeMult;
	/// time needed only for the computation of mult, interesting on accelerator boards
	double computeTimeMult;
	/// time needed for Mult transposed
	double completeTimeMultTrans;
	/// time needed only for the computation of mult transposed, interesting on accelerator boards
	double computeTimeMultTrans;
	/// Stopwatch needed to determine the durations of mult and mult transposed
	SGppStopwatch* myTimer;

public:
	/**
	 * Std-Constructor
	 *
	 * @param SparseGrid reference to the sparse grid
	 * @param trainData reference to DataMatrix that contains the training data
	 * @param C the regression functional
	 * @param lambda the lambda, the regression parameter
	 * @param vecMode vectorization mode, possible values are SSE, AVX, OCL, ArBB
	 */
	DMSystemMatrixSPVectorizedIdentity(Grid& SparseGrid, DataMatrixSP& trainData, float lambda, std::string vecMode);

	/**
	 * Std-Destructor
	 */
	virtual ~DMSystemMatrixSPVectorizedIdentity();

	virtual void mult(DataVectorSP& alpha, DataVectorSP& result);

	/**
	 * Generates the right hand side of the classification equation
	 *
	 * @param classes the class information of the training data
	 * @param b reference to the vector that will contain the result of the matrix vector multiplication on the rhs
	 */
	void generateb(DataVectorSP& classes, DataVectorSP& b);

	/**
	 * rebuilds the DataMatrix for Level and Index
	 */
	void rebuildLevelAndIndex();

	/**
	 * resets all timers to 0
	 */
	void resetTimers();

	/**
	 * gets the timer's values by saving them into call by reference values
	 *
	 * @param timeMult variable to store overall time needed for Mult
	 * @param computeMult variable to store compute time needed for Mult
	 * @param timeMultTrans variable to store everall time needed for Mult Transposed
	 * @param computeMultTrans variable to store compute time needed for Mult Transposed
	 */
	void getTimers(double& timeMult, double& computeMult, double& timeMultTrans, double& computeMultTrans);
};

}

#endif /* DMSYSTEMMATRIXSPVECTORIZEDIDENTITY_HPP */
