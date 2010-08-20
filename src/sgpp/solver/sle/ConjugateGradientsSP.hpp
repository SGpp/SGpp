/******************************************************************************
* Copyright (C) 2010 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef CONJUGATEGRADIENTSSP_HPP
#define CONJUGATEGRADIENTSSP_HPP

#include "operation/common/OperationMatrixSP.hpp"
#include "data/DataVectorSP.hpp"

#include <iostream>

namespace sg
{

class ConjugateGradientsSP
{
private:
	/// Number of Iterations needed for the solve
	size_t nIterations;
	/// Number of maximum iterations for cg
	size_t nMaxIterations;
	/// residuum
	float residuum;
	/// epsilon needed in the, e.g. final error in the iterative solver, or a timestep
	float myEpsilon;

public:
	/**
	 * Std-Constructor
	 */
	ConjugateGradientsSP(size_t imax, float epsilon);

	/**
	 * Std-Destructor
	 */
	virtual ~ConjugateGradientsSP();

	/**
	 * function that defines a solve method for an iterative solver. In contrast to the normal solve routine
	 * this method operates on sinlge precision data.
	 *
	 * @param SystemMatrix reference to an OperationMatrix Object that implements the matrix vector multiplication
	 * @param alpha the sparse grid's coefficients which have to be determined
	 * @param b the right hand side of the system of linear equations
	 * @param reuse identifies if the alphas, stored in alpha at calling time, should be reused
	 * @param verbose prints information during execution of the solver
	 * @param max_threshold additional abort criteria for solver
	 */
	void solve(OperationMatrixSP& SystemMatrix, DataVectorSP& alpha, DataVectorSP& b, bool reuse = false, bool verbose = false, float max_threshold = -1.0);

	/**
	 * function that returns the number of needed solve steps
	 *
	 * @return the number of needed solve steps of the sovler
	 */
	size_t getNumberIterations();

	/**
	 * function the returns the residuum (current or final), error of the solver
	 *
	 * @return the residuum
	 */
	float getResiduum();

	/**
	 * resets the number of maximum iterations
	 *
	 * @param nIterations the new number of maximum iterations
	 */
	void setMaxIterations(size_t nIterations);

	/**
	 * resets the epsilon, that is used in the SGSolver
	 *
	 * @param eps the new value of epsilon
	 */
	void setEpsilon(float eps);

	/**
	 * gets the the epsilon, that is used in the SGSolver
	 *
	 * @return the epsilon, used in the solver
	 */
	float getEpsilon();
};

}

#endif /* CONJUGATEGRADIENTSSP_HPP */
