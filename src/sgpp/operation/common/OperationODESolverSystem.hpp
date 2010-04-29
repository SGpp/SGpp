/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Jörg Blank (blankj@in.tum.de), Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef OPERATIONODESOLVERSYSTEM_HPP
#define OPERATIONODESOLVERSYSTEM_HPP

#include "grid/Grid.hpp"
#include "operation/common/OperationMatrix.hpp"
#include "data/DataVector.hpp"

namespace sg
{

/**
 * Defines a System that is used to solve parabolic partial
 * differential equations. So an instance of this class has to pass to
 * any ODE Solver used in SGpp.
 *
 * This class defines an elliptic problem in every timestep which is solved
 * using an iterative SLE solver, that solving step is integrated in the
 * ODE Solver.
 */
class OperationODESolverSystem : public OperationMatrix
{
public:
	/**
	 * Constructor
	 */
	OperationODESolverSystem() {}

	/**
	 * Destructor
	 */
	virtual ~OperationODESolverSystem() {}

	/**
	 * Multiplicates a vector with the matrix
	 *
	 * @param alpha DataVector that contains the ansatzfunctions' coefficients
	 * @param result DataVector into which the result of the Laplace operation is stored
	 */
	virtual void mult(DataVector& alpha, DataVector& result) = 0;

	/**
	 * generates the right hand side of the system
	 *
	 * @return returns the rhs
	 */
	virtual DataVector* generateRHS() = 0;

	/**
	 * performs some action that might be needed after a timestep has be finished in the ODE
	 * Solver, e.g. some boundary adjustments.
	 */
	virtual void finishTimestep() = 0;

	/**
	 * get the pointer to the underlying grid object
	 *
	 * @return returns a pointer to the underlying grid object
	 */
	virtual Grid* getGrid() = 0;

	/**
	 * gets a pointer to the sparse grids coefficients used in the CG method to solve
	 * one timestep. This is useful because (direchlet) boundaries can be skipped when
	 * solving the system.
	 *
	 * @return alpha vector for CG method
	 */
	virtual DataVector* getGridCoefficientsForCG() = 0;

	/**
	 * gets a pointer to the sparse grids coefficients with evtl. boundaries
	 *
	 * @return alpha vector of complete grid
	 */
	virtual DataVector* getGridCoefficients() = 0;
};

}

#endif /* OPERATIONODESOLVERMATRIX_HPP */
