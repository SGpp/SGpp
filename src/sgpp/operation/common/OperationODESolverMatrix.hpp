/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Jörg Blank (blankj@in.tum.de), Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef OPERATIONODESOLVERMATRIX_HPP
#define OPERATIONODESOLVERMATRIX_HPP

#include "grid/Grid.hpp"
#include "operation/common/OperationMatrix.hpp"
#include "data/DataVector.hpp"

namespace sg
{

/**
 * Defines a Systemmatrix that is used to solve parabolic partial
 * differential equations. So an instance of this class has to pass to
 * any ODE Solver used in SGpp.
 *
 * This class defines an elliptic problem in every timestep which is solved
 * using an iterative SLE solver, that solving step is integrated in the
 * ODE Solver.
 */
class OperationODESolverMatrix : public OperationMatrix
{
public:
	/**
	 * Constructor
	 */
	OperationODESolverMatrix() {}

	/**
	 * Destructor
	 */
	virtual ~OperationODESolverMatrix() {}

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
	 * @param alpha data that is needed to generate the RHS of the system of linear equations
	 * @param rhs DataVector the contains the RHS
	 */
	virtual void generateRHS(DataVector& alpha, DataVector& rhs) = 0;

	/**
	 * performs some action that might be needed after a timestep has be finished in the ODE
	 * Solver, e.g. some boundary adjustments.
	 *
	 * @param alpha DataVector that contains the ansatzfunctions' coefficients
	 */
	virtual void finishTimestep(DataVector& alpha) = 0;

	/**
	 * performs some action that might be needed before a timestep is executed in the ODE
	 * Solver, e.g. some boundary adjustments.
	 *
	 * @param alpha DataVector that contains the ansatzfunctions' coefficients
	 */
	virtual void startTimestep(DataVector& alpha) = 0;

	/**
	 * get the pointer to the underlying grid object
	 *
	 * @return returns a pointer to the underlying grid object
	 */
	virtual Grid* getGrid() = 0;
};

}

#endif /* OPERATIONODESOLVERMATRIX_HPP */
