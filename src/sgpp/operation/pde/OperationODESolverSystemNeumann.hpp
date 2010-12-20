/******************************************************************************
* Copyright (C) 2010 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef OPERATIONODESOLVERSYSTEMNEUMANN_HPP
#define OPERATIONODESOLVERSYSTEMNEUMANN_HPP

#include "operation/pde/OperationODESolverSystem.hpp"

namespace sg
{

/**
 * Defines a System that is used to solve parabolic partial
 * differential equations. So an instance of this class has to pass to
 * any ODE Solver used in SGpp.
 *
 * \f$A \dot{u} = L \vec{u}\f$
 *
 * A: mass matrix
 * L: space discretization (L-Operator)
 *
 * This class defines an elliptic problem in every timestep which is solved
 * using an iterative SLE solver, that solving step is integrated in the
 * ODE Solver.
 */
class OperationODESolverSystemNeumann : public OperationODESolverSystem
{
protected:
	/**
	 * applies the PDE's mass matrix, on complete grid - with boundaries
	 *
	 * @param alpha the coefficients of the sparse grid's ansatzfunctions
	 * @param return reference to the DataVector into which the result is written
	 */
	virtual void applyMassMatrix(DataVector& alpha, DataVector& result) = 0;

	/**
	 * applies the PDE's system matrix, on complete grid - with boundaries
	 *
	 * @param alpha the coefficients of the sparse grid's ansatzfunctions
	 * @param return reference to the DataVector into which the result is written
	 */
	virtual void applyLOperator(DataVector& alpha, DataVector& result) = 0;

public:
	/**
	 * Constructor
	 */
	OperationODESolverSystemNeumann();

	/**
	 * Destructor
	 */
	virtual ~OperationODESolverSystemNeumann();

	virtual void mult(DataVector& alpha, DataVector& result);

	virtual DataVector* generateRHS();

	virtual DataVector* getGridCoefficientsForCG();
};

}

#endif /* OPERATIONODESOLVERMATRIXNEUMANN_HPP */
