/******************************************************************************
* Copyright (C) 2010 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef OPERATIONPARABOLICPDESOLVERSYSTEMDIRICHLET_HPP
#define OPERATIONPARABOLICPDESOLVERSYSTEMDIRICHLET_HPP

#include "operation/pde/OperationParabolicPDESolverSystem.hpp"
#include "grid/common/DirichletUpdateVector.hpp"
#include "grid/common/DirichletGridConverter.hpp"
using namespace sg::base;

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
 *
 * This class is a specialized version of OperationParabolicPDESolverSystem which
 * exploit Dirichlet boundary conditions. Hence there are no degrees of freedom
 * on on the boundaries the iterative solver (CG or BiCGSTAB) has only to take
 * inner grid points into account.
 */
class OperationParabolicPDESolverSystemDirichlet : public OperationParabolicPDESolverSystem
{
protected:
	/// Pointer to the alphas (ansatzfunctions' coefficients; inner points only)
	DataVector* alpha_inner;
	/// Routine to modify the boundaries/inner points of the grid
	DirichletUpdateVector* BoundaryUpdate;
	/// Class that allows a simple conversion between a grid with and a without boundary points
	DirichletGridConverter* GridConverter;
	/// Pointer to the inner grid object
	Grid* InnerGrid;

	/**
	 * applies the PDE's mass matrix, on complete grid - with boundaries
	 *
	 * @param alpha the coefficients of the sparse grid's ansatzfunctions
	 * @param return reference to the DataVector into which the result is written
	 */
	virtual void applyMassMatrixComplete(DataVector& alpha, DataVector& result) = 0;

	/**
	 * applies the PDE's system matrix, on complete grid - with boundaries
	 *
	 * @param alpha the coefficients of the sparse grid's ansatzfunctions
	 * @param return reference to the DataVector into which the result is written
	 */
	virtual void applyLOperatorComplete(DataVector& alpha, DataVector& result) = 0;

	/**
	 * applies the PDE's mass matrix, on inner grid only
	 *
	 * @param alpha the coefficients of the sparse grid's ansatzfunctions
	 * @param return reference to the DataVector into which the result is written
	 */
	virtual void applyMassMatrixInner(DataVector& alpha, DataVector& result) = 0;

	/**
	 * applies the PDE's system matrix, on inner grid only
	 *
	 * @param alpha the coefficients of the sparse grid's ansatzfunctions
	 * @param return reference to the DataVector into which the result is written
	 */
	virtual void applyLOperatorInner(DataVector& alpha, DataVector& result) = 0;

public:
	/**
	 * Constructor
	 */
	OperationParabolicPDESolverSystemDirichlet();

	/**
	 * Destructor
	 */
	virtual ~OperationParabolicPDESolverSystemDirichlet();

	/**
	 * Multiplicates a vector with the matrix
	 *
	 * @param alpha DataVector that contains the ansatzfunctions' coefficients
	 * @param result DataVector into which the result of the space discretization operation is stored
	 */
	virtual void mult(DataVector& alpha, DataVector& result);

	/**
	 * generates the right hand side of the system
	 *
	 * @return returns the rhs
	 */
	virtual DataVector* generateRHS();

	virtual DataVector* getGridCoefficientsForCG();
};

}

#endif /* OPERATIONPARABOLICPDESOLVERSYSTEMDIRICHLET_HPP */
