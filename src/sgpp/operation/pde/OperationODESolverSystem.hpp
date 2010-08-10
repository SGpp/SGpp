/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef OPERATIONODESOLVERSYSTEM_HPP
#define OPERATIONODESOLVERSYSTEM_HPP

#include "grid/Grid.hpp"
#include "operation/common/OperationMatrix.hpp"
#include "data/DataVector.hpp"
#include "grid/common/DirichletUpdateVector.hpp"
#include "grid/common/DirichletGridConverter.hpp"

namespace sg
{

/**
 * Defines a System that is used to solve parabolic partial
 * differential equations. So an instance of this class has to pass to
 * any ODE Solver used in SGpp.
 *
 * \f$\A \dot{u} = L \vec{u}\f$
 *
 * A: mass matrix
 * L: space discretization (L-Operator)
 *
 * This class defines an elliptic problem in every timestep which is solved
 * using an iterative SLE solver, that solving step is integrated in the
 * ODE Solver.
 */
class OperationODESolverSystem : public OperationMatrix
{
protected:
	/// Pointer to the alphas (ansatzfunctions' coefficients)
	DataVector* alpha_complete;
	/// Pointer to the alphas (ansatzfunctions' coefficients; inner points only)
	DataVector* alpha_inner;



	DataVector* alpha_complete_old;
	DataVector* alpha_complete_tmp;

	/**
	 *  specifies in which solver this matrix is used, valid values are:
	 *  ExEul for explicit Euler
	 *  ImEul for implicit Euler
	 *  CrNic for Crank Nicolson solver
	 */
	std::string tOperationMode;
	/// the size of one timestep used in the ODE Solver
	double TimestepSize;

	double TimestepSize_old;

	/// Routine to modify the boundaries/inner points of the grid
	DirichletUpdateVector* BoundaryUpdate;
	/// Class that allows a simple conversion between a grid with and a without boundary points
	DirichletGridConverter* GridConverter;
	/// DateVector to store the right hand side
	DataVector* rhs;
	/// Pointer to the grid object
	Grid* BoundGrid;
	/// Pointer to the inner grid object
	Grid* InnerGrid;
	/// Stores number of average gridpoints, inner grid
	size_t numSumGridpointsInner;
	/// Stores number of average gridpoints, complete grid
	size_t numSumGridpointsComplete;

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
	OperationODESolverSystem();

	/**
	 * Destructor
	 */
	virtual ~OperationODESolverSystem();

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
	DataVector* generateRHS();

	/**
	 * performs some action that might be needed after a timestep has be finished in the ODE
	 * Solver, e.g. some boundary adjustments.
	 *
	 * @param isLastTimestep denotes of this is the clean up for the last time of solving the ODE
	 */
	virtual void finishTimestep(bool isLastTimestep = false);

	/**
	 * Implements some start jobs of every timestep, e.g.discounting boundaries
	 */
	virtual void startTimestep();

	/**
	 * get the pointer to the underlying grid object
	 *
	 * @return returns a pointer to the underlying grid object
	 */
	Grid* getGrid();

	/**
	 * gets a pointer to the sparse grids coefficients used in the CG method to solve
	 * one timestep. This is useful because (direchlet) boundaries can be skipped when
	 * solving the system.
	 *
	 * @return alpha vector for CG method
	 */
	DataVector* getGridCoefficientsForCG();

	/**
	 * gets a pointer to the sparse grids coefficients with evtl. boundaries
	 *
	 * @return alpha vector of complete grid
	 */
	DataVector* getGridCoefficients();

	/**
	 * defines the used ODE Solver for this instance, this is important because
	 * the implementation of mult and generateRHS depends on the used
	 * ODE solver
	 *
	 * @param ode the used ODESolver: ExEul, ImEul or CrNic
	 */
	void setODESolver(std::string ode);

	/**
	 * Returns the specified ODE solver for this instance
	 *
	 * @return the ODE solver: ExEul, ImEul or CrNic
	 */
	std::string getODESolver();

	/**
	 * Returns the number of average grid points for the complete grid
	 *
	 * @return the number of average grid points for the complete grid
	 */
	size_t getSumGridPointsComplete();

	/**
	 * Returns the number of average grid points for the inner grid
	 *
	 * @return the number of average grid points for the inner grid
	 */
	size_t getSumGridPointsInner();


	void setTimestepSize(double newTimestepSize);
	void abortTimestep();
	void saveAlpha();
};

}

#endif /* OPERATIONODESOLVERMATRIX_HPP */
