/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef BLACKSCHOLESODESOLVERSYSTEM_HPP
#define BLACKSCHOLESODESOLVERSYSTEM_HPP

#include "grid/Grid.hpp"
#include "data/DataVector.hpp"
#include "operation/common/OperationODESolverSystem.hpp"
#include "grid/common/DirichletUpdateVector.hpp"
#include "grid/common/DirichletGridConverter.hpp"

namespace sg
{

/**
 * This class implements the ODESolverSystem for the BlackScholes
 * Equation.
 */
class BlackScholesODESolverSystem : public OperationODESolverSystem
{
protected:
	/// the riskfree interest rate
	double r;
	/// the delta Operation, on boundary grid
	OperationMatrix* OpDeltaBound;
	/// the Gamma Operation, on boundary grid
	OperationMatrix* OpGammaBound;
	/// the LTwoDotProduct Operation (Mass Matrix), on boundary grid
	OperationMatrix* OpLTwoBound;
	/// the delta Operation, on Inner grid
	OperationMatrix* OpDeltaInner;
	/// the Gamma Operation, on Inner grid
	OperationMatrix* OpGammaInner;
	/// the LTwoDotProduct Operation (Mass Matrix), on Inner grid
	OperationMatrix* OpLTwoInner;
	/// Pointer to the alphas (ansatzfunctions' coefficients)
	DataVector* alpha_complete;
	/// Pointer to the alphas (ansatzfunctions' coefficients; inner points only)
	DataVector* alpha_inner;
	/// Pointer to the mus
	DataVector* mus;
	/// Pointer to the sigmas
	DataVector* sigmas;
	/// Pointer to the rhos;
	DataVector* rhos;
	/// Pointer to the coefficients of operation Delta
	DataVector* deltaCoef;
	/// Pointer to the coefficients ot operation Gamma
	DataVector* gammaCoef;
	/// Routine to modify the boundaries/inner points of the grid
	DirichletUpdateVector* BoundaryUpdate;
	/// Class that allows a simple conversion between a grid with and a without boundary points
	DirichletGridConverter* GridConverter;
	/// DateVector to store the right hand side
	DataVector* rhs;

	/**
	 *  specifies in which solver this matrix is used, valid values are:
	 *  ExEul for explicit Euler
	 *  ImEul for implicit Euler
	 *  CrNic for Crank Nicolson solver
	 */
	std::string tOperationMode;
	/// the size of one timestep used in the ODE Solver
	double TimestepSize;
	/// Pointer to the boundary grid object
	Grid* BoundGrid;
	/// Pointer to the inner grid object
	Grid* InnerGrid;

	/**
	 * Do Matrix mutlitplication with the Black Scholes Systemmatrix, inner grid points only
	 *
	 * @param alpha the coefficients of the sparse grid's ansatzfunctions (inner grid points)
	 * @param return reference to the DataVector into which the result is written (inner grid points)
	 */
	void applyLOperatorInner(DataVector& alpha, DataVector& result);

	/**
	 * Do Matrix mutlitplication with the Black Scholes Systemmatrix, complete boundary grid
	 *
	 * @param alpha the coefficients of the sparse grid's ansatzfunctions
	 * @param return reference to the DataVector into which the result is written
	 */
	void applyLOperatorComplete(DataVector& alpha, DataVector& result);

	/**
	 * Do Matrix mutlitplication with left hand side mass matrix, inner grid points only
	 *
	 * @param alpha the coefficients of the sparse grid's ansatzfunctions (inner grid points)
	 * @param return reference to the DataVector into which the result is written (inner grid points)
	 */
	void applyMassMatrixInner(DataVector& alpha, DataVector& result);

	/**
	 * Do Matrix mutlitplication with left hand side mass matrix
	 *
	 * @param alpha the coefficients of the sparse grid's ansatzfunctions
	 * @param return reference to the DataVector into which the result is written
	 */
	void applyMassMatrixComplete(DataVector& alpha, DataVector& result);

	/**
	 * Build the coefficients for the Gamma Operation, which
	 * are the assets' covariance matrix multiplied by 0.5
	 *
	 * this routine handles also the symmtrie of the
	 * gamma operation
	 */
	void buildGammaCoefficients();

	/**
	 * Build the coefficients for the combined Delta Operation
	 */
	void buildDeltaCoefficients();

	/**
	 * Implements some start jobs of every timestep, e.g.discounting boundaries
	 */
	void startTimestep();

public:
	/**
	 * Std-Constructor
	 *
	 * @param SparseGrid reference to the sparse grid
	 * @param alpha the ansatzfunctions' coefficients
	 * @param mu reference to the mus
	 * @param sigma reference to the sigmas
	 * @param rho reference to the rhos
	 * @param r the riskfree interest rate
	 * @param TimestepSize the size of one timestep used in the ODE Solver
	 * @param OperationMode specifies in which solver this matrix is used, valid values are: ExEul for explicit Euler,
	 *  							ImEul for implicit Euler, CrNic for Crank Nicolson solver
	 */
	BlackScholesODESolverSystem(Grid& SparseGrid, DataVector& alpha, DataVector& mu, DataVector& sigma, DataVector& rho, double r, double TimestepSize, std::string OperationMode = "ExEul");

	/**
	 * Std-Destructor
	 */
	virtual ~BlackScholesODESolverSystem();

	virtual void mult(DataVector& alpha, DataVector& result);

	virtual DataVector* generateRHS();

	virtual void finishTimestep();

	virtual Grid* getGrid();

	virtual DataVector* getGridCoefficientsForCG();

	virtual DataVector* getGridCoefficients();
};

}

#endif /* BLACKSCHOLESODESOLVERSYSTEM_HPP */
