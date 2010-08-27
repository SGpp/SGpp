/******************************************************************************
* Copyright (C) 2010 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef BLACKSCHOLESODESOLVERSYSTEMPARALLELOMP_HPP
#define BLACKSCHOLESODESOLVERSYSTEMPARALLELOMP_HPP

#include "algorithm/pde/BlackScholesODESolverSystem.hpp"

namespace sg
{

/**
 * This class implements the ODESolverSystem for the BlackScholes
 * Equation.
 *
 * It's derived from the existing BlackScholesODESolverSystem but uses
 * the OMP task concept to enable further parallelization possibilities
 * in the calculation of the space-discretization operator (L)
 */
class BlackScholesODESolverSystemParallelOMP : public BlackScholesODESolverSystem
{
protected:
	virtual void applyLOperatorInner(DataVector& alpha, DataVector& result);

	virtual void applyLOperatorComplete(DataVector& alpha, DataVector& result);

	virtual void applyMassMatrixInner(DataVector& alpha, DataVector& result);

	virtual void applyMassMatrixComplete(DataVector& alpha, DataVector& result);

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
	 * @param bLogTransform indicates that this system belongs to a log-transformed Black Scholes Equation
	 * @param useCoarsen specifies if the grid should be coarsened between timesteps
	 * @param coarsenThreshold Threshold to decide, if a grid point should be deleted
	 * @param adaptSolveMode adaptive mode during solving: coarsen, refine, coarsenNrefine
	 * @param numCoarsenPoints number of point that should be coarsened in one coarsening step !CURRENTLY UNUSED PARAMETER!
	 * @param refineThreshold Threshold to decide, if a grid point should be refined
	 * @param refineMode refineMode during solving Black Scholes Equation: classic or maxLevel
	 * @param refineMaxLevel max. level for refinement during solving
	 */
	BlackScholesODESolverSystemParallelOMP(Grid& SparseGrid, DataVector& alpha, DataVector& mu, DataVector& sigma,
			DataMatrix& rho, double r, double TimestepSize, std::string OperationMode = "ExEul",
			bool bLogTransform = false, bool useCoarsen = false, double coarsenThreshold = 0.0, std::string adaptSolveMode = "none",
			int numCoarsenPoints = -1, double refineThreshold = 0.0, std::string refineMode = "classic", size_t refineMaxLevel = 0);

	/**
	 * Std-Destructor
	 */
	virtual ~BlackScholesODESolverSystemParallelOMP();

	/**
	 * Multiplicates a vector with the matrix, parallel
	 *
	 * @param alpha DataVector that contains the ansatzfunctions' coefficients
	 * @param result DataVector into which the result of the space discretization operation is stored
	 */
	virtual void mult(DataVector& alpha, DataVector& result);

	/**
	 * generates the right hand side of the system, parallel
	 *
	 * @return returns the rhs
	 */
	virtual DataVector* generateRHS();
};

}

#endif /* BLACKSCHOLESODESOLVERSYSTEMPARALLELOMP_HPP */
