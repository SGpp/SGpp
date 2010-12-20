/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Chao qi (qic@in.tum.de)

#ifndef MODIFIEDBLACKSCHOLESODESOLVERSYSTEM_HPP
#define MODIFIEDBLACKSCHOLESODESOLVERSYSTEM_HPP

#include "grid/Grid.hpp"
#include "data/DataVector.hpp"
#include "data/DataMatrix.hpp"
#include "algorithm/pde/BlackScholesODESolverSystem.hpp"

namespace sg
{

/**
 * This class implements the Modified ODESolverSystem for the BlackScholes
 * Equation just use for combination of BlackScholes and HullWhite.
 */
class ModifiedBlackScholesODESolverSystem  : public BlackScholesODESolverSystem
{
protected:

	OperationMatrix* OpFBound;

	virtual void applyLOperator(DataVector& alpha, DataVector& result);

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
	 * @param refineMaxLevel max. level of refinement during solving
	 */
   	ModifiedBlackScholesODESolverSystem(Grid& SparseGrid, DataVector& alpha, DataVector& mu,
			DataVector& sigma, DataMatrix& rho, double r, double TimestepSize, std::string OperationMode,
			bool bLogTransform, bool useCoarsen, double coarsenThreshold, std::string adaptSolveMode,
			int numCoarsenPoints, double refineThreshold, std::string refineMode, size_t refineMaxLevel);
	/**
	 * Multiplies the corresponding r coordinates with the whole grid value
	 *
	 * @param updateVector the vector that should be updated
	 */
   	void multiplyrBSHW(DataVector& updateVector);

   	/**
	 * Std-Destructor
	 */
	virtual ~ModifiedBlackScholesODESolverSystem();

	virtual void finishTimestep(bool isLastTimestep = false);

	virtual void startTimestep();
};

}

#endif /* MODIFIEDBLACKSCHOLESODESOLVERSYSTEM_HPP */
