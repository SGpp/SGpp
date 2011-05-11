/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef BLACKSCHOLESPARABOLICPDESOLVERSYSTEM_HPP
#define BLACKSCHOLESPARABOLICPDESOLVERSYSTEM_HPP

#include "grid/Grid.hpp"
#include "data/DataVector.hpp"
#include "data/DataMatrix.hpp"
#include "grid/common/DirichletUpdateVector.hpp"
#include "operation/pde/OperationParabolicPDESolverSystemNeumann.hpp"
using namespace sg::pde;
using namespace sg::base;

namespace sg
{
namespace finance
{
/**
 * This class implements the ParabolicPDESolverSystem for the BlackScholes
 * Equation.
 */
class BlackScholesParabolicPDESolverSystem : public OperationParabolicPDESolverSystemNeumann
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
	/// Pointer to the mus
	DataVector* mus;
	/// Pointer to the sigmas
	DataVector* sigmas;
	/// Pointer to the rhos;
	DataMatrix* rhos;
	/// Pointer to the coefficients of operation Delta
	DataVector* deltaCoef;
	/// Pointer to the coefficients ot operation Gamma
	DataMatrix* gammaCoef;
	/// use coarsening between timesteps in order to reduce gridsize
	bool useCoarsen;
	/// adaptive mode during solving Black Scholes Equation: coarsen, refine, coarsenNrefine
	std::string adaptSolveMode;
	/// number of points the are coarsened in each coarsening-step !CURRENTLY UNUSED PARAMETER!
	int numCoarsenPoints;
	/// Threshold used to decide if a grid point should be deleted
	double coarsenThreshold;
	/// Threshold used to decide if a grid point should be refined
	double refineThreshold;
	/// refine mode during solving Black Scholes Equation: classic or maxLevel
	std::string refineMode;
	/// maxLevel max. Level of refinement
	size_t refineMaxLevel;
	/// the algorithmic dimensions used in this system
	std::vector<size_t> BSalgoDims;
	/// Routine to modify the boundaries/inner points of the grid
	DirichletUpdateVector* BoundaryUpdate;

	virtual void applyLOperator(DataVector& alpha, DataVector& result);

	virtual void applyMassMatrix(DataVector& alpha, DataVector& result);

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
	 * Build the coefficients for the Gamma Operation, which
	 * are the assets' covariance matrix multiplied by 0.5
	 *
	 * this routine handles also the symmtrie of the
	 * gamma operation
	 *
	 * This function builds the coefficients for the Log Transformed Black Scholes Equation
	 */
	void buildGammaCoefficientsLogTransform();

	/**
	 * Build the coefficients for the combined Delta Operation
	 *
	 * This function builds the coefficients for the Log Transformed Black Scholes Equation
	 */
	void buildDeltaCoefficientsLogTransform();

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
	BlackScholesParabolicPDESolverSystem(Grid& SparseGrid, DataVector& alpha, DataVector& mu, DataVector& sigma,
			DataMatrix& rho, double r, double TimestepSize, std::string OperationMode = "ExEul",
			bool bLogTransform = false, bool useCoarsen = false, double coarsenThreshold = 0.0, std::string adaptSolveMode ="none",
			int numCoarsenPoints = -1, double refineThreshold = 0.0, std::string refineMode = "classic", size_t refineMaxLevel = 0);

	/**
	 * Std-Destructor
	 */
	virtual ~BlackScholesParabolicPDESolverSystem();

	virtual void finishTimestep(bool isLastTimestep = false);

	virtual void startTimestep();
};
}
}

#endif /* BLACKSCHOLESPARABOLICPDESOLVERSYSTEM_HPP */
