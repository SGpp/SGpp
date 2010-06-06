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
#include "data/DataMatrix.hpp"
#include "operation/pde/OperationODESolverSystem.hpp"

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
	/// Pointer to the mus
	DataVector* mus;
	/// Pointer to the sigmas
	DataVector* sigmas;
	/// Pointer to the rhos;
	DataVector* rhos;
	/// Pointer to the coefficients of operation Delta
	DataVector* deltaCoef;
	/// Pointer to the coefficients ot operation Gamma
	DataMatrix* gammaCoef;

	void applyLOperatorInner(DataVector& alpha, DataVector& result);

	void applyLOperatorComplete(DataVector& alpha, DataVector& result);

	void applyMassMatrixInner(DataVector& alpha, DataVector& result);

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
	 * @param MPIRank indicates the MPI-Rank of this instance, 0 indicates the master rank
	 */
	BlackScholesODESolverSystem(Grid& SparseGrid, DataVector& alpha, DataVector& mu, DataVector& sigma, DataVector& rho, double r, double TimestepSize, std::string OperationMode = "ExEul", bool bLogTransform = false, size_t MPIRank = 0);

	/**
	 * Std-Destructor
	 */
	virtual ~BlackScholesODESolverSystem();

	void finishTimestep();

	void startTimestep();
};

}

#endif /* BLACKSCHOLESODESOLVERSYSTEM_HPP */
