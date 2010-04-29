/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de), Stefanie Schraufstetter (schraufs@in.tum.de)

#ifndef LOGBLACKSCHOLESODESOLVERSYSTEM_HPP
#define LOGBLACKSCHOLESODESOLVERSYSTEM_HPP

#include "grid/Grid.hpp"
#include "data/DataVector.hpp"
#include "algorithm/pde/BlackScholesODESolverSystem.hpp"
#include "grid/common/DirichletUpdateVector.hpp"


namespace sg
{

/**
 * implements the ODESolverSystem for the log-transformed Black-Scholes equation
 */
class LogBlackScholesODESolverSystem : public BlackScholesODESolverSystem
{
private:
	/**
	 * Build the coefficients for the Gamma Operation
	 * (term for second derivative in PDE), which are the
	 * assets' covariance matrix multiplied by 0.5
	 *
	 * this routine handles also the symmetry of the
	 * gamma operation
	 */
	void buildGammaCoefficients();

	/**
	 * Build the coefficients for the combined Delta Operation
	 * (term for first derivative in PDE), which are the drift factors
	 */
	void buildDeltaCoefficients();


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
	LogBlackScholesODESolverSystem(Grid& SparseGrid, DataVector& alpha, DataVector& mu, DataVector& sigma, DataVector& rho, double r, double TimestepSize, std::string OperationMode = "ExEul");
};

}

#endif /* LOGBLACKSCHOLESODESOLVERSYSTEM_HPP */
