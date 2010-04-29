/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#ifndef CRANKNICOLSON_HPP
#define CRANKNICOLSON_HPP

#include "application/common/ScreenOutput.hpp"
#include "solver/ODESolver.hpp"

namespace sg
{

/**
 * This class implements the Crank-Nicolson method
 * for solving ordinary partial equations
 *
 * For solving the system of linear equations the
 * already implemented CG-method is used
 *
 * @version $HEAD$
 */
class CrankNicolson : public ODESolver
{
private:
	/// the number of CG maximum CG iterations
	size_t maxCGIterations;
	/// the CG's epsilon
	double epsilonCG;
	/// Pointer to ScreenOutput object
	ScreenOutput* myScreen;

public:
	/**
	 * Std-Constructer
	 *
	 * @param nTimesteps number of maximum executed iterations
	 * @param timestepSize the size of one timestep
	 * @param iMaxCG maximum number of CG steps
	 * @param epsilonCG the epsilon used in CG
	 * @param screen possible pointer to a ScreenOutput object
	 */
	CrankNicolson(size_t nTimesteps, double timestepSize, size_t iMaxCG, double epsilonCG, ScreenOutput* screen = NULL);

	/**
	 * Std-Destructor
	 */
	virtual ~CrankNicolson();

	virtual void solve(OperationODESolverMatrix& SystemMatrix, DataVector& alpha, bool verbose = false);
};

}

#endif /* CRANKNICOLSON_HPP */
