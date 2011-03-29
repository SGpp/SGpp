/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/

#ifndef ADAMSBASHFORTH_HPP
#define ADAMSBASHFORTH_HPP

#include "application/common/ScreenOutput.hpp"
#include "solver/ODESolver.hpp"
#include <string>
using namespace sg::pde;
using namespace sg::base;

namespace sg
{
namespace solver
{

/**
 * This class implements the explicit Adams-Bashforth method
 * for solving ordinary partial equations
 *
 * @version $HEAD$
 */
class AdamsBashforth : public ODESolver
{
private:
	/// Pointer to ScreenOutput object
	ScreenOutput* myScreen;

public:
	/**
	 * Std-Constructer
	 *
	 * @param imax number of maximum executed iterations
	 * @param timestepSize the size of one timestep
	 * @param screen possible pointer to a ScreenOutput object
	 */
	AdamsBashforth(size_t imax, double timestepSize, ScreenOutput* screen = NULL);

	/**
	 * Std-Destructor
	 */
	virtual ~AdamsBashforth();

	virtual void solve(SLESolver& LinearSystemSolver, OperationParabolicPDESolverSystem& System, bool bIdentifyLastStep = false, bool verbose = false);
};

}
}

#endif /* ADAMSBASHFORTH_HPP */
