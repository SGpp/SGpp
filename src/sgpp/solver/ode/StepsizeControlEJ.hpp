/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/

#ifndef STEPSIZECONTROLEJ_HPP
#define STEPSIZECONTROLEJ_HPP

#include "application/common/ScreenOutput.hpp"
#include "solver/ODESolver.hpp"

namespace sg
{
namespace solver
{

/**
 * This class implements a time step size control using 1D-Diffusion
 * for solving ordinary partial equations
 *
 * For solving the system of linear equations the
 * already implemented CG-method is used
 *
 * @version $HEAD$
 */
class StepsizeControlEJ : public ODESolver
{
private:
	/// Pointer to sg::base::ScreenOutput object
	sg::base::ScreenOutput* myScreen;

	/// epsilon for the stepsize control
	double myEps;

	double mySC;

public:
	/**
	 * Std-Constructer
	 *
	 * @param nTimesteps number of maximum executed iterations
	 * @param timestepSize the size of one timestep
	 * @param eps the epsilon for the stepsize control
	 * @param sc
	 * @param screen possible pointer to a sg::base::ScreenOutput object
	 */
	StepsizeControlEJ(size_t nTimesteps, double timestepSize, double eps, double sc, sg::base::ScreenOutput* screen = NULL);

	/**
	 * Std-Destructor
	 */
	virtual ~StepsizeControlEJ();

	virtual void solve(SLESolver& LinearSystemSolver, sg::pde::OperationParabolicPDESolverSystem& System, bool bIdentifyLastStep = false, bool verbose = false);
};

}
}

#endif /* STEPSIZECONTROLEJ_HPP */
