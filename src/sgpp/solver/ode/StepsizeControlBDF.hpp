/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/

#ifndef STEPSIZECONTROLBDF_HPP
#define STEPSIZECONTROLBDF_HPP

#include "application/common/ScreenOutput.hpp"
#include "solver/ODESolver.hpp"
using namespace sg::pde;
using namespace sg::base;

namespace sg
{
namespace solver
{

/**
 * This class implements a step size control using the midpoint method and BDF2
 * for solving ordinary partial equations
 *
 * For solving the system of linear equations the
 * already implemented CG-method is used
 *
 * @version $HEAD$
 */
class StepsizeControlBDF : public ODESolver
{
private:
	/// Pointer to ScreenOutput object
	ScreenOutput* myScreen;

	/// epsilon for the step size control
	double myEps;

public:
	/**
	 * Std-Constructer
	 *
	 * @param nTimesteps number of maximum executed iterations
	 * @param timestepSize the size of one timestep
	 * @param eps the epsilon for the step size control
	 * @param screen possible pointer to a ScreenOutput object
	 */
	StepsizeControlBDF(size_t nTimesteps, double timestepSize, double eps, ScreenOutput* screen = NULL);

	/**
	 * Std-Destructor
	 */
	virtual ~StepsizeControlBDF();

	virtual void solve(SLESolver& LinearSystemSolver, OperationParabolicPDESolverSystem& System, bool bIdentifyLastStep = false, bool verbose = false);
};

}
}

#endif /* STEPSIZECONTROLBDF_HPP */
