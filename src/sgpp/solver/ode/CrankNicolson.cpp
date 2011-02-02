/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "grid/common/DirichletUpdateVector.hpp"
#include "solver/ode/CrankNicolson.hpp"

namespace sg
{

CrankNicolson::CrankNicolson(size_t nTimesteps, double timestepSize, ScreenOutput* screen) : ODESolver(nTimesteps, timestepSize), myScreen(screen)
{
	this->residuum = 0.0;
}

CrankNicolson::~CrankNicolson()
{
}

void CrankNicolson::solve(SLESolver& LinearSystemSolver, OperationParabolicPDESolverSystem& System, bool bIdentifyLastStep, bool verbose)
{
	size_t allIter = 0;
    DataVector* rhs = NULL;

	for (size_t i = 0; i < this->nMaxIterations; i++)
	{
		// generate right hand side
		rhs = System.generateRHS();

		// solve the system of the current timestep
		LinearSystemSolver.solve(System, *System.getGridCoefficientsForCG(), *rhs, true, false, -1.0);
	    allIter += LinearSystemSolver.getNumberIterations();

	    if (verbose == true)
	    {
	    	if (myScreen == NULL)
	    	{
	    		std::cout << "Final residuum " << LinearSystemSolver.residuum << "; with " << LinearSystemSolver.getNumberIterations() << " Iterations (Total Iter.: " << allIter << ")" << std::endl;
	    	}
	    }
	    if (myScreen != NULL)
    	{
    		std::stringstream soutput;
    		soutput << "Final residuum " << LinearSystemSolver.residuum << "; with " << LinearSystemSolver.getNumberIterations() << " Iterations (Total Iter.: " << allIter << ")";

    		if (i < this->nMaxIterations-1)
    		{
    			myScreen->update((size_t)(((double)(i+1)*100.0)/((double)this->nMaxIterations)), soutput.str());
    		}
    		else
    		{
    			myScreen->update(100, soutput.str());
    		}
    	}

	    // Do some adjustments on the boundaries if needed, copy values back
	    if (bIdentifyLastStep == false)
	    {
			System.finishTimestep(false);
	    }
	    else
	    {
			if (i < (this->nMaxIterations-1))
			{
				System.finishTimestep(false);
			}
			else
			{
				System.finishTimestep(true);
			}
	    }
	}

	// write some empty lines to console
    if (myScreen != NULL)
	{
    	myScreen->writeEmptyLines(2);
	}

    this->nIterations = allIter;
}

}
