/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "solver/sle/ConjugateGradients.hpp"
#include "grid/common/DirichletUpdateVector.hpp"
#include "solver/ode/Euler.hpp"
#include "operation/common/OperationEval.hpp"
#include "tools/common/GridPrinter.hpp"
#include "solver/sle/BiCGStab.hpp"
#include "exception/solver_exception.hpp"

#include <iostream>
#include <string>
#include <sstream>

namespace sg
{

Euler::Euler(std::string Mode, size_t imax, double timestepSize, size_t iMaxCG, double epsilonCG, bool generateAnimation, size_t numEvalsAnimation, ScreenOutput* screen) : ODESolver(imax, timestepSize), maxCGIterations(iMaxCG), epsilonCG(epsilonCG), bAnimation(generateAnimation), evalsAnimation(numEvalsAnimation), ExMode(Mode), myScreen(screen)
{
	this->residuum = 0.0;

	if (Mode != "ExEul" && Mode != "ImEul")
	{
		throw new solver_exception("Euler::Euler : An unknown Euler-Mode was specified!");
	}
}

Euler::~Euler()
{
}

void Euler::solve(OperationODESolverMatrix& SystemMatrix, DataVector& alpha, bool verbose)
{
	size_t allIter = 0;
	DataVector rhs(alpha.getSize());
    BiCGStab myCG(this->maxCGIterations, this->epsilonCG);

    // Do some animation creation exception handling
    size_t animationStep = this->nMaxIterations/1500;
    if (animationStep == 0)
    {
    	animationStep = 1;
    }

	for (size_t i = 0; i < this->nMaxIterations; i++)
	{
		rhs.setAll(0.0);


		// generate right hand side
		SystemMatrix.generateRHS(alpha, rhs);


	    // Do some adjustments on the boundaries if needed
		if (this->ExMode == "ImEul")
		{
		    SystemMatrix.startTimestep(alpha);
		}


		// solve the system of the current timestep
	    myCG.solve(SystemMatrix, alpha, rhs, true, false, -1.0);
	    allIter += myCG.getNumberIterations();
	    if (verbose == true)
	    {
	    	if (myScreen == NULL)
	    	{
	    		std::cout << "Final residuum " << myCG.residuum << "; with " << myCG.getNumberIterations() << " Iterations (Total Iter.: " << allIter << ")" << std::endl;
	    	}
	    }
	    if (myScreen != NULL)
    	{
    		std::stringstream soutput;
    		soutput << "Final residuum " << myCG.residuum << "; with " << myCG.getNumberIterations() << " Iterations (Total Iter.: " << allIter << ")";

    		if (i < this->nMaxIterations-1)
    		{
    			myScreen->update((size_t)(((double)(i+1)*100.0)/((double)this->nMaxIterations)), soutput.str());
    		}
    		else
    		{
    			myScreen->update(100, soutput.str());
    		}
    	}


	    // Do some adjustments on the boundaries if needed
		if (this->ExMode == "ExEul")
		{
			SystemMatrix.finishTimestep(alpha);
		}


		// Create pictures of the animation, if specified
	    if ((this->bAnimation == true) && (i%animationStep == 0))
		{
			// Build filename
			std::string tFilename = "00000000000000000000000000000000";
			std::stringstream number;
			number << i;
			tFilename.append(number.str());
			tFilename = tFilename.substr(tFilename.length()-14,14);
			tFilename.append(".gnuplot");

			// Print grid to file
			GridPrinter myPrinter(*SystemMatrix.getGrid());
			myPrinter.printGrid(alpha, tFilename, this->evalsAnimation);
		}
	}

	// write some empty lines to console
    if (myScreen != NULL)
	{
    	myScreen->writeEmptyLines(2);
	}
}

}
