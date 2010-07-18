/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/

#include "grid/common/DirichletUpdateVector.hpp"
#include "solver/ode/VarTimestep.hpp"
#include "operation/common/OperationEval.hpp"
#include "tools/common/GridPrinter.hpp"
#include "exception/solver_exception.hpp"

#include <iostream>
#include <string>
#include <sstream>
#include <fstream>

namespace sg
{

VarTimestep::VarTimestep(size_t imax, double timestepSize, bool generateAnimation, size_t numEvalsAnimation, ScreenOutput* screen) : ODESolver(imax, timestepSize), bAnimation(generateAnimation), evalsAnimation(numEvalsAnimation), myScreen(screen)
{
	this->residuum = 0.0;
	this->TimestepSize = timestepSize;

}

VarTimestep::~VarTimestep()
{
}

void VarTimestep::solve(SLESolver& LinearSystemSolver, OperationODESolverSystem& System, bool bIdentifyLastStep, bool verbose)
{
	size_t allIter = 0;
    DataVector* rhs;
    DataVector YkAdBas(System.getGridCoefficients()->getSize());
    DataVector YkImEul(System.getGridCoefficients()->getSize());

    double tmp_timestepsize = this->TimestepSize;
    double tmp_timestepsize_old = this->TimestepSize;
    double tmp_timestepsize_new = this->TimestepSize;

    // Do some animation creation exception handling
    size_t animationStep = this->nMaxIterations/1500;
    if (animationStep == 0)
    {
    	animationStep = 1;
    }

    double maxTimestep = this->nMaxIterations*this->TimestepSize;

    size_t maxIter = this->nMaxIterations*1000;

    double time = 0.0;


    std::ofstream fileout;

    fileout.open("time.gnuplot");



	for (size_t i = 0; i < maxIter && time < maxTimestep; i++)
	{

		System.setTimestepSize(tmp_timestepsize);
		YkAdBas.setAll(0.0);
		YkImEul.setAll(0.0);
		System.setODESolver("AdBas");
		// generate right hand side
		rhs = System.generateRHS();

		// solve the system of the current timestep
		LinearSystemSolver.solve(System, *System.getGridCoefficientsForCG(), *rhs, true, false, -1.0);

		System.finishTimestep();
		YkAdBas.add(*System.getGridCoefficients());
		System.abortTimestep();

		System.setODESolver("ImEul");

		rhs = System.generateRHS();
		LinearSystemSolver.solve(System, *System.getGridCoefficientsForCG(), *rhs, true, false, -1.0);

		YkImEul.add(*System.getGridCoefficients());

	    double epsilon = pow(10,-4);
	    YkImEul.sub(YkAdBas);
	    double tmp  = sqrt(YkImEul.dotProduct(YkImEul));
	    double deltaY = tmp/(3.0*(1.0+tmp_timestepsize/tmp_timestepsize_old));

	    tmp_timestepsize_new = tmp_timestepsize * pow(epsilon/deltaY,(double)1.0/(double)3.0);

	    fileout << i << " " << (tmp_timestepsize_new-tmp_timestepsize) << " " << tmp << std::endl;

	    if (0.8*tmp_timestepsize > tmp_timestepsize_new) {
	    	tmp_timestepsize = tmp_timestepsize_new;
	    	System.abortTimestep();
	    	allIter += LinearSystemSolver.getNumberIterations();

	    }
	    else {
	    	time += tmp_timestepsize;
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

	    		soutput << " Final residuum " << LinearSystemSolver.residuum << "; with " << LinearSystemSolver.getNumberIterations() << " Iterations (Total Iter.: " << allIter << ")";

	    		if (i < this->nMaxIterations-1)
	    		{
	    			myScreen->update((size_t)(((double)(time)*100.0)/((double)maxTimestep)), soutput.str());
	    		}
	    		else
	    		{
	    			myScreen->update(100, soutput.str());
	    		}
	    	}
	    	/** FINIS TS*/
	    	System.finishTimestep();
	    	System.saveAlpha();




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
				GridPrinter myPrinter(*System.getGrid());
				myPrinter.printGrid(*System.getGridCoefficients(), tFilename, this->evalsAnimation);
			}

			tmp_timestepsize_old = tmp_timestepsize;
			if(tmp_timestepsize_new > tmp_timestepsize)
			{
				tmp_timestepsize = tmp_timestepsize_new;
			}

			tmp_timestepsize = std::min(tmp_timestepsize,maxTimestep-time);

	    }

	}
	fileout.close();



	// write some empty lines to console
    if (myScreen != NULL)
	{
    	myScreen->writeEmptyLines(2);
	}
}

}
