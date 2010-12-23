/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/

#include "grid/common/DirichletUpdateVector.hpp"
#include "solver/ode/StepsizeControlH.hpp"
#include "operation/common/OperationEval.hpp"
#include "tools/common/GridPrinter.hpp"
#include "exception/solver_exception.hpp"

#include <iostream>
#include <string>
#include <sstream>
#include <fstream>

namespace sg
{

StepsizeControlH::StepsizeControlH(size_t imax, double timestepSize, double eps, ScreenOutput* screen) : ODESolver(imax, timestepSize), myScreen(screen)
{
	this->residuum = 0.0;
	this->myEps = eps;
}

StepsizeControlH::~StepsizeControlH()
{
}

void StepsizeControlH::solve(SLESolver& LinearSystemSolver, OperationODESolverSystem& System, bool bIdentifyLastStep, bool verbose)
{
	size_t allIter = 0;
    DataVector* rhs;
    DataVector YkCrNic(System.getGridCoefficients()->getSize());
    DataVector YkCrNicHalf(System.getGridCoefficients()->getSize());

    double tmp_timestepsize = this->myEpsilon;
    double tmp_timestepsize_old = tmp_timestepsize;
    double tmp_timestepsize_new = tmp_timestepsize;
    double epsilon = this->myEps;


    double maxTimestep = this->nMaxIterations*tmp_timestepsize;

    size_t maxIter = this->nMaxIterations*1000000;

    double time = 0.0;

    std::stringstream fnsstream;
    fnsstream << "Time_" << "SCH" << this->myEps << ".gnuplot";
    std::string filename = fnsstream.str();

    std::ofstream fileout;

    fileout.open(filename.c_str());


	System.setODESolver("CrNic");
	for (size_t i = 0; i < maxIter && time < maxTimestep; i++)
	{
		System.setTimestepSize(tmp_timestepsize);

		// generate right hand side
		rhs = System.generateRHS();

		// solve the system of the current timesteps
		LinearSystemSolver.solve(System, *System.getGridCoefficientsForCG(), *rhs, true, false, -1.0);

		System.finishTimestep();

		System.getGridCoefficientsForSC(YkCrNic);

		System.abortTimestep();

		System.setTimestepSize(tmp_timestepsize/2.0);

		// generate right hand side
		rhs = System.generateRHS();

		// solve the system of the current timesteps
		LinearSystemSolver.solve(System, *System.getGridCoefficientsForCG(), *rhs, true, false, -1.0);
		System.finishTimestep();

		rhs = System.generateRHS();
		LinearSystemSolver.solve(System, *System.getGridCoefficientsForCG(), *rhs, true, false, -1.0);
		System.finishTimestep();
		System.getGridCoefficientsForSC(YkCrNicHalf);

	    YkCrNicHalf.sub(YkCrNic);
	    double tmp  = sqrt(YkCrNicHalf.dotProduct(YkCrNicHalf));


	    double deltaY = 3.0*tmp/4.0;

	    tmp_timestepsize_new = tmp_timestepsize * pow(epsilon/deltaY,(double)1.0/(double)3.0);

	    if (tmp >= epsilon) {
	    	tmp_timestepsize = tmp_timestepsize*0.9;
	    	System.abortTimestep();
	    	allIter += LinearSystemSolver.getNumberIterations();

	    }
	    else {
	    	fileout << i << " " << (tmp_timestepsize_new-tmp_timestepsize) << " " << time << " " << tmp_timestepsize << std::endl;
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

	    	System.saveAlpha();
			tmp_timestepsize_old = tmp_timestepsize;
			tmp_timestepsize = tmp_timestepsize_new;
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
