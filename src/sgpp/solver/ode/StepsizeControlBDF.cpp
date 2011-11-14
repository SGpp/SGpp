/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/

#include "base/grid/common/DirichletUpdateVector.hpp"
#include "solver/ode/StepsizeControlBDF.hpp"

#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <cmath>

namespace sg
{
namespace solver
{

StepsizeControlBDF::StepsizeControlBDF(size_t nTimesteps, double timestepSize, double eps, sg::base::ScreenOutput* screen) : ODESolver(nTimesteps, timestepSize), myScreen(screen)
{
	this->residuum = 0.0;
	this->myEps = eps;
}

StepsizeControlBDF::~StepsizeControlBDF()
{
}

void StepsizeControlBDF::solve(SLESolver& LinearSystemSolver, sg::pde::OperationParabolicPDESolverSystem& System, bool bIdentifyLastStep, bool verbose)
{
	size_t allIter = 0;
	sg::base::DataVector* rhs;
	sg::base::DataVector YkBDF2(System.getGridCoefficients()->getSize());
	sg::base::DataVector YkF23(System.getGridCoefficients()->getSize());

	double tmp_timestepsize = this->myEpsilon;
	double tmp_timestepsize_old = tmp_timestepsize;
	double tmp_timestepsize_new = tmp_timestepsize;
	double epsilon = this->myEps;


	double maxTimestep = this->nMaxIterations*tmp_timestepsize;

	size_t maxIter = this->nMaxIterations*100000;

	double time = 0.0;


    std::stringstream fnsstream;
    fnsstream << "Time_" << "SCBDF" << this->myEps << ".gnuplot";
    std::string filename = fnsstream.str();

    std::ofstream fileout;

    //fileout.open(filename.c_str());
    fileout.open(filename.c_str(), std::ofstream::app); // apend to file
	fileout << std::endl;


	for (size_t i = 0; i < maxIter && time < maxTimestep; i++)
	{

		System.setTimestepSize(tmp_timestepsize);


		if(i > 0)
			System.setODESolver("AdBas");
		else
			System.setODESolver("ExEul");

		// generate right hand side
		rhs = System.generateRHS();

		// solve the system of the current timestep
		LinearSystemSolver.solve(System, *System.getGridCoefficientsForCG(), *rhs, true, false);
		System.finishTimestep();

		System.getGridCoefficientsForSC(YkF23);
		System.abortTimestep();

		if(i > 0)
			System.setODESolver("BDF2");
		else
			System.setODESolver("ImEul");

		// generate right hand side
		rhs = System.generateRHS();

		// solve the system of the current timestep
		LinearSystemSolver.solve(System, *System.getGridCoefficientsForCG(), *rhs, true, false);
		System.finishTimestep();
		System.getGridCoefficientsForSC(YkBDF2);


/*
		double epsilon = 0.001;
		YkBDF2.sub(YkF23);
		double u  = std::max(YkBDF2.max(),-YkBDF2.min());
		double tD = tmp_timestepsize/tmp_timestepsize_old;
		double C1 = (1.0+1.0/tD)/6.0;
		double Cp = -(1.0+tD)*(1.0+tD)/(6*tD*(1.0+2.0*tD));
		double alpha0 = (1+2*tD)/(1+tD); //??
		double deltaY = alpha0*C1*u/(tmp_timestepsize*(C1-Cp));
		tmp_timestepsize_new = tmp_timestepsize * sqrt(epsilon/deltaY);
*/
	    YkBDF2.sub(YkF23);
	    double u  = sqrt(YkBDF2.dotProduct(YkBDF2));

	    // double deltaY = u/(3.0*(1.0+tmp_timestepsize/tmp_timestepsize_old));

	    tmp_timestepsize_new = tmp_timestepsize * pow(epsilon/u,(double)1.0/(double)3.0);



		if (0.8*tmp_timestepsize > tmp_timestepsize_new) {
			tmp_timestepsize = tmp_timestepsize_new;
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
			if(tmp_timestepsize_new > tmp_timestepsize)
			{
				tmp_timestepsize = tmp_timestepsize_new;
			}

			// avoid small last time steps
			if(maxTimestep-time<1.3*tmp_timestepsize){
				tmp_timestepsize = maxTimestep-time;
			}
			// adapt size of last time step
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
}
