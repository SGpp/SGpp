/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/

#include "grid/common/DirichletUpdateVector.hpp"
#include "solver/ode/StepsizeControlEJ.hpp"
#include "solver/ode/StepsizeControlEJ.hpp"
#include "grid/Grid.hpp"
#include "operation/common/OperationEval.hpp"
#include "tools/common/GridPrinter.hpp"
#include "exception/solver_exception.hpp"


#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <cmath>

namespace sg
{
namespace solver
{

StepsizeControlEJ::StepsizeControlEJ(size_t nTimesteps, double timestepSize, double eps, double sc, sg::base::ScreenOutput* screen) : ODESolver(nTimesteps, timestepSize), myScreen(screen)
{
	this->residuum = 0.0;
	this->myEps = eps;
	this->mySC = sc;
}

StepsizeControlEJ::~StepsizeControlEJ()
{
}

void StepsizeControlEJ::solve(SLESolver& LinearSystemSolver, sg::pde::OperationParabolicPDESolverSystem& System, bool bIdentifyLastStep, bool verbose)
{
	size_t allIter = 0;
	sg::base::DataVector* rhs;
	sg::base::DataVector YkImEul(System.getGridCoefficients()->getSize());
	sg::base::DataVector YkImEulOld(System.getGridCoefficients()->getSize());

	double tmp_timestepsize = this->myEpsilon;
	double tmp_timestepsize_old = tmp_timestepsize;
	double tmp_timestepsize_new = tmp_timestepsize;
	double epsilon = this->myEps;


	double maxTimestep = this->nMaxIterations*tmp_timestepsize;

	size_t maxIter = this->nMaxIterations*100000;

	double time = 0.0;

    std::stringstream fnsstream;
    fnsstream << "Time_" << "SCEJ" << this->myEps << "_" << this->mySC << ".gnuplot";
    std::string filename = fnsstream.str();

    std::ofstream fileout;

    //fileout.open(filename.c_str());
    fileout.open(filename.c_str(), std::ofstream::app); // apend to file
    fileout << std::endl;

	System.getGridCoefficientsForSC(YkImEul);

       double sc = this->mySC;

	System.setODESolver("CrNic");

	for (size_t i = 0; i < maxIter && time < maxTimestep; i++)
	{

		YkImEulOld.setAll(0.0);
		YkImEulOld.add(YkImEul);

		System.setTimestepSize(tmp_timestepsize);

		// generate right hand side
		rhs = System.generateRHS();

		// solve the system of the current timestep
		LinearSystemSolver.solve(System, *System.getGridCoefficientsForCG(), *rhs, true, false);

		System.finishTimestep();

		System.getGridCoefficientsForSC(YkImEul);

		double *OldData = YkImEulOld.getPointer();
		double *Data = YkImEul.getPointer();

	       double max = 0.0;
              double diff = 0.0;
              int idx = -1;
              double divisor = 0.0;

	       // calculate the max norm
		for(int j=0;j<System.getGridCoefficientsForCG()->getSize();j++)
		{
			double t2 = std::max(fabs(Data[j]),fabs(OldData[j]));
			double tmpData = fabs(Data[j]-OldData[j])/std::max(sc,t2);
			if (max < fabs(tmpData)){
				max = fabs(tmpData);
                            // just for output
                            idx=j;  
                            diff = fabs(Data[j]-OldData[j]);
                            divisor = std::max(sc,t2);
                     }
		}

		double u  = max;

		tmp_timestepsize_new = tmp_timestepsize * epsilon/u;

		if (u >= epsilon) {
			tmp_timestepsize = tmp_timestepsize/2.0;
			System.abortTimestep();
			allIter += LinearSystemSolver.getNumberIterations();
		}
		else {
			fileout << i << " " << (tmp_timestepsize_new-tmp_timestepsize) << " " << time << " " << tmp_timestepsize  << " " <<  u << " " << diff << " " << divisor << " " << idx << std::endl;
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
			//	soutput << " asd" << i << " " << (tmp_timestepsize_new-tmp_timestepsize)  << " " <<  index<< " " <<  u << " " << tmp_timestepsize << " ";

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

			// avoid small last time steps
			if(maxTimestep-time<1.3*tmp_timestepsize){
				tmp_timestepsize = maxTimestep-time;
			}
			// adapt size of last timestep
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
