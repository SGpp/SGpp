/*****************************************************************************/
/* This file is part of sgpp, a program package making use of spatially      */
/* adaptive sparse grids to solve numerical problems                         */
/*                                                                           */
/* Copyright (C) 2009-2010 Alexander Heinecke (Alexander.Heinecke@mytum.de)  */
/*                                                                           */
/* sgpp is free software; you can redistribute it and/or modify              */
/* it under the terms of the GNU Lesser General Public License as published  */
/* by the Free Software Foundation; either version 3 of the License, or      */
/* (at your option) any later version.                                       */
/*                                                                           */
/* sgpp is distributed in the hope that it will be useful,                   */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of            */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             */
/* GNU Lesser General Public License for more details.                       */
/*                                                                           */
/* You should have received a copy of the GNU Lesser General Public License  */
/* along with sgpp; if not, write to the Free Software                       */
/* Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA */
/* or see <http://www.gnu.org/licenses/>.                                    */
/*****************************************************************************/

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

void CrankNicolson::solve(SLESolver& LinearSystemSolver, OperationODESolverSystem& System, bool verbose)
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
		System.finishTimestep();
	}

	// write some empty lines to console
    if (myScreen != NULL)
	{
    	myScreen->writeEmptyLines(2);
	}
}

}
