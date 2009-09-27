/*****************************************************************************/
/* This file is part of sgpp, a program package making use of spatially      */
/* adaptive sparse grids to solve numerical problems                         */
/*                                                                           */
/* Copyright (C) 2009 Alexander Heinecke (Alexander.Heinecke@mytum.de)       */
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

#include "solver/ode/Euler.hpp"
#include "operation/common/OperationEval.hpp"
#include "tools/common/GridPrinter.hpp"
#include "solver/sle/BiCGStab.hpp"
#include "solver/sle/ConjugateGradients.hpp"
#include "grid/common/DirichletUpdateVector.hpp"

#include <iostream>
#include <string>
#include <sstream>

namespace sg
{

Euler::Euler(std::string Mode, size_t imax, double timestepSize, size_t iMaxCG, double epsilonCG, bool generateAnimation) : ODESolver(imax, timestepSize), maxCGIterations(iMaxCG), epsilonCG(epsilonCG), bAnimation(generateAnimation), ExMode(Mode)
{
	this->residuum = 0.0;

	if (Mode != "ExEul" || Mode != "ImEul")
	{
		// @todo (heinecke) throw solver exception
	}
}

Euler::~Euler()
{
}

void Euler::solve(OperationSolverMatrix& SystemMatrix, DataVector& alpha, bool verbose)
{
	DataVector rhs(alpha.getSize());
	DataVector saveAlpha(alpha.getSize());
    BiCGStab myCG(this->maxCGIterations, this->epsilonCG);
    DirichletUpdateVector myDirichletUpdate(SystemMatrix.getGrid()->getStorage());

	for (size_t i = 0; i < this->nMaxIterations; i++)
	{
		rhs.setAll(0.0);

		saveAlpha = alpha;

		SystemMatrix.generateRHS(alpha, rhs);

		//myDirichletUpdate.setBoundariesToZero(alpha);
		//myDirichletUpdate.setBoundariesToZero(rhs);

	    myCG.solve(SystemMatrix, alpha, rhs, true, false, -1.0);
	    if (verbose)
	    {
	    	std::cout << "Final residuum " << myCG.residuum << "; with " << myCG.getNumberIterations() << " Iterations" << std::endl;
	    }

	    //myDirichletUpdate.applyDirichletConditions(alpha, saveAlpha);
	    //std::cout << "New alpha: " << alpha.toString() << std::endl;

		if (this->bAnimation == true && (i%(this->nMaxIterations/1500)) == 0)
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
			myPrinter.printGrid(alpha, tFilename, 1000);
		}
	}
}

}
