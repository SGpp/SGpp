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

#include "solver/ode/CrankNicolson.hpp"
#include "solver/sle/ConjugateGradients.hpp"

namespace sg
{

CrankNicolson::CrankNicolson(size_t nTimesteps, double timestepSize, size_t iMaxCG, double epsilonCG) : ODESolver(nTimesteps, timestepSize), maxCGIterations(iMaxCG), epsilonCG(epsilonCG)
{
	this->residuum = 0.0;
}

CrankNicolson::~CrankNicolson()
{
}

void CrankNicolson::solve(OperationSolverMatrix& SystemMatrix, DataVector& alpha, bool verbose)
{
	DataVector rhs(alpha.getSize());
    ConjugateGradients myCG(this->maxCGIterations, this->epsilonCG);

	for (size_t i = 0; i < this->nMaxIterations; i++)
	{
		rhs.setAll(0.0);

		SystemMatrix.generateRHS(alpha, rhs);

	    myCG.solve(SystemMatrix, alpha, rhs, false, true, -1.0);
	}
}

}
