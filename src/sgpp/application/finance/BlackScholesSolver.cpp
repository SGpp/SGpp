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

#include "application/finance/BlackScholesSolver.hpp"

namespace sg
{

BlackScholesSolver::BlackScholesSolver()
{
	bStochasticDataAlloc = false;
	bGridConstructed = false;
}

BlackScholesSolver::~BlackScholesSolver()
{
	if (bStochasticDataAlloc)
	{
		delete mus;
		delete sigmas;
		delete rhos;
	}
	if (bGridConstructed)
	{
		delete myGrid;
	}
}

void BlackScholesSolver::constructGrid(BoundingBox& myBoundingBox, size_t level)
{
	dim = myBoundingBox.getDimensions();
	levels = level;

	myGrid = new LinearTrapezoidBoundaryGrid(dim, true);

	GridGenerator* myGenerator = myGrid->createGridGenerator();
	myGenerator->regular(levels);
	delete myGenerator;



	bGridConstructed = true;
}

void BlackScholesSolver::constructGrid(std::string tfilename)
{
	// @todo (heinecke) implement this
}

void BlackScholesSolver::setStochasticData(DataVector& mus, DataVector& sigmas, DataVector& rhos, double r)
{
	this->mus = DataVector(mus);
	this->sigmas = DataVector(sigmas);
	this->rhos = DataVector(rhos);
	this->r = r;

	bStochasticDataAlloc = true;
}

void BlackScholesSolver::solveEuler(size_t numTimesteps, double timestepsize, DataVector& alpha)
{
	// @todo (heinecke) implement this
}

void BlackScholesSolver::printGrid(DataVector& alpha, double resolution, std::string filename)
{
	DimensionBoundary dimOne;
	DimensionBoundary dimTwo;
	std::ofstream fileout;

	if (bGridConstructed)
	{
		if (dim > 2)
		{
			// @todo (heinecke) thrown an application exception
		}
		else
		{
			// Open filehandle
			fileout.open(tfilename.str().c_str());
			OperationEval* myEval = myGrid->createOperationEvalBB();

			if (dim == 1)
			{
				dimOne = myGrid->getBoundingBox()->getBoundary(0);

				for (double i = dimOne.leftBoundary; i < dimOne.rightBoundary; i+=resolution)
				{
					std::vector<double> point;
					point.push_back(i);
					fileout << i << " " << myEval->eval(alpha,point) << std::endl;
				}
			}
			else
			{
				dimOne = myGrid->getBoundingBox()->getBoundary(0);
				dimTwo = myGrid->getBoundingBox()->getBoundary(1);

				for (double i = dimOne.leftBoundary; i < dimOne.rightBoundary; i+=resolution)
				{
					for (double j = dimTwo.leftBoundary; j < dimTwo.rightBoundary; j+=resolution)
					{
						std::vector<double> point;
						point.push_back(i);
						point.push_back(j);
						fileout << i << " " << j << " " << myEval->eval(alpha,point) << std::endl;
					}
				}
			}

			delete myEval;
			// close filehandle
			fileout.close();
		}
	}
	else
	{
		// @todo (heinecke) thrown an application exception
	}
}

}
