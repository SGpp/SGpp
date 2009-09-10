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

#include "algorithm/finance/BlackScholesTimestepMatrix.hpp"
#include "application/finance/BlackScholesSolver.hpp"
#include "solver/ode/ExplicitEuler.hpp"
#include "solver/ode/CrankNicolson.hpp"
#include "grid/Grid.hpp"
#include "stdlib.h"
#include <sstream>

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

void BlackScholesSolver::constructGrid(BoundingBox& BoundingBox, size_t level)
{
	dim = BoundingBox.getDimensions();
	levels = level;

	myGrid = new LinearTrapezoidBoundaryGrid(BoundingBox);

	GridGenerator* myGenerator = myGrid->createGridGenerator();
	myGenerator->regular(levels);
	delete myGenerator;

	myBoundingBox = myGrid->getBoundingBox();
	myGridStorage = myGrid->getStorage();

	bGridConstructed = true;
}

void BlackScholesSolver::constructGrid(std::string tfilename)
{
	// @todo (heinecke) implement this
}

void BlackScholesSolver::setStochasticData(DataVector& mus, DataVector& sigmas, DataVector& rhos, double r)
{
	this->mus = new DataVector(mus);
	this->sigmas = new DataVector(sigmas);
	this->rhos = new DataVector(rhos);
	this->r = r;

	bStochasticDataAlloc = true;
}

void BlackScholesSolver::solveEuler(size_t numTimesteps, double timestepsize, DataVector& alpha)
{
	if (bGridConstructed && bStochasticDataAlloc)
	{
		ExplicitEuler* myEuler = new ExplicitEuler(numTimesteps, timestepsize);
		BlackScholesTimestepMatrix* myBSMatrix = new BlackScholesTimestepMatrix(*myGrid, *this->mus, *this->sigmas, *this->rhos, r, numTimesteps, false);

		myEuler->solve(*myBSMatrix, alpha, false);

		delete myBSMatrix;
		delete myEuler;
	}
	else
	{
		// @todo (heinecke) through an application exception here
	}
}

void BlackScholesSolver::solveCrankNicolson(size_t numTimesteps, double timestepsize, size_t maxCGIterations, double epsilonCG, DataVector& alpha)
{
	if (bGridConstructed && bStochasticDataAlloc)
	{
		CrankNicolson* myCN = new CrankNicolson(numTimesteps, timestepsize, maxCGIterations, epsilonCG);
		BlackScholesTimestepMatrix* myBSMatrix = new BlackScholesTimestepMatrix(*myGrid, *this->mus, *this->sigmas, *this->rhos, r, numTimesteps, true);

		myCN->solve(*myBSMatrix, alpha, false);

		delete myBSMatrix;
		delete myCN;
	}
	else
	{
		// @todo (heinecke) through an application exception here
	}
}

void BlackScholesSolver::printGrid(DataVector& alpha, double resolution, std::string tfilename)
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
			fileout.open(tfilename.c_str());
			OperationEval* myEval = myGrid->createOperationEval();

			if (dim == 1)
			{
				dimOne = myGrid->getBoundingBox()->getBoundary(0);

				for (double i = dimOne.leftBoundary; i <= dimOne.rightBoundary; i+=resolution)
				{
					std::vector<double> point;
					point.push_back(i);
					fileout << i << " " << myEval->eval(alpha,point) << std::endl;
				}
			}
			else if (dim == 2)
			{
				dimOne = myGrid->getBoundingBox()->getBoundary(0);
				dimTwo = myGrid->getBoundingBox()->getBoundary(1);

				for (double i = dimOne.leftBoundary; i <= dimOne.rightBoundary; i+=resolution)
				{
					for (double j = dimTwo.leftBoundary; j <= dimTwo.rightBoundary; j+=resolution)
					{
						std::vector<double> point;
						point.push_back(i);
						point.push_back(j);
						fileout << i << " " << j << " " << myEval->eval(alpha,point) << std::endl;
					}
					fileout << std::endl;
				}
			}
			else
			{
				// @todo (heinecke) thrown an application exception
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

void BlackScholesSolver::initGridWithPayoff(DataVector& alpha, double* strike)
{
	double tmp;
	double tmp2;

	if (bGridConstructed)
	{
		if (dim == 1)
		{
			for (size_t i = 0; i < myGrid->getStorage()->size(); i++)
			{
				tmp = atof(myGridStorage->get(i)->getCoordsStringBB(*myBoundingBox).c_str());
				alpha[i] = get1DPayoffValue(tmp, strike[0]);
			}

			OperationHierarchisation* myHierarchisation = myGrid->createOperationHierarchisation();
			myHierarchisation->doHierarchisation(alpha);
			delete myHierarchisation;
		}
		else if (dim == 2)
		{
			for (size_t i = 0; i < myGrid->getStorage()->size(); i++)
			{
					std::string coords = myGridStorage->get(i)->getCoordsStringBB(*myBoundingBox);
					std::stringstream coordsStream(coords);

					coordsStream >> tmp;
					coordsStream >> tmp2;

					alpha[i] = max(get1DPayoffValue(tmp, strike[0]),get1DPayoffValue(tmp2, strike[1]));
			}

			OperationHierarchisation* myHierarchisation = myGrid->createOperationHierarchisation();
			myHierarchisation->doHierarchisation(alpha);
			delete myHierarchisation;
		}
		else
		{
			// @todo (heinecke) thrown an application exception
		}
	}
	else
	{
		// @todo (heinecke) throw an application exception
	}
}

size_t BlackScholesSolver:: getNumberGridPoints()
{
	if (bGridConstructed)
	{
		return myGridStorage->size();
	}
	else
	{
		// @todo (heinecke) throw an application exception
		return 0;
	}
}

double BlackScholesSolver::get1DPayoffValue(double assetValue, double strike)
{
	if (assetValue <= strike)
	{
		return 0.0;
	}
	else
	{
		return assetValue - strike;
	}
}

}
