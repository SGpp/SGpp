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

#include "algorithm/pde/HeatEquationTimestepMatrix.hpp"
#include "application/pde/HeatEquationSolver.hpp"
#include "solver/ode/Euler.hpp"
#include "solver/ode/CrankNicolson.hpp"
#include "grid/Grid.hpp"
#include "stdlib.h"
#include <sstream>

namespace sg
{

HeatEquationSolver::HeatEquationSolver()
{
	bGridConstructed = false;
}

HeatEquationSolver::~HeatEquationSolver()
{
	if (bGridConstructed)
	{
		delete myGrid;
	}
}

void HeatEquationSolver::constructGrid(BoundingBox& BoundingBox, size_t level)
{
	dim = BoundingBox.getDimensions();
	levels = level;

	//myGrid = new LinearTrapezoidBoundaryGrid(BoundingBox);
	myGrid = new LinearGrid(dim);

	GridGenerator* myGenerator = myGrid->createGridGenerator();
	myGenerator->regular(levels);
	delete myGenerator;

	myBoundingBox = myGrid->getBoundingBox();
	myGridStorage = myGrid->getStorage();

	bGridConstructed = true;
}

void HeatEquationSolver::solveExplicitEuler(size_t numTimesteps, double timestepsize, size_t maxCGIterations, double epsilonCG, double a, DataVector& alpha, bool generateAnimation)
{
	if (bGridConstructed)
	{
		Euler* myEuler = new Euler("ExEul", numTimesteps, timestepsize, maxCGIterations, epsilonCG, generateAnimation);
		HeatEquationTimestepMatrix* myHEMatrix = new HeatEquationTimestepMatrix(*myGrid, a, timestepsize, "ExEul");

		myEuler->solve(*myHEMatrix, alpha, true);

		delete myHEMatrix;
		delete myEuler;
	}
	else
	{
		// @todo (heinecke) through an application exception here
	}
}

void HeatEquationSolver::solveImplicitEuler(size_t numTimesteps, double timestepsize, size_t maxCGIterations, double epsilonCG, double a, DataVector& alpha, bool generateAnimation)
{
	if (bGridConstructed)
	{
		Euler* myEuler = new Euler("ImEul", numTimesteps, timestepsize, maxCGIterations, epsilonCG, generateAnimation);
		HeatEquationTimestepMatrix* myHEMatrix = new HeatEquationTimestepMatrix(*myGrid, a, timestepsize, "ImEul");

		myEuler->solve(*myHEMatrix, alpha, false);

		delete myHEMatrix;
		delete myEuler;
	}
	else
	{
		// @todo (heinecke) through an application exception here
	}
}

void HeatEquationSolver::solveCrankNicolson(size_t numTimesteps, double timestepsize, size_t maxCGIterations, double epsilonCG, double a, DataVector& alpha)
{
	if (bGridConstructed)
	{
		CrankNicolson* myCN = new CrankNicolson(numTimesteps, timestepsize, maxCGIterations, epsilonCG);
		HeatEquationTimestepMatrix* myHEMatrix = new HeatEquationTimestepMatrix(*myGrid, a, timestepsize, "CrNic");

		myCN->solve(*myHEMatrix, alpha, false);

		delete myHEMatrix;
		delete myCN;
	}
	else
	{
		// @todo (heinecke) through an application exception here
	}
}

void HeatEquationSolver::printGrid(DataVector& alpha, double PointesPerDimension, std::string tfilename)
{
	GridPrinter myPrinter(*this->myGrid);
	myPrinter.printGrid(alpha, tfilename, PointesPerDimension);
}

void HeatEquationSolver::initGridWithSingleHeat(DataVector& alpha, double heat)
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

				if (tmp == 0.5)
				{
					alpha[i] = heat;
				}
				else
				{
					alpha[i] = 0.0;
				}
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

					if (tmp == 0.5 && tmp2 == 0.5)
					{
						alpha[i] = heat;
					}
					else
					{
						alpha[i] = 0.0;
					}
			}

			//std::cout << alpha.toString() << std::endl;

			OperationHierarchisation* myHierarchisation = myGrid->createOperationHierarchisation();
			myHierarchisation->doHierarchisation(alpha);

			//std::cout << alpha.toString() << std::endl;

			//myHierarchisation->doDehierarchisation(alpha);

			//std::cout << alpha.toString() << std::endl;
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

size_t HeatEquationSolver:: getNumberGridPoints()
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

void HeatEquationSolver::initGridWithSmoothHeat(DataVector& alpha, double mu, double sigma)
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

				alpha[i] = (1.0/(sigma*2.0*3.145))*exp((-0.5)*((tmp-mu)/sigma)*((tmp-mu)/sigma));
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

				alpha[i] = (1.0/(sigma*2.0*3.145))*exp((-0.5)*((tmp-mu)/sigma)*((tmp-mu)/sigma)) * (1.0/(sigma*2.0*3.145))*exp((-0.5)*((tmp2-mu)/sigma)*((tmp2-mu)/sigma));
			}

			OperationHierarchisation* myHierarchisation = myGrid->createOperationHierarchisation();
			myHierarchisation->doHierarchisation(alpha);
			delete myHierarchisation;
		}
		else
		{

		}
	}
}

void HeatEquationSolver::initGridWithConstantHeat(DataVector& alpha, double constHeat)
{
	double tmp;
	//double tmp2;

	if (bGridConstructed)
	{
		if (dim == 1)
		{
			for (size_t i = 0; i < myGrid->getStorage()->size(); i++)
			{
				tmp = atof(myGridStorage->get(i)->getCoordsStringBB(*myBoundingBox).c_str());

				alpha[i] = constHeat;
			}

			OperationHierarchisation* myHierarchisation = myGrid->createOperationHierarchisation();
			myHierarchisation->doHierarchisation(alpha);
			delete myHierarchisation;
		}
		else
		{

		}
	}
}

}
