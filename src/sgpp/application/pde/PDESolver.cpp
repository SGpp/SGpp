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

#include "application/pde/PDESolver.hpp"
#include "grid/Grid.hpp"
#include "exception/application_exception.hpp"
#include <sstream>
#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <cmath>

namespace sg
{

PDESolver::PDESolver()
{
	bGridConstructed = false;
}

PDESolver::~PDESolver()
{
	if (bGridConstructed)
	{
		delete myGrid;
	}
}

void PDESolver::deleteGrid()
{
	if (bGridConstructed)
	{
		delete myGrid;
		bGridConstructed = false;
		myBoundingBox = NULL;
		myGridStorage = NULL;
	}
	else
	{
		throw new application_exception("PDESolver::deleteGrid : The grid wasn't initialized before!");
	}
}

void PDESolver::setGrid(const std::string& serializedGrid)
{
	if (bGridConstructed)
	{
		delete myGrid;
		bGridConstructed = false;
		myBoundingBox = NULL;
		myGridStorage = NULL;
	}

	myGrid = Grid::unserialize(serializedGrid);

	myBoundingBox = myGrid->getBoundingBox();
	myGridStorage = myGrid->getStorage();

	dim = myGrid->getStorage()->dim();
	levels = 0;

	bGridConstructed = true;
}

std::string PDESolver::getGrid()
{
	std::string gridSer = "";

	if (bGridConstructed)
	{
		// Serialize the grid
		myGrid->serialize(gridSer);
	}
	else
	{
		throw new application_exception("PDESolver::getGrid : The grid wasn't initialized before!");
	}

	return gridSer;
}

void PDESolver::refineInitialGridSurplus(DataVector& alpha, double dPercentage)
{
	size_t nRefinements = static_cast<size_t>(static_cast<double>(alpha.getNumberNonZero())*(dPercentage/100));

	if (bGridConstructed)
	{
		SurplusRefinementFunctor* myRefineFunc = new SurplusRefinementFunctor(&alpha, nRefinements, 0.0);

		myGrid->createGridGenerator()->refine(myRefineFunc);

		delete myRefineFunc;

		alpha.resize(myGridStorage->size());
	}
	else
	{
		throw new application_exception("PDESolver::refineIntialGridSurplus : The grid wasn't initialized before!");
	}
}

void PDESolver::constructGridBonn(std::string tfilename, DataVector& emptyAlpha, bool& ishierarchized)
{
	IOToolBonnSG* myImporter = new IOToolBonnSG();
	std::string serGrid;

	// test if emptyAlpha is really empty
	if (emptyAlpha.getSize() != 0)
	{
		throw new application_exception("PDESolver::constructGrid : A non-empty coefficients' vector was used in gird construction from file!");
	}

	myImporter->readFile(tfilename, serGrid, emptyAlpha, ishierarchized);

	myGrid = Grid::unserialize(serGrid);

	myBoundingBox = myGrid->getBoundingBox();
	myGridStorage = myGrid->getStorage();

	// Set every boundary to dirichlet boundaries
	for (size_t i = 0; i < myGridStorage->dim(); i++)
	{
		DimensionBoundary myDimBound = myBoundingBox->getBoundary(i);

		myDimBound.bDirichletLeft = true;
		myDimBound.bDirichletRight = true;

		myBoundingBox->setBoundary(i, myDimBound);
	}

	if (ishierarchized == false)
	{
		OperationHierarchisation* myHierarchisation = myGrid->createOperationHierarchisation();
		myHierarchisation->doHierarchisation(emptyAlpha);
		delete myHierarchisation;
	}

	bGridConstructed = true;

	delete myImporter;
}

void PDESolver::storeGridBonn(std::string tfilename, DataVector& alpha, bool ishierarchized)
{
	IOToolBonnSG* myExporter = new IOToolBonnSG();
	DataVector copyAlpha(alpha);

	if (bGridConstructed)
	{
		if (ishierarchized == false)
		{
			OperationHierarchisation* myHierarchisation = myGrid->createOperationHierarchisation();
			myHierarchisation->doDehierarchisation(copyAlpha);
			delete myHierarchisation;
		}

		myExporter->writeFile(tfilename, *myGrid, copyAlpha, ishierarchized);
	}
	else
	{
		throw new application_exception("PDESolver::storeGrid : A grid wasn't constructed before!");
	}

	delete myExporter;
}

void PDESolver::printGrid(DataVector& alpha, double PointesPerDimension, std::string tfilename)
{
	GridPrinter myPrinter(*this->myGrid);
	myPrinter.printGrid(alpha, tfilename, PointesPerDimension);
}

void PDESolver::printSparseGrid(DataVector& alpha, std::string tfilename, bool bSurplus)
{
	GridPrinter myPrinter(*this->myGrid);
	myPrinter.printSparseGrid(alpha, tfilename, bSurplus);
}

double PDESolver::evaluatePoint(std::vector<double>& evalPoint, DataVector& alpha)
{
	double result = 0.0;

	if (bGridConstructed)
	{
		OperationEval* myEval = myGrid->createOperationEval();
		result = myEval->eval(alpha, evalPoint);
		delete myEval;
	}
	else
	{
		throw new application_exception("PDESolver::evaluatePoint : A grid wasn't constructed before!");
	}

	return result;
}

void PDESolver::getCuboidEvalPoints(std::vector<DataVector>& evalPoints, DataVector& curPoint, std::vector<double>& center, double size, size_t points, size_t curDim)
{
	if (curDim == 0)
	{
		if (points > 1)
		{
			for (size_t i = 0; i < points; i++)
			{
				curPoint.set(curDim, std::min(std::max(center[curDim]-(myBoundingBox->getIntervalWidth(curDim)*size)+
						((myBoundingBox->getIntervalWidth(curDim)*size*2/(points-1))*static_cast<double>(i)),
						myBoundingBox->getIntervalOffset(curDim)),
						myBoundingBox->getIntervalOffset(curDim)+myBoundingBox->getIntervalWidth(curDim)));

				evalPoints.push_back(curPoint);
			}
		}
		else
		{
			curPoint.set(curDim, center[curDim]);

			evalPoints.push_back(curPoint);
		}
	}
	else
	{
		if (points > 1)
		{
			for (size_t i = 0; i < points; i++)
			{
				curPoint.set(curDim, std::min(std::max(center[curDim]-(myBoundingBox->getIntervalWidth(curDim)*size)+
						((myBoundingBox->getIntervalWidth(curDim)*size*2/(points-1))*static_cast<double>(i)),
						myBoundingBox->getIntervalOffset(curDim)),
						myBoundingBox->getIntervalOffset(curDim)+myBoundingBox->getIntervalWidth(curDim)));

				getCuboidEvalPoints(evalPoints, curPoint, center, size, points, curDim-1);
			}
		}
		else
		{
			curPoint.set(curDim, center[curDim]);

			getCuboidEvalPoints(evalPoints, curPoint, center, size, points, curDim-1);
		}
	}
}

void PDESolver::evaluateCuboid(DataVector& alpha, DataVector& OptionPrices, DataVector& EvaluationPoints)
{
	if (bGridConstructed)
	{
		if (OptionPrices.getSize() != EvaluationPoints.getSize())
		{
			throw new application_exception("PDESolver::evaluateCuboid : The size of the price vector doesn't match the size of the evaluation points' vector!");
		}

		OperationB* myOpB = myGrid->createOperationB();
		myOpB->multTranspose(alpha, EvaluationPoints, OptionPrices);
		delete myOpB;
	}
	else
	{
		throw new application_exception("PDESolver::evaluateCuboid : A grid wasn't constructed before!");
	}
}

void PDESolver::getEvaluationCuboid(DataVector& EvaluationPoints, std::vector<double>& center, double size, size_t points)
{
	std::vector<DataVector> evalPoints;
	DataVector curPoint(getNumberDimensions());

	getCuboidEvalPoints(evalPoints, curPoint, center, size, points, getNumberDimensions()-1);

	size_t numEvalPoints = evalPoints.size();
	EvaluationPoints.resize(numEvalPoints);

	for (size_t i = 0; i < numEvalPoints; i++)
	{
		EvaluationPoints.setRow(i, evalPoints[i]);
	}
}

size_t PDESolver::getNumberGridPoints()
{
	if (bGridConstructed)
	{
		return myGridStorage->size();
	}
	else
	{
		throw new application_exception("PDESolver::getNumberGridPoints : A grid wasn't constructed before!");
	}
}

size_t PDESolver::getNumberInnerGridPoints()
{
	if (bGridConstructed)
	{
		return myGridStorage->getNumInnerPoints();
	}
	else
	{
		throw new application_exception("PDESolver::getNumberGridPoints : A grid wasn't constructed before!");
	}
}

size_t PDESolver:: getNumberDimensions()
{
	if (bGridConstructed)
	{
		return myGridStorage->dim();
	}
	else
	{
		throw new application_exception("PDESolver::getNumberDimensions : A grid wasn't constructed before!");
	}
}

}
