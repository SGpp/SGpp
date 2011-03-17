/******************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)

#include "application/pde/PDESolver.hpp"
#include "grid/Grid.hpp"
#include "exception/application_exception.hpp"
#include "tools/common/StdNormalDistribution.hpp"
#include "grid/generation/SurplusRefinementFunctor.hpp"
#include "grid/generation/SurplusCoarseningFunctor.hpp"

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

void PDESolver::getGridNormalDistribution(DataVector& alpha, std::vector<double>& norm_mu, std::vector<double>& norm_sigma)
{
	if (bGridConstructed)
	{
		double tmp;
		double value;
		StdNormalDistribution myNormDistr;

		for (size_t i = 0; i < this->myGrid->getStorage()->size(); i++)
		{
			std::string coords = this->myGridStorage->get(i)->getCoordsStringBB(*(this->myBoundingBox));
			std::stringstream coordsStream(coords);

			value = 1.0;
			for (size_t j = 0; j < this->dim; j++)
			{
				coordsStream >> tmp;

				value *= myNormDistr.getDensity(tmp, norm_mu[j], norm_sigma[j]);
			}

			alpha[i] = value;
		}
	}
	else
	{
		throw new application_exception("PDESolver::getGridNormalDistribution : The grid wasn't initialized before!");
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

std::string PDESolver::getGrid() const
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

void PDESolver::refineInitialGridSurplus(DataVector& alpha, int numRefinePoints, double dThreshold)
{
	size_t nRefinements;

	if (numRefinePoints < 0)
	{
		nRefinements = myGrid->createGridGenerator()->getNumberOfRefinablePoints();
	}
	else
	{
		nRefinements = numRefinePoints;
	}

	if (bGridConstructed)
	{
		SurplusRefinementFunctor* myRefineFunc = new SurplusRefinementFunctor(&alpha, nRefinements, dThreshold);

		myGrid->createGridGenerator()->refine(myRefineFunc);

		delete myRefineFunc;

		alpha.resize(myGridStorage->size());
	}
	else
	{
		throw new application_exception("PDESolver::refineIntialGridSurplus : The grid wasn't initialized before!");
	}
}

void PDESolver::refineInitialGridSurplusSubDomain(DataVector& alpha, int numRefinePoints, double dThreshold, std::vector<double>& norm_mu, std::vector<double>& norm_sigma)
{
	size_t nRefinements;

	if (numRefinePoints < 0)
	{
		nRefinements = myGrid->createGridGenerator()->getNumberOfRefinablePoints();
	}
	else
	{
		nRefinements = numRefinePoints;
	}

	if (bGridConstructed)
	{
		DataVector stdNormDist(alpha.getSize());

		// calculate multidimensional normal distribution and apply to alpha on it
		this->getGridNormalDistribution(stdNormDist, norm_mu, norm_sigma);
		//printSparseGrid(stdNormDist, "normalDistribution.grid.gnuplot", true);
		stdNormDist.componentwise_mult(alpha);
		//printSparseGrid(stdNormDist, "normalDistribution_refine.grid.gnuplot", true);

		SurplusRefinementFunctor* myRefineFunc = new SurplusRefinementFunctor(&stdNormDist, nRefinements, dThreshold);

		myGrid->createGridGenerator()->refine(myRefineFunc);

		delete myRefineFunc;

		alpha.resize(myGridStorage->size());
	}
	else
	{
		throw new application_exception("PDESolver::refineIntialGridSurplusSubDomain : The grid wasn't initialized before!");
	}
}

void PDESolver::refineInitialGridSurplusToMaxLevel(DataVector& alpha, double dThreshold, size_t maxLevel)
{
	if (bGridConstructed)
	{
		size_t nRefinements = myGrid->createGridGenerator()->getNumberOfRefinablePointsToMaxLevel(maxLevel);

		SurplusRefinementFunctor* myRefineFunc = new SurplusRefinementFunctor(&alpha, nRefinements, dThreshold);

		myGrid->createGridGenerator()->refineMaxLevel(myRefineFunc, maxLevel);

		delete myRefineFunc;

		alpha.resize(myGridStorage->size());
	}
	else
	{
		throw new application_exception("PDESolver::refineInitialGridSurplusToMaxLevel : The grid wasn't initialized before!");
	}
}

void PDESolver::refineInitialGridSurplusToMaxLevelSubDomain(DataVector& alpha, double dThreshold, size_t maxLevel, std::vector<double>& norm_mu, std::vector<double>& norm_sigma)
{
	if (bGridConstructed)
	{
		size_t nRefinements = myGrid->createGridGenerator()->getNumberOfRefinablePointsToMaxLevel(maxLevel);

		DataVector stdNormDist(alpha.getSize());

		// calculate multidimensional normal distribution and apply to alpha on it
		this->getGridNormalDistribution(stdNormDist, norm_mu, norm_sigma);
		//printSparseGrid(stdNormDist, "normalDistribution.grid.gnuplot", true);
		stdNormDist.componentwise_mult(alpha);
		//printSparseGrid(stdNormDist, "normalDistribution_refine.grid.gnuplot", true);

		SurplusRefinementFunctor* myRefineFunc = new SurplusRefinementFunctor(&stdNormDist, nRefinements, dThreshold);

		myGrid->createGridGenerator()->refineMaxLevel(myRefineFunc, maxLevel);

		delete myRefineFunc;

		alpha.resize(myGridStorage->size());
	}
	else
	{
		throw new application_exception("PDESolver::refineInitialGridSurplusToMaxLevelSubDomain : The grid wasn't initialized before!");
	}
}

void PDESolver::coarsenInitialGridSurplus(DataVector& alpha, double dThreshold)
{
	if (bGridConstructed)
	{
		GridGenerator* myGenerator =  myGrid->createGridGenerator();
		size_t numCoarsen = myGenerator->getNumberOfRemoveablePoints();
		size_t originalGridSize = myGrid->getStorage()->size();
		SurplusCoarseningFunctor* myCoarsenFunctor = new SurplusCoarseningFunctor(&alpha, numCoarsen, dThreshold);

		myGenerator->coarsenNFirstOnly(myCoarsenFunctor, &alpha, originalGridSize);

		delete myCoarsenFunctor;
		delete myGenerator;
	}
	else
	{
		throw new application_exception("PDESolver::coarsenInitialGridSurplus : The grid wasn't initialized before!");
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

void PDESolver::printGrid(DataVector& alpha, double PointesPerDimension, std::string tfilename) const
{
	GridPrinter myPrinter(*this->myGrid);
	myPrinter.printGrid(alpha, tfilename, PointesPerDimension);
}

void PDESolver::printGridDomain(DataVector& alpha, double PointesPerDimension, BoundingBox& GridArea, std::string tfilename) const
{
	GridPrinter myPrinter(*this->myGrid);
	myPrinter.printGridDomain(alpha, tfilename, GridArea, PointesPerDimension);
}

void PDESolver::printSparseGrid(DataVector& alpha, std::string tfilename, bool bSurplus) const
{
	GridPrinter myPrinter(*this->myGrid);
	myPrinter.printSparseGrid(alpha, tfilename, bSurplus);
}

void PDESolver::printSparseGridExpTransform(DataVector& alpha, std::string tfilename, bool bSurplus) const
{
	GridPrinter myPrinter(*this->myGrid);
	myPrinter.printSparseGridExpTransform(alpha, tfilename, bSurplus);
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

void PDESolver::evaluateCuboid(DataVector& alpha, DataVector& OptionPrices, DataMatrix& EvaluationPoints)
{
	if (bGridConstructed)
	{
		if (OptionPrices.getSize() != EvaluationPoints.getNrows())
		{
			throw new application_exception("PDESolver::evaluateCuboid : The size of the price vector doesn't match the size of the evaluation points' vector!");
		}

		OperationMultipleEval* myOpMultEval = myGrid->createOperationMultipleEval(&EvaluationPoints);
		myOpMultEval->mult(alpha, OptionPrices);
		delete myOpMultEval;
	}
	else
	{
		throw new application_exception("PDESolver::evaluateCuboid : A grid wasn't constructed before!");
	}
}

size_t PDESolver::getNumberGridPoints() const
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

size_t PDESolver::getNumberInnerGridPoints() const
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

size_t PDESolver::getNumberDimensions() const
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
