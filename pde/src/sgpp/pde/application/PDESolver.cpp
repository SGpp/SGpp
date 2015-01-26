/******************************************************************************
 * Copyright (C) 2009 Technische Universitaet Muenchen                         *
 * This file is part of the SG++ project. For conditions of distribution and   *
 * use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
 ******************************************************************************/
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de)
#include <sgpp/pde/application/PDESolver.hpp>
#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/exception/application_exception.hpp>
#include <sgpp/base/tools/StdNormalDistribution.hpp>
#include <sgpp/base/grid/generation/functors/SurplusRefinementFunctor.hpp>
#include <sgpp/base/grid/generation/functors/SurplusCoarseningFunctor.hpp>
#include <sgpp/base/operation/BaseOpFactory.hpp>

#include <sstream>
#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <cmath>

namespace sg {
namespace pde {

PDESolver::PDESolver(): levels(0), dim(0), myBoundingBox(nullptr), myGridStorage(nullptr), myGrid(nullptr) {
	//initializers may be wrong - David
	bGridConstructed = false;
}

PDESolver::~PDESolver() {
	if (bGridConstructed) {
		delete myGrid;
	}
}

void PDESolver::getGridNormalDistribution(sg::base::DataVector& alpha, std::vector<double>& norm_mu,
		std::vector<double>& norm_sigma) {
	if (bGridConstructed) {
		double tmp;
		double value;
		sg::base::StdNormalDistribution myNormDistr;

		for (size_t i = 0; i < this->myGrid->getStorage()->size(); i++) {
			std::string coords = this->myGridStorage->get(i)->getCoordsStringBB(*(this->myBoundingBox));
			std::stringstream coordsStream(coords);

			value = 1.0;

			for (size_t j = 0; j < this->dim; j++) {
				coordsStream >> tmp;

				value *= myNormDistr.getDensity(tmp, norm_mu[j], norm_sigma[j]);
			}

			alpha[i] = value;
		}
	} else {
		throw new sg::base::application_exception(
				"PDESolver::getGridNormalDistribution : The grid wasn't initialized before!");
	}
}

void PDESolver::deleteGrid() {
	if (bGridConstructed) {
		delete myGrid;
		bGridConstructed = false;
		myBoundingBox = NULL;
		myGridStorage = NULL;
	} else {
		throw new sg::base::application_exception("PDESolver::deleteGrid : The grid wasn't initialized before!");
	}
}

void PDESolver::setGrid(const std::string& serializedGrid) {
	if (bGridConstructed) {
		delete myGrid;
		bGridConstructed = false;
		myBoundingBox = NULL;
		myGridStorage = NULL;
	}

	myGrid = sg::base::Grid::unserialize(serializedGrid);

	myBoundingBox = myGrid->getBoundingBox();
	myGridStorage = myGrid->getStorage();

	dim = myGrid->getStorage()->dim();
	levels = 0;

	bGridConstructed = true;
}

std::string PDESolver::getGrid() const {
	std::string gridSer = "";

	if (bGridConstructed) {
		// Serialize the grid
		myGrid->serialize(gridSer);
	} else {
		throw new sg::base::application_exception("PDESolver::getGrid : The grid wasn't initialized before!");
	}

	return gridSer;
}

void PDESolver::refineInitialGridSurplus(sg::base::DataVector& alpha, int numRefinePoints, double dThreshold) {
	size_t nRefinements;

	if (numRefinePoints < 0) {
		nRefinements = myGrid->createGridGenerator()->getNumberOfRefinablePoints();
	} else {
		nRefinements = numRefinePoints;
	}

	if (bGridConstructed) {
		sg::base::SurplusRefinementFunctor* myRefineFunc = new sg::base::SurplusRefinementFunctor(&alpha, nRefinements,
				dThreshold);

		myGrid->createGridGenerator()->refine(myRefineFunc);

		delete myRefineFunc;

		alpha.resize(myGridStorage->size());
	} else {
		throw new sg::base::application_exception(
				"PDESolver::refineIntialGridSurplus : The grid wasn't initialized before!");
	}
}

void PDESolver::refineInitialGridSurplusSubDomain(sg::base::DataVector& alpha, int numRefinePoints, double dThreshold,
		std::vector<double>& norm_mu, std::vector<double>& norm_sigma) {
	size_t nRefinements;

	if (numRefinePoints < 0) {
		nRefinements = myGrid->createGridGenerator()->getNumberOfRefinablePoints();
	} else {
		nRefinements = numRefinePoints;
	}

	if (bGridConstructed) {
		sg::base::DataVector stdNormDist(alpha.getSize());

		// calculate multidimensional normal distribution and apply to alpha on it
		this->getGridNormalDistribution(stdNormDist, norm_mu, norm_sigma);
		//printSparseGrid(stdNormDist, "normalDistribution.grid.gnuplot", true);
		stdNormDist.componentwise_mult(alpha);
		//printSparseGrid(stdNormDist, "normalDistribution_refine.grid.gnuplot", true);

		sg::base::SurplusRefinementFunctor* myRefineFunc = new sg::base::SurplusRefinementFunctor(&stdNormDist,
				nRefinements, dThreshold);

		myGrid->createGridGenerator()->refine(myRefineFunc);

		delete myRefineFunc;

		alpha.resize(myGridStorage->size());
	} else {
		throw new sg::base::application_exception(
				"PDESolver::refineIntialGridSurplusSubDomain : The grid wasn't initialized before!");
	}
}

void PDESolver::refineInitialGridSurplusToMaxLevel(sg::base::DataVector& alpha, double dThreshold,
		sg::base::GridStorage::index_type::level_type maxLevel) {
	if (bGridConstructed) {
		size_t nRefinements = myGrid->createGridGenerator()->getNumberOfRefinablePointsToMaxLevel(maxLevel);

		sg::base::SurplusRefinementFunctor* myRefineFunc = new sg::base::SurplusRefinementFunctor(&alpha, nRefinements,
				dThreshold);

		myGrid->createGridGenerator()->refineMaxLevel(myRefineFunc, maxLevel);

		delete myRefineFunc;

		alpha.resize(myGridStorage->size());
	} else {
		throw new sg::base::application_exception(
				"PDESolver::refineInitialGridSurplusToMaxLevel : The grid wasn't initialized before!");
	}
}

void PDESolver::refineInitialGridSurplusToMaxLevelSubDomain(sg::base::DataVector& alpha, double dThreshold,
		sg::base::GridStorage::index_type::level_type maxLevel, std::vector<double>& norm_mu,
		std::vector<double>& norm_sigma) {
	if (bGridConstructed) {
		size_t nRefinements = myGrid->createGridGenerator()->getNumberOfRefinablePointsToMaxLevel(maxLevel);

		sg::base::DataVector stdNormDist(alpha.getSize());

		// calculate multidimensional normal distribution and apply to alpha on it
		this->getGridNormalDistribution(stdNormDist, norm_mu, norm_sigma);
		//printSparseGrid(stdNormDist, "normalDistribution.grid.gnuplot", true);
		stdNormDist.componentwise_mult(alpha);
		//printSparseGrid(stdNormDist, "normalDistribution_refine.grid.gnuplot", true);

		sg::base::SurplusRefinementFunctor* myRefineFunc = new sg::base::SurplusRefinementFunctor(&stdNormDist,
				nRefinements, dThreshold);

		myGrid->createGridGenerator()->refineMaxLevel(myRefineFunc, maxLevel);

		delete myRefineFunc;

		alpha.resize(myGridStorage->size());
	} else {
		throw new sg::base::application_exception(
				"PDESolver::refineInitialGridSurplusToMaxLevelSubDomain : The grid wasn't initialized before!");
	}
}

void PDESolver::coarsenInitialGridSurplus(sg::base::DataVector& alpha, double dThreshold) {
	if (bGridConstructed) {
		sg::base::GridGenerator* myGenerator = myGrid->createGridGenerator();
		size_t numCoarsen = myGenerator->getNumberOfRemovablePoints();
		size_t originalGridSize = myGrid->getStorage()->size();
		sg::base::SurplusCoarseningFunctor* myCoarsenFunctor = new sg::base::SurplusCoarseningFunctor(&alpha,
				numCoarsen, dThreshold);

		myGenerator->coarsenNFirstOnly(myCoarsenFunctor, &alpha, originalGridSize);

		delete myCoarsenFunctor;
		delete myGenerator;
	} else {
		throw new sg::base::application_exception(
				"PDESolver::coarsenInitialGridSurplus : The grid wasn't initialized before!");
	}
}

void PDESolver::printLevelIndexGrid(std::string tfilename) const {
	sg::base::GridPrinter myPrinter(*this->myGrid);
	myPrinter.printLevelIndexGrid(tfilename);
}

void PDESolver::printGrid(sg::base::DataVector& alpha, double PointesPerDimension, std::string tfilename) const {
	sg::base::GridPrinter myPrinter(*this->myGrid);
	myPrinter.printGrid(alpha, tfilename, static_cast<size_t>(PointesPerDimension));
}

void PDESolver::printGridDomain(sg::base::DataVector& alpha, double PointesPerDimension,
		sg::base::BoundingBox& GridArea, std::string tfilename) const {
	sg::base::GridPrinter myPrinter(*this->myGrid);
	myPrinter.printGridDomain(alpha, tfilename, GridArea, static_cast<size_t>(PointesPerDimension));
}

void PDESolver::printSparseGrid(sg::base::DataVector& alpha, std::string tfilename, bool bSurplus) const {
	sg::base::GridPrinter myPrinter(*this->myGrid);
	myPrinter.printSparseGrid(alpha, tfilename, bSurplus);
}

void PDESolver::printSparseGridExpTransform(sg::base::DataVector& alpha, std::string tfilename, bool bSurplus) const {
	sg::base::GridPrinter myPrinter(*this->myGrid);
	myPrinter.printSparseGridExpTransform(alpha, tfilename, bSurplus);
}

double PDESolver::evaluatePoint(std::vector<double>& evalPoint, sg::base::DataVector& alpha) {
	double result = 0.0;

	if (bGridConstructed) {
		sg::base::OperationEval* myEval = sg::op_factory::createOperationEval(*myGrid);
		result = myEval->eval(alpha, evalPoint);
		delete myEval;
	} else {
		throw new sg::base::application_exception("PDESolver::evaluatePoint : A grid wasn't constructed before!");
	}

	return result;
}

void PDESolver::evaluateCuboid(sg::base::DataVector& alpha, sg::base::DataVector& OptionPrices,
		sg::base::DataMatrix& EvaluationPoints) {
	if (bGridConstructed) {
		if (OptionPrices.getSize() != EvaluationPoints.getNrows()) {
			throw new sg::base::application_exception(
					"PDESolver::evaluateCuboid : The size of the price vector doesn't match the size of the evaluation points' vector!");
		}

		sg::base::OperationMultipleEval* myOpMultEval = sg::op_factory::createOperationMultipleEval(*myGrid,
				EvaluationPoints);
		myOpMultEval->mult(alpha, OptionPrices);
		delete myOpMultEval;
	} else {
		throw new sg::base::application_exception("PDESolver::evaluateCuboid : A grid wasn't constructed before!");
	}
}

size_t PDESolver::getNumberGridPoints() const {
	if (bGridConstructed) {
		return myGridStorage->size();
	} else {
		throw new sg::base::application_exception("PDESolver::getNumberGridPoints : A grid wasn't constructed before!");
	}
}

size_t PDESolver::getNumberInnerGridPoints() const {
	if (bGridConstructed) {
		return myGridStorage->getNumInnerPoints();
	} else {
		throw new sg::base::application_exception("PDESolver::getNumberGridPoints : A grid wasn't constructed before!");
	}
}

size_t PDESolver::getNumberDimensions() const {
	if (bGridConstructed) {
		return myGridStorage->dim();
	} else {
		throw new sg::base::application_exception("PDESolver::getNumberDimensions : A grid wasn't constructed before!");
	}
}

}
}
