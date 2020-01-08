// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/exception/application_exception.hpp>
#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/grid/generation/functors/SurplusCoarseningFunctor.hpp>
#include <sgpp/base/grid/generation/functors/SurplusRefinementFunctor.hpp>
#include <sgpp/base/operation/BaseOpFactory.hpp>
#include <sgpp/base/tools/StdNormalDistribution.hpp>
#include <sgpp/pde/application/PDESolver.hpp>

#include <sgpp/globaldef.hpp>

#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace sgpp {
namespace pde {

PDESolver::PDESolver()
    : levels(0), dim(0), myBoundingBox(nullptr), myGridStorage(nullptr), myGrid(nullptr) {
  // initializers may be wrong - David
  bGridConstructed = false;
}

PDESolver::~PDESolver() {
  if (bGridConstructed) {
    delete myGrid;
  }
}

void PDESolver::getGridNormalDistribution(sgpp::base::DataVector& alpha,
                                          std::vector<double>& norm_mu,
                                          std::vector<double>& norm_sigma) {
  if (bGridConstructed) {
    double tmp;
    double value;
    sgpp::base::StdNormalDistribution myNormDistr;

    for (size_t i = 0; i < this->myGrid->getSize(); i++) {
      std::string coords =
          this->myGridStorage->getCoordinates(this->myGridStorage->getPoint(i)).toString();
      std::stringstream coordsStream(coords);

      value = 1.0;

      for (size_t j = 0; j < this->dim; j++) {
        coordsStream >> tmp;

        value *= myNormDistr.getDensity(tmp, norm_mu[j], norm_sigma[j]);
      }

      alpha[i] = value;
    }
  } else {
    throw sgpp::base::application_exception(
        "PDESolver::getGridNormalDistribution : The grid wasn't initialized before!");
  }
}

void PDESolver::deleteGrid() {
  if (bGridConstructed) {
    delete myGrid;
    bGridConstructed = false;
    myBoundingBox = nullptr;
    myGridStorage = nullptr;
  } else {
    throw sgpp::base::application_exception(
        "PDESolver::deleteGrid : The grid wasn't initialized before!");
  }
}

void PDESolver::setGrid(const std::string& serializedGrid) {
  if (bGridConstructed) {
    delete myGrid;
    bGridConstructed = false;
    myBoundingBox = nullptr;
    myGridStorage = nullptr;
  }

  myGrid = sgpp::base::Grid::unserialize(serializedGrid);

  myBoundingBox = &myGrid->getBoundingBox();
  myGridStorage = &myGrid->getStorage();

  dim = myGrid->getDimension();
  levels = 0;

  bGridConstructed = true;
}

std::string PDESolver::getGrid() const {
  std::string gridSer = "";

  if (bGridConstructed) {
    // Serialize the grid
    myGrid->serialize(gridSer);
  } else {
    throw sgpp::base::application_exception(
        "PDESolver::getGrid : The grid wasn't initialized before!");
  }

  return gridSer;
}

void PDESolver::refineInitialGridSurplus(sgpp::base::DataVector& alpha, int numRefinePoints,
                                         double dThreshold) {
  size_t nRefinements;

  if (numRefinePoints < 0) {
    nRefinements = myGrid->getGenerator().getNumberOfRefinablePoints();
  } else {
    nRefinements = numRefinePoints;
  }

  if (bGridConstructed) {
    sgpp::base::SurplusRefinementFunctor myRefineFunc(alpha, nRefinements, dThreshold);
    myGrid->getGenerator().refine(myRefineFunc);

    alpha.resize(myGridStorage->getSize());
  } else {
    throw sgpp::base::application_exception(
        "PDESolver::refineIntialGridSurplus : The grid wasn't initialized before!");
  }
}

void PDESolver::refineInitialGridSurplusSubDomain(sgpp::base::DataVector& alpha,
                                                  int numRefinePoints, double dThreshold,
                                                  std::vector<double>& norm_mu,
                                                  std::vector<double>& norm_sigma) {
  size_t nRefinements;

  if (numRefinePoints < 0) {
    nRefinements = myGrid->getGenerator().getNumberOfRefinablePoints();
  } else {
    nRefinements = numRefinePoints;
  }

  if (bGridConstructed) {
    sgpp::base::DataVector stdNormDist(alpha.getSize());

    // calculate multidimensional normal distribution and apply to alpha on it
    this->getGridNormalDistribution(stdNormDist, norm_mu, norm_sigma);
    // printSparseGrid(stdNormDist, "normalDistribution.grid.gnuplot", true);
    stdNormDist.componentwise_mult(alpha);
    // printSparseGrid(stdNormDist, "normalDistribution_refine.grid.gnuplot", true);

    sgpp::base::SurplusRefinementFunctor myRefineFunc(stdNormDist, nRefinements, dThreshold);
    myGrid->getGenerator().refine(myRefineFunc);

    alpha.resize(myGridStorage->getSize());
  } else {
    throw sgpp::base::application_exception(
        "PDESolver::refineIntialGridSurplusSubDomain : The grid wasn't initialized before!");
  }
}

void PDESolver::refineInitialGridSurplusToMaxLevel(sgpp::base::DataVector& alpha, double dThreshold,
                                                   sgpp::base::level_t maxLevel) {
  if (bGridConstructed) {
    size_t nRefinements = myGrid->getGenerator().getNumberOfRefinablePointsToMaxLevel(maxLevel);

    sgpp::base::SurplusRefinementFunctor myRefineFunc(alpha, nRefinements, dThreshold);
    myGrid->getGenerator().refineMaxLevel(myRefineFunc, maxLevel);

    alpha.resize(myGridStorage->getSize());
  } else {
    throw sgpp::base::application_exception(
        "PDESolver::refineInitialGridSurplusToMaxLevel : The grid wasn't initialized before!");
  }
}

void PDESolver::refineInitialGridSurplusToMaxLevelSubDomain(sgpp::base::DataVector& alpha,
                                                            double dThreshold,
                                                            sgpp::base::level_t maxLevel,
                                                            std::vector<double>& norm_mu,
                                                            std::vector<double>& norm_sigma) {
  if (bGridConstructed) {
    size_t nRefinements = myGrid->getGenerator().getNumberOfRefinablePointsToMaxLevel(maxLevel);

    sgpp::base::DataVector stdNormDist(alpha.getSize());

    // calculate multidimensional normal distribution and apply to alpha on it
    this->getGridNormalDistribution(stdNormDist, norm_mu, norm_sigma);
    // printSparseGrid(stdNormDist, "normalDistribution.grid.gnuplot", true);
    stdNormDist.componentwise_mult(alpha);
    // printSparseGrid(stdNormDist, "normalDistribution_refine.grid.gnuplot", true);

    sgpp::base::SurplusRefinementFunctor myRefineFunc(stdNormDist, nRefinements, dThreshold);
    myGrid->getGenerator().refineMaxLevel(myRefineFunc, maxLevel);

    alpha.resize(myGridStorage->getSize());
  } else {
    throw sgpp::base::application_exception(
        "PDESolver::refineInitialGridSurplusToMaxLevelSubDomain : The grid wasn't initialized "
        "before!");
  }
}

void PDESolver::coarsenInitialGridSurplus(sgpp::base::DataVector& alpha, double dThreshold) {
  if (bGridConstructed) {
    sgpp::base::GridGenerator& myGenerator = myGrid->getGenerator();
    size_t numCoarsen = myGenerator.getNumberOfRemovablePoints();
    size_t originalGridSize = myGrid->getSize();
    sgpp::base::SurplusCoarseningFunctor myCoarsenFunctor(alpha, numCoarsen, dThreshold);
    myGenerator.coarsenNFirstOnly(myCoarsenFunctor, alpha, originalGridSize, nullptr);
  } else {
    throw sgpp::base::application_exception(
        "PDESolver::coarsenInitialGridSurplus : The grid wasn't initialized before!");
  }
}

void PDESolver::printLevelIndexGrid(std::string tfilename) const {
  sgpp::base::GridPrinter myPrinter(*this->myGrid);
  myPrinter.printLevelIndexGrid(tfilename);
}

void PDESolver::printGrid(sgpp::base::DataVector& alpha, size_t PointesPerDimension,
                          std::string tfilename) const {
  sgpp::base::GridPrinter myPrinter(*this->myGrid);
  myPrinter.printGrid(alpha, tfilename, PointesPerDimension);
}

void PDESolver::printGridDomain(sgpp::base::DataVector& alpha, size_t PointesPerDimension,
                                sgpp::base::BoundingBox& GridArea, std::string tfilename) const {
  sgpp::base::GridPrinter myPrinter(*this->myGrid);
  myPrinter.printGridDomain(alpha, tfilename, GridArea, PointesPerDimension);
}

void PDESolver::printSparseGrid(sgpp::base::DataVector& alpha, std::string tfilename,
                                bool bSurplus) const {
  sgpp::base::GridPrinter myPrinter(*this->myGrid);
  myPrinter.printSparseGrid(alpha, tfilename, bSurplus);
}

void PDESolver::printSparseGridExpTransform(sgpp::base::DataVector& alpha, std::string tfilename,
                                            bool bSurplus) const {
  sgpp::base::GridPrinter myPrinter(*this->myGrid);
  myPrinter.printSparseGridExpTransform(alpha, tfilename, bSurplus);
}

double PDESolver::evaluatePoint(sgpp::base::DataVector& evalPoint, sgpp::base::DataVector& alpha) {
  double result = 0.0;

  if (bGridConstructed) {
    result = sgpp::op_factory::createOperationEval(*myGrid)->eval(alpha, evalPoint);
  } else {
    throw sgpp::base::application_exception(
        "PDESolver::evaluatePoint : A grid wasn't constructed before!");
  }

  return result;
}

void PDESolver::evaluateCuboid(sgpp::base::DataVector& alpha, sgpp::base::DataVector& OptionPrices,
                               sgpp::base::DataMatrix& EvaluationPoints) {
  if (bGridConstructed) {
    if (OptionPrices.getSize() != EvaluationPoints.getNrows()) {
      throw sgpp::base::application_exception(
          "PDESolver::evaluateCuboid : The size of the price vector doesn't match the size of the "
          "evaluation points' vector!");
    }

    sgpp::op_factory::createOperationMultipleEval(*myGrid, EvaluationPoints)
        ->mult(alpha, OptionPrices);
  } else {
    throw sgpp::base::application_exception(
        "PDESolver::evaluateCuboid : A grid wasn't constructed before!");
  }
}

size_t PDESolver::getNumberGridPoints() const {
  if (bGridConstructed) {
    return myGridStorage->getSize();
  } else {
    throw sgpp::base::application_exception(
        "PDESolver::getNumberGridPoints : A grid wasn't constructed before!");
  }
}

size_t PDESolver::getNumberInnerGridPoints() const {
  if (bGridConstructed) {
    return myGridStorage->getNumberOfInnerPoints();
  } else {
    throw sgpp::base::application_exception(
        "PDESolver::getNumberGridPoints : A grid wasn't constructed before!");
  }
}

size_t PDESolver::getNumberDimensions() const {
  if (bGridConstructed) {
    return myGridStorage->getDimension();
  } else {
    throw sgpp::base::application_exception(
        "PDESolver::getNumberDimensions : A grid wasn't constructed before!");
  }
}
}  // namespace pde
}  // namespace sgpp
