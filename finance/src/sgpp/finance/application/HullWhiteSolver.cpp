// Copyright (C) 2008-today The SG++ Project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/finance/algorithm/HullWhiteParabolicPDESolverSystem.hpp>
#include <sgpp/finance/application/HullWhiteSolver.hpp>
#include <sgpp/solver/ode/Euler.hpp>
#include <sgpp/solver/ode/CrankNicolson.hpp>
#include <sgpp/solver/sle/BiCGStab.hpp>
#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/exception/application_exception.hpp>
#include <sgpp/base/operation/BaseOpFactory.hpp>

#include <sgpp/globaldef.hpp>
#include <sgpp/base/tools/SGppStopwatch.hpp>

#include <cstdlib>
#include <sstream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <string>
#include <vector>

namespace sgpp {
namespace finance {

HullWhiteSolver::HullWhiteSolver() : ParabolicPDESolver() {
  this->bStochasticDataAlloc = false;
  this->bGridConstructed = false;
  this->myScreen = NULL;
  this->useCoarsen = false;
  this->coarsenThreshold = 0.0;
  this->adaptSolveMode = "none";
  this->refineMode = "classic";
  this->numCoarsenPoints = -1;
  this->refineMaxLevel = 0;

  this->a = 0.0;
  this->refineThreshold = 0.0;
  this->sigma = 0.0;
  this->theta = 0.0;
}

HullWhiteSolver::~HullWhiteSolver() {
  /*if (this->bStochasticDataAlloc)
  {
    delete this->mus;
    delete this->sigmas;
    delete this->rhos;
  }*/
  if (this->myScreen != NULL) {
    delete this->myScreen;
  }
}

void HullWhiteSolver::constructGrid(base::BoundingBox& BoundingBox, size_t level) {
  this->dim = BoundingBox.getDimensions();
  this->levels = level;

  this->myGrid = new base::LinearBoundaryGrid(BoundingBox);

  this->myGrid->getGenerator().regular(this->levels);

  this->myBoundingBox = &this->myGrid->getBoundingBox();
  this->myGridStorage = &this->myGrid->getStorage();

  // std::string serGrid;
  // myGrid->serialize(serGrid);
  // std::cout << serGrid << std::endl;

  this->bGridConstructed = true;
}

void HullWhiteSolver::setStochasticData(double theta, double sigma, double a) {
  this->theta = theta;
  this->sigma = sigma;
  this->a = a;

  bStochasticDataAlloc = true;
}

void HullWhiteSolver::solveExplicitEuler(size_t numTimesteps, double timestepsize,
                                         size_t maxCGIterations, double epsilonCG,
                                         base::DataVector& alpha, bool verbose,
                                         bool generateAnimation) {
  if (this->bGridConstructed && this->bStochasticDataAlloc) {
    solver::Euler* myEuler = new solver::Euler("ExEul", numTimesteps, timestepsize,
                                               generateAnimation, myScreen);
    solver::BiCGStab* myCG = new solver::BiCGStab(maxCGIterations, epsilonCG);
    HullWhiteParabolicPDESolverSystem* myHWSystem = new HullWhiteParabolicPDESolverSystem(
        *this->myGrid, alpha, this->sigma, this->theta, this->a, timestepsize, "ExEul",
        this->useCoarsen, this->coarsenThreshold, this->adaptSolveMode, this->numCoarsenPoints,
        this->refineThreshold, this->refineMode, this->refineMaxLevel);
    base::SGppStopwatch* myStopwatch = new base::SGppStopwatch();
    double execTime;

    std::cout << "Using Explicit Euler to solve " << numTimesteps << " timesteps:" << std::endl;
    myStopwatch->start();
    myEuler->solve(*myCG, *myHWSystem, true, verbose);
    execTime = myStopwatch->stop();

    std::cout << std::endl << "Final Grid size: " << getNumberGridPoints() << std::endl;
    std::cout << "Final Grid size (inner): " << getNumberInnerGridPoints() << std::endl
              << std::endl
              << std::endl;

    std::cout << "Average Grid size: "
              << static_cast<double>(myHWSystem->getSumGridPointsComplete()) /
                     static_cast<double>(numTimesteps)
              << std::endl;
    std::cout << "Average Grid size (Inner): "
              << static_cast<double>(myHWSystem->getSumGridPointsInner()) /
                     static_cast<double>(numTimesteps)
              << std::endl
              << std::endl
              << std::endl;

    if (this->myScreen != NULL) {
      std::cout << "Time to solve: " << execTime << " seconds" << std::endl;
      this->myScreen->writeEmptyLines(2);
    }

    delete myHWSystem;
    delete myCG;
    delete myEuler;
    delete myStopwatch;
  } else {
    throw base::application_exception(
        "HullWhiteSolver::solveExplicitEuler : A grid wasn't constructed before or stochastic "
        "parameters weren't set!");
  }
}

void HullWhiteSolver::solveImplicitEuler(size_t numTimesteps, double timestepsize,
                                         size_t maxCGIterations, double epsilonCG,
                                         base::DataVector& alpha, bool verbose,
                                         bool generateAnimation) {
  if (this->bGridConstructed && this->bStochasticDataAlloc) {
    solver::Euler* myEuler = new solver::Euler("ImEul", numTimesteps, timestepsize,
                                               generateAnimation, myScreen);
    solver::BiCGStab* myCG = new solver::BiCGStab(maxCGIterations, epsilonCG);
    HullWhiteParabolicPDESolverSystem* myHWSystem = new HullWhiteParabolicPDESolverSystem(
        *this->myGrid, alpha, this->sigma, this->theta, this->a, timestepsize, "ImEul",
        this->useCoarsen, this->coarsenThreshold, this->adaptSolveMode, this->numCoarsenPoints,
        this->refineThreshold, this->refineMode, this->refineMaxLevel);
    base::SGppStopwatch* myStopwatch = new base::SGppStopwatch();
    double execTime;

    std::cout << "Using Implicit Euler to solve " << numTimesteps << " timesteps:" << std::endl;
    myStopwatch->start();
    myEuler->solve(*myCG, *myHWSystem, true, verbose);
    execTime = myStopwatch->stop();

    std::cout << std::endl << "Final Grid size: " << getNumberGridPoints() << std::endl;
    std::cout << "Final Grid size (inner): " << getNumberInnerGridPoints() << std::endl
              << std::endl
              << std::endl;

    std::cout << "Average Grid size: "
              << static_cast<double>(myHWSystem->getSumGridPointsComplete()) /
                     static_cast<double>(numTimesteps)
              << std::endl;
    std::cout << "Average Grid size (Inner): "
              << static_cast<double>(myHWSystem->getSumGridPointsInner()) /
                     static_cast<double>(numTimesteps)
              << std::endl
              << std::endl
              << std::endl;

    if (this->myScreen != NULL) {
      std::cout << "Time to solve: " << execTime << " seconds" << std::endl;
      this->myScreen->writeEmptyLines(2);
    }

    delete myHWSystem;
    delete myCG;
    delete myEuler;
    delete myStopwatch;
  } else {
    throw base::application_exception(
        "HullWhiteSolver::solveImplicitEuler : A grid wasn't constructed before or stochastic "
        "parameters weren't set!");
  }
}

void HullWhiteSolver::solveCrankNicolson(size_t numTimesteps, double timestepsize,
                                         size_t maxCGIterations, double epsilonCG,
                                         base::DataVector& alpha, size_t NumImEul) {
  throw base::application_exception(
      "HullWhiteSolver::solveCrankNicolson : Crank-Nicolson is not supported for "
      "HullWhiteSolver!!");
}

void HullWhiteSolver::initGridWithPayoff(base::DataVector& alpha, double strike,
                                         std::string payoffType, double sigma, double a,
                                         double t, double T) {
  double tmp;

  if (this->bGridConstructed) {
    for (size_t i = 0; i < this->myGrid->getSize(); i++) {
      std::string coords = this->myGridStorage->get(i)->getCoordsStringBB(*this->myBoundingBox);
      std::stringstream coordsStream(coords);
      coordsStream >> tmp;
      double* dblFuncValues = new double[1];

      for (size_t j = 0; j < 1; j++) {
        coordsStream >> tmp;

        dblFuncValues[j] = exp((0.04 * (t - T)) + 0.04 * (1 - exp(-a * (T - t))) / a -
                               1 / (4 * pow(a, 3)) * pow(sigma, 2) *
                                   pow((exp(-a * T) - exp(-a * t)), 2) * (exp(2 * a * t) - 1) -
                               tmp / a * (1 - exp(-a * (T - t))));
      }

      if (payoffType == "std_euro_call") {
        tmp = 0.0;

        for (size_t j = 0; j < 1; j++) {
          tmp += dblFuncValues[j];
        }

        alpha[i] = std::max<double>(((tmp)-strike), 0.0);
      } else if (payoffType == "std_euro_put") {
        tmp = 0.0;

        for (size_t j = 0; j < dim; j++) {
          tmp += dblFuncValues[j];
        }

        alpha[i] = std::max<double>(strike - ((tmp)), 0.0);
      } else {
        throw base::application_exception(
            "HullWhiteSolver::initGridWithPayoff : An unknown payoff-type was specified!");
      }

      delete[] dblFuncValues;

      // delete dblFuncValues;
    }

    op_factory::createOperationHierarchisation(*this->myGrid)->doHierarchisation(alpha);
  } else {
    throw base::application_exception(
        "HullWhiteSolver::initGridWithPayoff : A grid wasn't constructed before!");
  }
}

std::vector<size_t> HullWhiteSolver::getAlgorithmicDimensions() {
  return this->myGrid->getAlgorithmicDimensions();
}

void HullWhiteSolver::setAlgorithmicDimensions(std::vector<size_t> newAlgoDims) {
  this->myGrid->setAlgorithmicDimensions(newAlgoDims);
}

void HullWhiteSolver::initScreen() {
  this->myScreen = new base::ScreenOutput();
  this->myScreen->writeTitle("SGpp - Hull White Solver, 1.3.0",
                             "The SG++ Project (C) 2009-2010, by Chao qi");
  this->myScreen->writeStartSolve("One dimensional Hull White Solver");
}
}  // namespace finance
}  // namespace sgpp
