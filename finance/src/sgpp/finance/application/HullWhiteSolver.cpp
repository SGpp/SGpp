/* ****************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Chao qi (qic@in.tum.de), Stefanie Schraufstetter 8schraufs@in.tum.de)

#include <sgpp/finance/algorithm/HullWhiteParabolicPDESolverSystem.hpp>
#include <sgpp/finance/application/HullWhiteSolver.hpp>
#include <sgpp/solver/ode/Euler.hpp>
#include <sgpp/solver/ode/CrankNicolson.hpp>
#include <sgpp/solver/sle/BiCGStab.hpp>
#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/exception/application_exception.hpp>
#include <sgpp/base/operation/BaseOpFactory.hpp>
#include <cstdlib>
#include <sstream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <sgpp/base/tools/SGppStopwatch.hpp>

using namespace sg::pde;
using namespace sg::solver;
using namespace sg::base;

namespace sg {
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

      // @todo set to random value to remove compiler warnings due to uninitialized members
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

    void HullWhiteSolver::constructGrid(BoundingBox& BoundingBox, int level) {
      this->dim = BoundingBox.getDimensions();
      this->levels = level;

      this->myGrid = new LinearTrapezoidBoundaryGrid(BoundingBox);

      GridGenerator* myGenerator = this->myGrid->createGridGenerator();
      myGenerator->regular(this->levels);
      delete myGenerator;

      this->myBoundingBox = this->myGrid->getBoundingBox();
      this->myGridStorage = this->myGrid->getStorage();

      //std::string serGrid;
      //myGrid->serialize(serGrid);
      //std::cout << serGrid << std::endl;

      this->bGridConstructed = true;
    }

    void HullWhiteSolver::setStochasticData(double theta, double sigma, double a) {
      this->theta = theta;
      this->sigma = sigma;
      this->a = a;

      bStochasticDataAlloc = true;
    }

    void HullWhiteSolver::solveExplicitEuler(size_t numTimesteps, double timestepsize, size_t maxCGIterations, double epsilonCG, DataVector& alpha, bool verbose, bool generateAnimation, size_t numEvalsAnimation) {
      if (this->bGridConstructed && this->bStochasticDataAlloc) {
        Euler* myEuler = new Euler("ExEul", numTimesteps, timestepsize, generateAnimation, numEvalsAnimation, myScreen);
        BiCGStab* myCG = new BiCGStab(maxCGIterations, epsilonCG);
        HullWhiteParabolicPDESolverSystem* myHWSystem = new HullWhiteParabolicPDESolverSystem(*this->myGrid, alpha, this->sigma, this->theta, this->a, timestepsize, "ExEul", this->useCoarsen, this->coarsenThreshold, this->adaptSolveMode, this->numCoarsenPoints, this->refineThreshold, this->refineMode, this->refineMaxLevel);
        SGppStopwatch* myStopwatch = new SGppStopwatch();
        double execTime;

        std::cout << "Using Explicit Euler to solve " << numTimesteps << " timesteps:" << std::endl;
        myStopwatch->start();
        myEuler->solve(*myCG, *myHWSystem, true, verbose);
        execTime = myStopwatch->stop();

        std::cout << std::endl << "Final Grid size: " << getNumberGridPoints() << std::endl;
        std::cout << "Final Grid size (inner): " << getNumberInnerGridPoints() << std::endl << std::endl << std::endl;

        std::cout << "Average Grid size: " << static_cast<double>(myHWSystem->getSumGridPointsComplete()) / static_cast<double>(numTimesteps) << std::endl;
        std::cout << "Average Grid size (Inner): " << static_cast<double>(myHWSystem->getSumGridPointsInner()) / static_cast<double>(numTimesteps) << std::endl << std::endl << std::endl;

        if (this->myScreen != NULL) {
          std::cout << "Time to solve: " << execTime << " seconds" << std::endl;
          this->myScreen->writeEmptyLines(2);
        }

        delete myHWSystem;
        delete myCG;
        delete myEuler;
        delete myStopwatch;
      } else {
        throw new application_exception("HullWhiteSolver::solveExplicitEuler : A grid wasn't constructed before or stochastic parameters weren't set!");
      }
    }


    void HullWhiteSolver::solveImplicitEuler(size_t numTimesteps, double timestepsize, size_t maxCGIterations, double epsilonCG, DataVector& alpha, bool verbose, bool generateAnimation, size_t numEvalsAnimation) {
      if (this->bGridConstructed && this->bStochasticDataAlloc) {
        Euler* myEuler = new Euler("ImEul", numTimesteps, timestepsize, generateAnimation, numEvalsAnimation, myScreen);
        BiCGStab* myCG = new BiCGStab(maxCGIterations, epsilonCG);
        HullWhiteParabolicPDESolverSystem* myHWSystem = new HullWhiteParabolicPDESolverSystem(*this->myGrid, alpha, this->sigma, this->theta, this->a, timestepsize, "ImEul", this->useCoarsen, this->coarsenThreshold, this->adaptSolveMode, this->numCoarsenPoints, this->refineThreshold, this->refineMode, this->refineMaxLevel);
        SGppStopwatch* myStopwatch = new SGppStopwatch();
        double execTime;

        std::cout << "Using Implicit Euler to solve " << numTimesteps << " timesteps:" << std::endl;
        myStopwatch->start();
        myEuler->solve(*myCG, *myHWSystem, true, verbose);
        execTime = myStopwatch->stop();

        std::cout << std::endl << "Final Grid size: " << getNumberGridPoints() << std::endl;
        std::cout << "Final Grid size (inner): " << getNumberInnerGridPoints() << std::endl << std::endl << std::endl;

        std::cout << "Average Grid size: " << static_cast<double>(myHWSystem->getSumGridPointsComplete()) / static_cast<double>(numTimesteps) << std::endl;
        std::cout << "Average Grid size (Inner): " << static_cast<double>(myHWSystem->getSumGridPointsInner()) / static_cast<double>(numTimesteps) << std::endl << std::endl << std::endl;

        if (this->myScreen != NULL) {
          std::cout << "Time to solve: " << execTime << " seconds" << std::endl;
          this->myScreen->writeEmptyLines(2);
        }

        delete myHWSystem;
        delete myCG;
        delete myEuler;
        delete myStopwatch;
      } else {
        throw new application_exception("HullWhiteSolver::solveImplicitEuler : A grid wasn't constructed before or stochastic parameters weren't set!");
      }
    }

    void HullWhiteSolver::solveCrankNicolson(size_t numTimesteps, double timestepsize, size_t maxCGIterations, double epsilonCG, DataVector& alpha, size_t NumImEul) {
      throw new application_exception("HullWhiteSolver::solveCrankNicolson : Crank-Nicolson is not supported for HullWhiteSolver!!");
    }


    void HullWhiteSolver::initGridWithPayoff(DataVector& alpha, double strike, std::string payoffType, double sigma, double a, double t, double T) {
      double tmp;

      if (this->bGridConstructed) {
        for (size_t i = 0; i < this->myGrid->getStorage()->size(); i++) {
          std::string coords = this->myGridStorage->get(i)->getCoordsStringBB(*this->myBoundingBox);
          std::stringstream coordsStream(coords);
          coordsStream >> tmp;
          double* dblFuncValues = new double[1];

          for (size_t j = 0; j < 1; j++) {
            coordsStream >> tmp;

            dblFuncValues[j] = exp((0.04 * (t - T)) + 0.04 * (1 - exp(-a * (T - t))) / a - 1 / (4 * pow(a, 3)) * pow(sigma, 2) * pow((exp(-a * T) - exp(-a * t)), 2) * (exp(2 * a * t) - 1) - tmp / a * (1 - exp(-a * (T - t))));
          }

          if (payoffType == "std_euro_call") {
            tmp = 0.0;

            for (size_t j = 0; j < 1; j++) {
              tmp += dblFuncValues[j];
            }

            alpha[i] = std::max<double>(((tmp) - strike), 0.0);
          } else if (payoffType == "std_euro_put") {
            tmp = 0.0;

            for (size_t j = 0; j < dim; j++) {
              tmp += dblFuncValues[j];
            }

            alpha[i] = std::max<double>(strike - ((tmp)), 0.0);
          } else {
            throw new application_exception("HullWhiteSolver::initGridWithPayoff : An unknown payoff-type was specified!");
          }

          delete[] dblFuncValues;

          //delete dblFuncValues;
        }

        OperationHierarchisation* myHierarchisation = sg::op_factory::createOperationHierarchisation(*this->myGrid);
        myHierarchisation->doHierarchisation(alpha);
        delete myHierarchisation;
      } else {
        throw new application_exception("HullWhiteSolver::initGridWithPayoff : A grid wasn't constructed before!");
      }
    }

    std::vector<size_t> HullWhiteSolver::getAlgorithmicDimensions() {
      return this->myGrid->getAlgorithmicDimensions();
    }


    void HullWhiteSolver::setAlgorithmicDimensions(std::vector<size_t> newAlgoDims) {
      this->myGrid->setAlgorithmicDimensions(newAlgoDims);
    }

    void HullWhiteSolver::initScreen() {
      this->myScreen = new ScreenOutput();
      this->myScreen->writeTitle("SGpp - Hull White Solver, 1.3.0", "TUM (C) 2009-2010, by Chao qi");
      this->myScreen->writeStartSolve("One dimensional Hull White Solver");
    }
  }
}
