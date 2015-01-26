// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at 
// sgpp.sparsegrids.org

#include <sgpp/finance/algorithm/ModifiedBlackScholesParabolicPDESolverSystem.hpp>
#include <sgpp/finance/algorithm/HullWhiteParabolicPDESolverSystem.hpp>
#include <sgpp/finance/application/BlackScholesSolver.hpp>
#include <sgpp/finance/application/BlackScholesHullWhiteSolver.hpp>
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
using namespace SGPP::pde;
using namespace SGPP::solver;
using namespace SGPP::base;

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace finance {

    BlackScholesHullWhiteSolver::BlackScholesHullWhiteSolver(bool useLogTransform) : ParabolicPDESolver() {
      this->bStochasticDataAlloc = false;
      this->bGridConstructed = false;
      this->myScreen = NULL;
      this->useCoarsen = false;
      this->coarsenThreshold = 0.0;
      this->adaptSolveMode = "none";
      this->refineMode = "classic";
      this->numCoarsenPoints = -1;
      this->useLogTransform = useLogTransform;
      this->refineMaxLevel = 0;
      this->useLogTransform = useLogTransform;
      this->refineMaxLevel = 0;
      this->nNeededIterations = 0;
      this->dNeededTime = 0.0;
      this->staInnerGridSize = 0;
      this->finInnerGridSize = 0;
      this->avgInnerGridSize = 0;
      this->dim_BS = 0;
      this->dim_HW = 1;

      // @todo set to random value to remove compiler warnings due to uninitialized members
      this->a = 0.0;
      this->coarsenPercent = 0.0;
      this->mus = NULL;
      this->numExecCoarsen = 0;
      this->r = 0.0;
      this->refineThreshold = 0.0;
      this->rhos = NULL;
      this->sigma = 0.0;
      this->sigmas = NULL;
      this->theta = 0.0;
    }

    BlackScholesHullWhiteSolver::~BlackScholesHullWhiteSolver() {
      if (this->bStochasticDataAlloc) {
        delete this->mus;
        delete this->sigmas;
        delete this->rhos;
      }

      if (this->myScreen != NULL) {
        delete this->myScreen;
      }
    }

    void BlackScholesHullWhiteSolver::constructGrid(BoundingBox& BoundingBox, int level) {
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


    void BlackScholesHullWhiteSolver::setStochasticData(DataVector& mus, DataVector& sigmas, DataMatrix& rhos, double r, double theta, double sigma, double a) {
      this->mus = new DataVector(mus);
      this->sigmas = new DataVector(sigmas);
      this->rhos = new DataMatrix(rhos);
      this->r = r;
      this->theta = theta;
      this->sigma = sigma;
      this->a = a;

      bStochasticDataAlloc = true;
    }


    void BlackScholesHullWhiteSolver::setProcessDimensions(int dim_BS, int dim_HW) {
      this->dim_BS = dim_BS;
      this->dim_HW = dim_HW;
    }

    void BlackScholesHullWhiteSolver::solveExplicitEuler(size_t numTimesteps, double timestepsize, size_t maxCGIterations, double epsilonCG, DataVector& alpha, bool verbose, bool generateAnimation, size_t numEvalsAnimation) {
      throw new application_exception("BlackScholesHullWhiteSolver::solveExplicitEuler : explicit Euler is not supported for BlackScholesHullWhiteSolver!");
    }

    void BlackScholesHullWhiteSolver::solveImplicitEuler(size_t numTimesteps, double timestepsize, size_t maxCGIterations, double epsilonCG, DataVector& alpha, bool verbose, bool generateAnimation, size_t numEvalsAnimation) {
      if (this->bGridConstructed && this->bStochasticDataAlloc) {
        Euler* myEuler = new Euler("ImEul", numTimesteps, timestepsize, generateAnimation, numEvalsAnimation, myScreen);
        BiCGStab* myCG = new BiCGStab(maxCGIterations, epsilonCG);

        SGppStopwatch* myStopwatch = new SGppStopwatch();
        this->staInnerGridSize = getNumberInnerGridPoints();

        std::cout << "Using Implicit Euler to solve " << numTimesteps << " timesteps:" << std::endl;
        myStopwatch->start();

        //DimensionBoundary* myBoundaries = new DimensionBoundary[2];
        BoundingBox* t = this->myGrid->getBoundingBox();

        DimensionBoundary* myBoundaries = new DimensionBoundary[dim];

        myBoundaries[0].leftBoundary = t->getIntervalOffset(0);
        myBoundaries[0].rightBoundary = t->getIntervalOffset(0) + t->getIntervalWidth(0);
        myBoundaries[1].leftBoundary = t->getIntervalOffset(1);
        myBoundaries[1].rightBoundary = t->getIntervalOffset(1) + t->getIntervalWidth(1);

        // step 1: do BS along dimension dim_BS
        myBoundaries[this->dim_BS].bDirichletLeft = true;
        myBoundaries[this->dim_BS].bDirichletRight = true;
        myBoundaries[this->dim_HW].bDirichletLeft = false;
        myBoundaries[this->dim_HW].bDirichletRight = false;

        t->setBoundary(0, myBoundaries[0]);
        t->setBoundary(1, myBoundaries[1]);
        //this->myGrid->setBoundingBox(*t);*/

        std::vector<size_t> newAlgoDimsBS(1);
        newAlgoDimsBS[0] = this->dim_BS;
        setAlgorithmicDimensions(newAlgoDimsBS);

        ModifiedBlackScholesParabolicPDESolverSystem* myBSSystem = new ModifiedBlackScholesParabolicPDESolverSystem(*this->myGrid, alpha, *this->mus, *this->sigmas, *this->rhos, this->r, timestepsize, "ImEul", this->useLogTransform, this->useCoarsen, this->coarsenThreshold, this->adaptSolveMode, this->numCoarsenPoints, this->refineThreshold, this->refineMode, this->refineMaxLevel, this->dim_HW);
        //std::cout << alpha.toString() << std::endl;
        myEuler->solve(*myCG, *myBSSystem, true, verbose);

        // step 2: do HW along dim_HW
        myBoundaries[this->dim_BS].bDirichletLeft = false;
        myBoundaries[this->dim_BS].bDirichletRight = false;
        myBoundaries[this->dim_HW].bDirichletLeft = true;
        myBoundaries[this->dim_HW].bDirichletRight = true;

        t->setBoundary(0, myBoundaries[0]);
        t->setBoundary(1, myBoundaries[1]);

        //this->myGrid->setBoundingBox(*t);*/

        std::vector<size_t> newAlgoDimsHW(1);
        newAlgoDimsHW[0] = this->dim_HW;
        setAlgorithmicDimensions(newAlgoDimsHW);

        HullWhiteParabolicPDESolverSystem* myHWSystem = new HullWhiteParabolicPDESolverSystem(*this->myGrid, alpha, this->sigma, this->theta, this->a, timestepsize, "ImEul", this->useCoarsen, this->coarsenThreshold, this->adaptSolveMode, this->numCoarsenPoints, this->refineThreshold, this->refineMode, this->refineMaxLevel, this->dim_HW);

        myEuler->solve(*myCG, *myHWSystem, true, verbose);
        this->dNeededTime = myStopwatch->stop();

        std::cout << std::endl << "Final Grid size: " << getNumberGridPoints() << std::endl;
        std::cout << "Final Grid size (inner): " << getNumberInnerGridPoints() << std::endl << std::endl << std::endl;

        std::cout << "Average Grid size: " << static_cast<double>(myBSSystem->getSumGridPointsComplete()) / static_cast<double>(numTimesteps) << std::endl;
        std::cout << "Average Grid size (Inner): " << static_cast<double>(myBSSystem->getSumGridPointsInner()) / static_cast<double>(numTimesteps) << std::endl << std::endl << std::endl;

        if (this->myScreen != NULL) {
          std::cout << "Time to solve: " << this->dNeededTime << " seconds" << std::endl;
          this->myScreen->writeEmptyLines(2);
        }

        this->finInnerGridSize = getNumberInnerGridPoints();
        this->avgInnerGridSize = static_cast<size_t>((static_cast<double>(myBSSystem->getSumGridPointsInner()) / static_cast<double>(numTimesteps)) + 0.5);
        this->nNeededIterations = myEuler->getNumberIterations();
        delete myBSSystem;
        delete myHWSystem;
        delete myCG;
        delete myEuler;
        delete myStopwatch;
        delete myBoundaries;
      } else {
        throw new application_exception("BlackScholesHullWhiteSolver::solveImplicitEuler : A grid wasn't constructed before or stochastic parameters weren't set!");
      }

    }

    void BlackScholesHullWhiteSolver::solveCrankNicolson(size_t numTimesteps, double timestepsize, size_t maxCGIterations, double epsilonCG, DataVector& alpha, size_t NumImEul) {
      throw new application_exception("BlackScholesHullWhiteSolver::solveCrankNicolson : Crank-Nicloson is not supported for BlackScholesHullWhiteSolver!");
    }

    std::vector<size_t> BlackScholesHullWhiteSolver::getAlgorithmicDimensions() {
      return this->myGrid->getAlgorithmicDimensions();
    }

    void BlackScholesHullWhiteSolver::setAlgorithmicDimensions(std::vector<size_t> newAlgoDims) {
      this->myGrid->setAlgorithmicDimensions(newAlgoDims);
    }

    void BlackScholesHullWhiteSolver::initScreen() {
      this->myScreen = new ScreenOutput();
      this->myScreen->writeTitle("SGpp - combine black scholes and Hull White Solver, 1.3.0" ,  "TUM (C) 2009-2010, by Chao qi");
      this->myScreen->writeStartSolve("combine Black Scholes and Hull White Solver");
    }

    void BlackScholesHullWhiteSolver::setEnableCoarseningData(std::string adaptSolveMode, std::string refineMode, SGPP::base::GridIndex::level_type refineMaxLevel, int numCoarsenPoints, double coarsenThreshold, double refineThreshold) {
      this->useCoarsen = true;
      this->coarsenThreshold = coarsenThreshold;
      this->refineThreshold = refineThreshold;
      this->refineMaxLevel = refineMaxLevel;
      this->adaptSolveMode = adaptSolveMode;
      this->refineMode = refineMode;
      this->numCoarsenPoints = numCoarsenPoints;
    }

    size_t BlackScholesHullWhiteSolver::getGridPointsAtMoney(std::string payoffType, double strike, double eps) {
      size_t nPoints = 0;

      if (this->useLogTransform == false) {
        if (this->bGridConstructed) {
          for (size_t i = 0; i < this->myGrid->getStorage()->size(); i++) {
            bool isAtMoney = true;
            DataVector coords(this->dim);
            this->myGridStorage->get(i)->getCoordsBB(coords, *this->myBoundingBox);

            if (payoffType == "std_euro_call" || payoffType == "std_euro_put" || payoffType == "GMIB") {
              for (size_t d = 0; d < this->dim; d++) {
                if ( ((coords.sum() / static_cast<double>(this->dim)) < (strike - eps)) || ((coords.sum() / static_cast<double>(this->dim)) > (strike + eps)) ) {
                  isAtMoney = false;
                }

              }
            } else {
              throw new application_exception("BlackScholesHullWhiteSolver::getGridPointsAtMoney : An unknown payoff-type was specified!");
            }

            if (isAtMoney == true) {
              nPoints++;
            }
          }
        } else {
          throw new application_exception("BlackScholesHullWhiteSolver::getGridPointsAtMoney : A grid wasn't constructed before!");
        }
      }

      return nPoints;
    }


    void BlackScholesHullWhiteSolver::initGridWithPayoffBSHW(DataVector& alpha, double strike, std::string payoffType, double a, double sigma) {
      double tmp;
      int timeT = 12;
      int endtime = 32;
      double c = 0.06;



      if (this->bGridConstructed) {
        double* dblFuncValues = new double[dim];

        for (size_t i = 0; i < this->myGrid->getStorage()->size(); i++) {
          std::string coords = this->myGridStorage->get(i)->getCoordsStringBB(*this->myBoundingBox);
          std::stringstream coordsStream(coords);

          for (size_t j = 0; j < this->dim; j++) {
            coordsStream >> tmp;

            dblFuncValues[j] = tmp;
          }

          if (payoffType == "std_euro_call") {
            alpha[i] = std::max<double>(dblFuncValues[this->dim_BS] - strike, 0.0);
          } else if (payoffType == "std_euro_put") {
            alpha[i] = std::max<double>(strike - dblFuncValues[this->dim_BS], 0.0);
          } else if (payoffType == "GMIB") {
            double PB = 0;
            double PT = 0;

            for (int k = (timeT + 1); k <= endtime; k++) {
              PT = exp(0.04 * (timeT - k) + 0.04 * (1 - exp(-a * (k - timeT))) / a - pow(sigma, 2.0) * pow((exp(-a * k) - exp(-a * timeT)), 2.0) * (exp(2 * a * timeT) - 1) / (4 * pow(a, 3.0)) - (1 - exp(-a * (k - timeT))) * dblFuncValues[this->dim_HW] / a);
              PB = PB + c * PT;
            }

            //std::cout << "r=" << dblFuncValues[this->dim_HW] << " PB=" <<PB <<std::endl;
            alpha[i] = std::max<double>(PB - dblFuncValues[this->dim_BS], 0.0);
          } else {
            throw new application_exception("BlackScholesSolver::initGridWithPayoffBSHW : An unknown payoff-type was specified!");
          }



        }

        delete[] dblFuncValues;

        OperationHierarchisation* myHierarchisation = SGPP::op_factory::createOperationHierarchisation(*this->myGrid);
        myHierarchisation->doHierarchisation(alpha);
        delete myHierarchisation;
      } else {
        throw new application_exception("BlackScholesSolver::initGridWithPayoffBSHW : A grid wasn't constructed before!");
      }
    }

    size_t BlackScholesHullWhiteSolver::getNeededIterationsToSolve() {
      return this->nNeededIterations;
    }

    double BlackScholesHullWhiteSolver::getNeededTimeToSolve() {
      return this->dNeededTime;
    }

    size_t BlackScholesHullWhiteSolver::getStartInnerGridSize() {
      return this->staInnerGridSize;
    }

    size_t BlackScholesHullWhiteSolver::getFinalInnerGridSize() {
      return this->finInnerGridSize;
    }

    size_t BlackScholesHullWhiteSolver::getAverageInnerGridSize() {
      return this->avgInnerGridSize;
    }

  }
}