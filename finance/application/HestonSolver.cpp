/* ****************************************************************************
* Copyright (C) 2009 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Sam Maurus (MA thesis)

#include "finance/algorithm/HestonParabolicPDESolverSystemEuroAmer.hpp"
#include "finance/application/HestonSolver.hpp"
#include "finance/application/BlackScholesSolver.hpp"
#include "solver/ode/Euler.hpp"
#include "solver/ode/CrankNicolson.hpp"
#include "solver/ode/AdamsBashforth.hpp"
#include "solver/ode/VarTimestep.hpp"
#include "solver/ode/StepsizeControlH.hpp"
#include "solver/ode/StepsizeControlBDF.hpp"
#include "solver/ode/StepsizeControlEJ.hpp"
#include "solver/sle/BiCGStab.hpp"
#include "base/grid/Grid.hpp"
#include "base/exception/application_exception.hpp"

#include "solver/sle/ConjugateGradients.hpp"

#include "base/operation/BaseOpFactory.hpp"
#include "pde/operation/PdeOpFactory.hpp"

#include <cstdlib>
#include <sstream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <complex>
#include <limits>

using namespace sg::pde;
using namespace sg::solver;
using namespace sg::base;
using namespace std;

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace sg {
  namespace finance {


    HestonSolver::HestonSolver(bool useLogTransform) : ParabolicPDESolver() {
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
      this->nNeededIterations = 0;
      this->dNeededTime = 0.0;
      this->staInnerGridSize = 0;
      this->finInnerGridSize = 0;
      this->avgInnerGridSize = 0;
      this->current_time = 0.0;
      this->tBoundaryType = "freeBoundaries";

      // @todo set to random value to remove compiler warnings due to uninitialized members
      this->dStrike = 0.0;
      this->hMatrix = NULL;
      this->kappas = NULL;
      this->numAssets = 0;
      this->r = 0.0;
      this->refineThreshold = 0.0;
      this->thetas = NULL;
      this->volvols = NULL;
    }

    HestonSolver::~HestonSolver() {
      if (this->bStochasticDataAlloc) {
        delete this->thetas;
        delete this->kappas;
        delete this->volvols;
        delete this->hMatrix;
      }

      if (this->myScreen != NULL) {
        delete this->myScreen;
      }
    }

    void HestonSolver::constructGrid(BoundingBox& BoundingBox, int level) {
      this->dim = BoundingBox.getDimensions();

      if ((dim % 2) != 0)
        throw new application_exception("HestonSolver::constructGrid : The number of dimensions in the grid is not an even number! This doesn't correspond to an integer number of assets. The number of dimensions in the grid must be divisible by two.");

      this->numAssets = this->dim / 2;
      this->levels = level;

      this->myGrid = new LinearTrapezoidBoundaryGrid(BoundingBox);

      GridGenerator* myGenerator = this->myGrid->createGridGenerator();
      myGenerator->regular(this->levels);
      delete myGenerator;

      this->myBoundingBox = this->myGrid->getBoundingBox();
      this->myGridStorage = this->myGrid->getStorage();

      this->bGridConstructed = true;
    }

    void HestonSolver::refineInitialGridWithPayoff(DataVector& alpha, double strike, std::string payoffType, double dStrikeDistance) {
      size_t nRefinements = 0;

      this->dStrike = strike;
      this->payoffType = payoffType;

      if (this->useLogTransform == false) {
        if (this->bGridConstructed) {

          DataVector refineVector(alpha.getSize());

          if (payoffType == "std_euro_call" || payoffType == "std_euro_put" || payoffType == "std_amer_put") {
            this->tBoundaryType = "Dirichlet";
            double tmp;
            double* dblFuncValues = new double[dim];
            double dDistance = 0.0;

            for (size_t i = 0; i < this->myGrid->getStorage()->size(); i++) {
              std::string coords = this->myGridStorage->get(i)->getCoordsStringBB(*(this->myBoundingBox));
              std::stringstream coordsStream(coords);

              for (size_t j = 0; j < this->dim; j++) {
                coordsStream >> tmp;

                dblFuncValues[j] = tmp;
              }

              tmp = 0.0;

              for (size_t j = 0; j < this->dim; j++) {
                tmp += dblFuncValues[j];
              }

              if (payoffType == "std_euro_call") {
                dDistance = fabs(((tmp / static_cast<double>(this->dim)) - strike));
              }

              if (payoffType == "std_euro_put" || payoffType == "std_amer_put") {
                dDistance = fabs((strike - (tmp / static_cast<double>(this->dim))));
              }

              if (dDistance <= dStrikeDistance) {
                refineVector[i] = dDistance;
                nRefinements++;
              } else {
                refineVector[i] = 0.0;
              }
            }

            delete[] dblFuncValues;

            SurplusRefinementFunctor* myRefineFunc = new SurplusRefinementFunctor(&refineVector, nRefinements, 0.0);

            this->myGrid->createGridGenerator()->refine(myRefineFunc);

            delete myRefineFunc;

            alpha.resize(this->myGridStorage->size());

            // reinit the grid with the payoff function
            initGridWithPayoff(alpha, strike, payoffType);
          } else {
            throw new application_exception("HestonSolver::refineInitialGridWithPayoff : An unsupported payoffType was specified!");
          }
        } else {
          throw new application_exception("HestonSolver::refineInitialGridWithPayoff : The grid wasn't initialized before!");
        }
      }
    }

    void HestonSolver::refineInitialGridWithPayoffToMaxLevel(DataVector& alpha, double strike, std::string payoffType, double dStrikeDistance, sg::base::GridIndex::level_type maxLevel) {
      size_t nRefinements = 0;

      this->dStrike = strike;
      this->payoffType = payoffType;

      if (this->useLogTransform == false) {
        if (this->bGridConstructed) {

          DataVector refineVector(alpha.getSize());

          if (payoffType == "std_euro_call" || payoffType == "std_euro_put" || payoffType == "std_amer_put") {
            double tmp;
            double* dblFuncValues = new double[dim];
            double dDistance = 0.0;

            this->tBoundaryType = "Dirichlet";

            for (size_t i = 0; i < this->myGrid->getStorage()->size(); i++) {
              std::string coords = this->myGridStorage->get(i)->getCoordsStringBB(*this->myBoundingBox);
              std::stringstream coordsStream(coords);

              for (size_t j = 0; j < this->dim; j++) {
                coordsStream >> tmp;

                dblFuncValues[j] = tmp;
              }

              tmp = 0.0;

              for (size_t j = 0; j < this->dim; j++) {
                tmp += dblFuncValues[j];
              }

              if (payoffType == "std_euro_call") {
                dDistance = fabs(((tmp / static_cast<double>(this->dim)) - strike));
              }

              if (payoffType == "std_euro_put" || payoffType == "std_amer_put") {
                dDistance = fabs((strike - (tmp / static_cast<double>(this->dim))));
              }

              if (dDistance <= dStrikeDistance) {
                refineVector[i] = dDistance;
                nRefinements++;
              } else {
                refineVector[i] = 0.0;
              }
            }

            delete[] dblFuncValues;

            SurplusRefinementFunctor* myRefineFunc = new SurplusRefinementFunctor(&refineVector, nRefinements, 0.0);

            this->myGrid->createGridGenerator()->refineMaxLevel(myRefineFunc, maxLevel);

            delete myRefineFunc;

            alpha.resize(this->myGridStorage->size());

            // reinit the grid with the payoff function
            initGridWithPayoff(alpha, strike, payoffType);
          } else {
            throw new application_exception("HestonSolver::refineInitialGridWithPayoffToMaxLevel : An unsupported payoffType was specified!");
          }
        } else {
          throw new application_exception("HestonSolver::refineInitialGridWithPayoffToMaxLevel : The grid wasn't initialized before!");
        }
      }
    }

    void HestonSolver::setStochasticData(DataVector& thetas_arg, DataVector& kappas_arg, DataVector& volvols_arg, DataMatrix& hMatrix_arg, double r) {
      this->thetas = new sg::base::DataVector(thetas_arg);
      this->kappas = new sg::base::DataVector(kappas_arg);
      this->volvols = new sg::base::DataVector(volvols_arg);
      this->hMatrix = new sg::base::DataMatrix(hMatrix_arg);
      this->r = r;

      bStochasticDataAlloc = true;
    }

    void HestonSolver::solveCrankNicolson(size_t numTimesteps, double timestepsize, size_t maxCGIterations, double epsilonCG, DataVector& alpha, size_t NumImEul) {
      if (this->bGridConstructed && this->bStochasticDataAlloc) {
        SLESolver* myCG = NULL;
        OperationParabolicPDESolverSystem* myHestonSystem = NULL;

        if (this->tBoundaryType == "Dirichlet") {
          myCG = new BiCGStab(maxCGIterations, epsilonCG);
          myHestonSystem = new HestonParabolicPDESolverSystemEuroAmer(*this->myGrid, alpha, *this->thetas, *this->volvols, *this->kappas, *this->hMatrix, this->r, timestepsize, "CrNic", this->dStrike, this->payoffType, this->useLogTransform, this->useCoarsen, this->coarsenThreshold, this->adaptSolveMode, this->numCoarsenPoints, this->refineThreshold, this->refineMode, this->refineMaxLevel);
        } else {
        }

        SGppStopwatch* myStopwatch = new SGppStopwatch();
        this->staInnerGridSize = getNumberInnerGridPoints();

        size_t numCNSteps;
        size_t numIESteps;

        numCNSteps = numTimesteps;

        if (numTimesteps > NumImEul) {
          numCNSteps = numTimesteps - NumImEul;
        }

        numIESteps = NumImEul;

        Euler* myEuler = new Euler("ImEul", numIESteps, timestepsize, false, 0, this->myScreen);
        CrankNicolson* myCN = new CrankNicolson(numCNSteps, timestepsize, this->myScreen);

        myStopwatch->start();

        if (numIESteps > 0) {
          std::cout << "Using Implicit Euler to solve " << numIESteps << " timesteps:" << std::endl;
          myHestonSystem->setODESolver("ImEul");
          myEuler->solve(*myCG, *myHestonSystem, false, false);
        }

        myHestonSystem->setODESolver("CrNic");
        std::cout << "Using Crank Nicolson to solve " << numCNSteps << " timesteps:" << std::endl << std::endl << std::endl << std::endl;
        myCN->solve(*myCG, *myHestonSystem, true, false);
        this->dNeededTime = myStopwatch->stop();

        std::cout << std::endl << "Final Grid size: " << getNumberGridPoints() << std::endl;
        std::cout << "Final Grid size (inner): " << getNumberInnerGridPoints() << std::endl << std::endl << std::endl;

        std::cout << "Average Grid size: " << static_cast<double>(myHestonSystem->getSumGridPointsComplete()) / static_cast<double>(numTimesteps) << std::endl;
        std::cout << "Average Grid size (Inner): " << static_cast<double>(myHestonSystem->getSumGridPointsInner()) / static_cast<double>(numTimesteps) << std::endl << std::endl << std::endl;

        if (this->myScreen != NULL) {
          std::cout << "Time to solve: " << this->dNeededTime << " seconds" << std::endl;
          this->myScreen->writeEmptyLines(2);
        }

        this->finInnerGridSize = getNumberInnerGridPoints();
        this->avgInnerGridSize = static_cast<size_t>((static_cast<double>(myHestonSystem->getSumGridPointsInner()) / static_cast<double>(numTimesteps)) + 0.5);
        this->nNeededIterations = myEuler->getNumberIterations() + myCN->getNumberIterations();

        delete myHestonSystem;
        delete myCG;
        delete myCN;
        delete myEuler;
        delete myStopwatch;

        this->current_time += (static_cast<double>(numTimesteps) * timestepsize);
      } else {
        throw new application_exception("HestonSolver::solveCrankNicolson : A grid wasn't constructed before or stochastic parameters weren't set!");
      }
    }


    void HestonSolver::initGridWithPayoff(DataVector& alpha, double strike, std::string payoffType) {
      this->dStrike = strike;
      this->payoffType = payoffType;

      if (payoffType == "std_euro_call" || payoffType == "std_euro_put" || payoffType == "std_amer_put") {
        this->tBoundaryType = "Dirichlet";
      }

      if (this->useLogTransform) {
        initLogTransformedGridWithPayoff(alpha, strike, payoffType);
      } else {
        initCartesianGridWithPayoff(alpha, strike, payoffType);
      }
    }

    double HestonSolver::get1DEuroCallPayoffValue(double assetValue, double strike) {
      if (assetValue <= strike) {
        return 0.0;
      } else {
        return assetValue - strike;
      }
    }

    void HestonSolver::solveImplicitEuler(size_t numTimesteps, double timestepsize, size_t maxCGIterations, double epsilonCG, sg::base::DataVector& alpha, bool verbose, bool generateAnimation, size_t numEvalsAnimation) {
      throw new application_exception("This scheme is not implemented for the Heston solver!");
    }

    void HestonSolver::solveExplicitEuler(size_t numTimesteps, double timestepsize, size_t maxCGIterations, double epsilonCG, sg::base::DataVector& alpha, bool verbose, bool generateAnimation, size_t numEvalsAnimation) {
      throw new application_exception("This scheme is not implemented for the Heston solver!");
    }

    std::vector<size_t> HestonSolver::getAlgorithmicDimensions() {
      return this->myGrid->getAlgorithmicDimensions();
    }

    void HestonSolver::setAlgorithmicDimensions(std::vector<size_t> newAlgoDims) {
      if (this->tBoundaryType == "freeBoundaries") {
        this->myGrid->setAlgorithmicDimensions(newAlgoDims);
      } else {
        throw new application_exception("HestonSolver::setAlgorithmicDimensions : Set algorithmic dimensions is only supported when choosing option type all!");
      }
    }

    void HestonSolver::initScreen() {
      this->myScreen = new ScreenOutput();
      this->myScreen->writeTitle("SGpp - Heston Solver, 1.0.0", "TUM (C) 2009-2010, by Sam Maurus (CSE Master's Thesis)");
      this->myScreen->writeStartSolve("Multidimensional Heston Solver");
    }

    void HestonSolver::setEnableCoarseningData(std::string adaptSolveMode, std::string refineMode, sg::base::GridIndex::level_type refineMaxLevel, int numCoarsenPoints, double coarsenThreshold, double refineThreshold) {
      this->useCoarsen = true;
      this->coarsenThreshold = coarsenThreshold;
      this->refineThreshold = refineThreshold;
      this->refineMaxLevel = refineMaxLevel;
      this->adaptSolveMode = adaptSolveMode;
      this->refineMode = refineMode;
      this->numCoarsenPoints = numCoarsenPoints;
    }

    size_t HestonSolver::getGridPointsAtMoney(std::string payoffType, double strike, double eps) {
      size_t nPoints = 0;

      if (this->useLogTransform == false) {
        if (this->bGridConstructed) {
          for (size_t i = 0; i < this->myGrid->getStorage()->size(); i++) {
            bool isAtMoney = true;
            DataVector coords(this->dim);
            this->myGridStorage->get(i)->getCoordsBB(coords, *this->myBoundingBox);

            // The stock prices come from the coordinates of the 0th, 2nd, 4th, 6th... dimensions. The other dimensions are the volatilities...they don't count for adding up stock prices.
            // So...in order to determine whether the point is at the money or not, we have to add up only the 0th, 2nd, 4th, 6th...dimension values.
            if (payoffType == "std_euro_call" || payoffType == "std_euro_put" || payoffType == "std_amer_put") {
              double stockPriceSum = 0.0;

              for (size_t j = 0; j < this->dim; j = j + 2) {
                stockPriceSum += coords[j];
              }

              if ( ((stockPriceSum / static_cast<double>(this->dim)) < (strike - eps)) || ((stockPriceSum / static_cast<double>(this->dim)) > (strike + eps)) ) {
                isAtMoney = false;
              }
            } else {
              throw new application_exception("HestonSolver::getGridPointsAtMoney : An unknown payoff-type was specified!");
            }

            if (isAtMoney == true) {
              nPoints++;
            }
          }
        } else {
          throw new application_exception("HestonSolver::getGridPointsAtMoney : A grid wasn't constructed before!");
        }
      }

      return nPoints;
    }

    void HestonSolver::initCartesianGridWithPayoff(DataVector& alpha, double strike, std::string payoffType) {
      double tmp;

      //BlackScholesSolver* bsSolver = new BlackScholesSolver();

      if (this->bGridConstructed) {
        std::ofstream fileout;

        for (size_t i = 0; i < this->myGrid->getStorage()->size(); i++) {
          std::string coords = this->myGridStorage->get(i)->getCoordsStringBB(*this->myBoundingBox);

          GridIndex* curPoint = (*myGridStorage)[i];

          std::stringstream coordsStream(coords);
          double* dblFuncValues = new double[dim];

          for (size_t j = 0; j < this->dim; j++) {
            coordsStream >> tmp;

            dblFuncValues[j] = tmp;
          }

          // print the values to a file

          if (payoffType == "std_euro_call") {
            // So now we have in alpha[i] the standard payoff value for the point S_i, v_i.
            // What we now want to do is solve the exact BS equation for every point S_i, v_i;
            //        alpha[i] = bsSolver->getAnalyticSolution1D(dblFuncValues[0], true, 0.05, pow(dblFuncValues[1], 2.0), this->r, this->dStrike);

            tmp = 0.0;

            for (size_t j = 0; j < numAssets; j++) {
              tmp += dblFuncValues[2 * j];
            }

            if (!curPoint->isInnerPoint() && dblFuncValues[1] == this->myBoundingBox->getBoundary(1).rightBoundary) {
              // Dirichlet condition when v -> inf is that U = S
              alpha[i] = dblFuncValues[0]; //(this->myBoundingBox->getBoundary(0).rightBoundary-strike)/(this->myBoundingBox->getBoundary(0).rightBoundary)*dblFuncValues[0];
            } else if (!curPoint->isInnerPoint() && dblFuncValues[0] == this->myBoundingBox->getBoundary(0).rightBoundary) {
              // Set boundary to be the linear function at s_max
              double normalPayoff = std::max<double>(((tmp / static_cast<double>(numAssets)) - strike), 0.0);
              double vRange = this->myBoundingBox->getBoundary(1).rightBoundary - this->myBoundingBox->getBoundary(1).leftBoundary;
              double sPayoffDiff = dblFuncValues[0] - normalPayoff;
              alpha[i] = normalPayoff + ((dblFuncValues[1] - this->myBoundingBox->getBoundary(1).leftBoundary) / vRange) * sPayoffDiff;
            } else {
              // Payoff function
              alpha[i] = std::max<double>(((tmp / static_cast<double>(numAssets)) - strike), 0.0);
            }
          } else if (payoffType == "std_euro_put") {
            if (!curPoint->isInnerPoint()) {
              if (numAssets == 1 && dblFuncValues[1] == this->myBoundingBox->getBoundary(1).rightBoundary) {
                // Vmax boundary for a single asset. Strike price
                alpha[i] = strike;
              } else if (numAssets == 1 && dblFuncValues[0] == this->myBoundingBox->getBoundary(0).rightBoundary) {
                // Smax boundary for a single asset. exponential function
                //double constantC = 20.0;
                double constantB = (strike * pow(exp(this->myBoundingBox->getBoundary(1).leftBoundary - this->myBoundingBox->getBoundary(1).rightBoundary), 20.0)) / (1 - exp(this->myBoundingBox->getBoundary(1).leftBoundary - this->myBoundingBox->getBoundary(1).rightBoundary));
                double constantA = strike + constantB;
                alpha[i] = constantA * pow(exp(dblFuncValues[1] - this->myBoundingBox->getBoundary(1).rightBoundary), 20.0) - constantB;
              } else {

                // Get the payoff function value for this point
                double payoffFuncVal = 0.0;

                for (size_t k = 0; k < numAssets; k++) {
                  payoffFuncVal += dblFuncValues[2 * k];
                }

                payoffFuncVal = std::max<double>(strike - payoffFuncVal / (static_cast<double>(numAssets)), 0.0);

                // Get the max-volatility function value for this point
                double maxVolFuncVal = 0.0;

                for (size_t k = 0; k < numAssets; k++) {
                  double stockPrice = dblFuncValues[2 * k];
                  double sMaxValue = this->myBoundingBox->getBoundary(2 * k).rightBoundary;
                  maxVolFuncVal +=  (strike / (2.0 * sMaxValue)) * stockPrice;
                }

                maxVolFuncVal = strike - maxVolFuncVal / (static_cast<double>(numAssets));

                // Get the fraction that we are away from the v=vMin value in the direction of vMax
                double volFraction = 0.0;

                for (size_t k = 0; k < numAssets; k++) {
                  double vValue = dblFuncValues[2 * k + 1];
                  double vRightValue = this->myBoundingBox->getBoundary(2 * k + 1).rightBoundary;
                  double vLeftValue = this->myBoundingBox->getBoundary(2 * k + 1).leftBoundary;
                  double vRange = vRightValue - vLeftValue;
                  volFraction += (vValue - vLeftValue) / vRange;
                }

                volFraction = volFraction / (static_cast<double>(numAssets));

                // Now set the value of the point as the payoff function + fraction*maxVol function
                alpha[i] = std::max<double>(payoffFuncVal + volFraction * (maxVolFuncVal - payoffFuncVal), 0.0);
              }
            } else {
              // Non-boundary point. Just set the payoff value.
              tmp = 0.0;

              for (size_t j = 0; j < numAssets; j++) {
                tmp += dblFuncValues[2 * j];
              }

              alpha[i] = std::max<double>(strike - ((tmp / static_cast<double>(numAssets))), 0.0);
            }
          } else {
            throw new application_exception("HestonSolver::initCartesianGridWithPayoff : An unknown payoff-type was specified!");
          }

          delete[] dblFuncValues;
        }

        //    fileout.close();

        // determine the number of grid points for both grids
        size_t numTotalGridPoints = myGridStorage->size();
        //size_t numInnerGridPoints = myGridStorage->getNumInnerPoints();

        size_t numNonZeroInner = 0;

        for (size_t i = 0; i < numTotalGridPoints; i++) {
          sg::base::GridIndex* curPoint = (*myGridStorage)[i];

          if (curPoint->isInnerPoint() && alpha.get(i) != 0)
            numNonZeroInner++;
        }

        int k = 0;
        k++;

        OperationHierarchisation* myHierarchisation = sg::op_factory::createOperationHierarchisation(*this->myGrid);
        myHierarchisation->doHierarchisation(alpha);
        delete myHierarchisation;
      } else {
        throw new application_exception("HestonSolver::initCartesianGridWithPayoff : A grid wasn't constructed before!");
      }
    }

    void HestonSolver::initLogTransformedGridWithPayoff(DataVector& alpha, double strike, std::string payoffType) {
      double tmp;

      //BlackScholesSolver* bsSolver = new BlackScholesSolver();

      if (this->bGridConstructed) {
        for (size_t i = 0; i < this->myGrid->getStorage()->size(); i++) {
          std::string coords = this->myGridStorage->get(i)->getCoordsStringBB(*this->myBoundingBox);
          std::stringstream coordsStream(coords);
          double* dblFuncValues = new double[dim];

          GridIndex* curPoint = (*myGridStorage)[i];

          for (size_t j = 0; j < this->dim; j++) {
            coordsStream >> tmp;

            dblFuncValues[j] = tmp;
          }

          if (payoffType == "std_euro_call") {
            // Here we have to be a little careful in comparison to the Black-Scholes solver.
            // In the Black-Scholes solver, every dimension represents an asset price, so determining the payoff value at this particular point
            // is simply a case of adding up all the dimension (asset) values and inserting that into the payoff function.
            // In the Heston model, however, only every second dimension represents an asset price. Every other dimension represents a variance, which
            // we MUST NOT use to determine the payoff value.
            // Thus, in oder to determine the payoff value, we have to only iterate over the first, third, fifth... dimensions.

            if (!curPoint->isInnerPoint()) {
              double sumStockPrices = 0.0;
              double kMultiplier = 0;

              for (size_t k = 0; k < numAssets; k++) {
                // Get the sum of all stock prices
                sumStockPrices += exp(dblFuncValues[2 * k]);

                // At the same time, figure out the fraction of the way we are away from the max-v value for this stock
                double vValue = dblFuncValues[2 * k + 1];
                double vRightValue = this->myBoundingBox->getBoundary(2 * k + 1).rightBoundary;
                double vLeftValue = this->myBoundingBox->getBoundary(2 * k + 1).leftBoundary;
                double vRange = vRightValue - vLeftValue;
                kMultiplier += (vRightValue - vValue) / vRange;
              }

              alpha[i] = std::max<double>((sumStockPrices - strike * kMultiplier) / static_cast<double>(numAssets), 0.0);
            } else {
              // Non-boundary point. Just set the payoff value.
              tmp = 0.0;

              for (size_t j = 0; j < numAssets; j++) {
                tmp += exp(dblFuncValues[2 * j]);
              }

              alpha[i] = std::max<double>(((tmp / static_cast<double>(numAssets)) - strike), 0.0);
            }
          } else if (payoffType == "std_euro_put") {
            if (!curPoint->isInnerPoint()) {
              if (numAssets == 1 && dblFuncValues[1] == this->myBoundingBox->getBoundary(1).rightBoundary) {
                // Vmax boundary for a single asset. Strike price
                alpha[i] = strike;
              } else if (numAssets == 1 && dblFuncValues[0] == this->myBoundingBox->getBoundary(0).rightBoundary) {
                // Smax boundary for a single asset. exponential function
                //double constantC = 20.0;
                double constantB = (strike * pow(exp(this->myBoundingBox->getBoundary(1).leftBoundary - this->myBoundingBox->getBoundary(1).rightBoundary), 20.0)) / (1 - exp(this->myBoundingBox->getBoundary(1).leftBoundary - this->myBoundingBox->getBoundary(1).rightBoundary));
                double constantA = strike + constantB;
                alpha[i] = constantA * pow(exp(dblFuncValues[1] - this->myBoundingBox->getBoundary(1).rightBoundary), 20.0) - constantB;
              } else {
                // Get the payoff function value for this point
                double payoffFuncVal = 0.0;

                for (size_t k = 0; k < numAssets; k++) {
                  payoffFuncVal += exp(dblFuncValues[2 * k]);
                }

                payoffFuncVal = std::max<double>(strike - payoffFuncVal / (static_cast<double>(numAssets)), 0.0);

                // Get the max-volatility function value for this point
                double maxVolFuncVal = 0.0;

                for (size_t k = 0; k < numAssets; k++) {
                  double stockPrice = exp(dblFuncValues[2 * k]);
                  double sMinValue = exp(this->myBoundingBox->getBoundary(2 * k).leftBoundary);
                  double sMaxValue = exp(this->myBoundingBox->getBoundary(2 * k).rightBoundary);
                  maxVolFuncVal +=  (strike / (2.0 * (sMaxValue - sMinValue))) * (stockPrice - sMinValue);
                }

                maxVolFuncVal = strike - maxVolFuncVal / (static_cast<double>(numAssets));

                // Get the fraction that we are away from the v=vMin value in the direction of vMax
                double volFraction = 0.0;

                for (size_t k = 0; k < numAssets; k++) {
                  double vValue = dblFuncValues[2 * k + 1];
                  double vRightValue = this->myBoundingBox->getBoundary(2 * k + 1).rightBoundary;
                  double vLeftValue = this->myBoundingBox->getBoundary(2 * k + 1).leftBoundary;
                  double vRange = vRightValue - vLeftValue;
                  volFraction += (vValue - vLeftValue) / vRange;
                }

                volFraction = volFraction / (static_cast<double>(numAssets));

                // Now set the value of the point as the payoff function + fraction*maxVol function
                alpha[i] = std::max<double>(payoffFuncVal + volFraction * (maxVolFuncVal - payoffFuncVal), 0.0);
              }
            } else {
              // Non-boundary point. Just set the payoff value.
              tmp = 0.0;

              for (size_t j = 0; j < numAssets; j++) {
                tmp += exp(dblFuncValues[2 * j]);
              }

              alpha[i] = std::max<double>(strike - ((tmp / static_cast<double>(numAssets))), 0.0);
            }
          } else {
            throw new application_exception("HestonSolver::initLogTransformedGridWithPayoff : An unknown payoff-type was specified!");
          }

          delete[] dblFuncValues;
        }

        OperationHierarchisation* myHierarchisation = sg::op_factory::createOperationHierarchisation(*this->myGrid);
        myHierarchisation->doHierarchisation(alpha);
        delete myHierarchisation;
      } else {
        throw new application_exception("HestonSolver::initLogTransformedGridWithPayoff : A grid wasn't constructed before!");
      }
    }

    double HestonSolver::evalOption(std::vector<double>& eval_point, sg::base::DataVector& alpha) {
      std::vector<double> trans_eval = eval_point;

      // apply needed coordinate transformations
      if (this->useLogTransform) {
        for (size_t i = 0; i < eval_point.size(); i = i + 2) { // Here we know that the variance dimension is not log-transformed, so we shouldn't log its value
          trans_eval[i] = log(trans_eval[i]);
        }
      }

      sg::base::OperationEval* myEval = sg::op_factory::createOperationEval(*this->myGrid);
      double result = myEval->eval(alpha, trans_eval);
      delete myEval;

      return result;
    }

    void HestonSolver::transformPoint(sg::base::DataVector& point) {
      sg::base::DataVector tmp_point(point);

      // apply needed coordinate transformations
      if (this->useLogTransform) {
        for (size_t i = 0; i < point.getSize(); i++) {
          tmp_point[i] = log(point[i]);
        }
      }

      point = tmp_point;
    }

    void HestonSolver::resetSolveTime() {
      this->current_time = 0.0;
    }

    size_t HestonSolver::getNeededIterationsToSolve() {
      return this->nNeededIterations;
    }

    double HestonSolver::getNeededTimeToSolve() {
      return this->dNeededTime;
    }

    size_t HestonSolver::getStartInnerGridSize() {
      return this->staInnerGridSize;
    }

    size_t HestonSolver::getFinalInnerGridSize() {
      return this->finInnerGridSize;
    }

    size_t HestonSolver::getAverageInnerGridSize() {
      return this->avgInnerGridSize;
    }

    double EvaluateHestonClosedFormIntegralFunction(double phi, double xi, double theta, double kappa, double rho, double r, double T, double K, double S, double v, int type) {
      double a, b, u, x;

      if (type == 1) {
        u = 0.5;
        b = kappa - rho * xi;
      } else {
        u = -0.5;
        b = kappa;
      }

      a = kappa * theta;
      x = log(S);

      //  Build d
      complex<double> dBuilder1;
      complex<double> dBuilder2;
      complex<double> d;
      dBuilder1 = std::complex<double>(0.0, 1.0);
      dBuilder1 = pow(dBuilder1 * phi * rho * xi - b, 2.0);
      dBuilder2 = std::complex<double>(0.0, 1.0);
      dBuilder2 = dBuilder2 * phi;
      dBuilder2 = dBuilder2 * 2.0 * u - pow(phi, 2.0);
      dBuilder2 = dBuilder2 * pow(xi, 2.0);
      d = sqrt(dBuilder1 - dBuilder2);

      // Build g
      complex<double> g;
      g = std::complex<double>(0.0, 1.0);
      g = (b - rho * xi * phi * g + d) / (b - rho * xi * phi * g - d);


      // Build C
      complex<double> C;
      C = std::complex<double>(0.0, 1.0);
      C = (r * phi * C * T) + (a / pow(xi, 2.0)) * ((b - rho * xi * phi * C + d) * T - 2.0 * log( (1.0 - g * exp(d * T)) / (1.0 - g)  ));

      // Build D
      complex<double> D;
      D = std::complex<double>(0.0, 1.0);
      D = ((b - rho * xi * phi * D + d) / (pow(xi, 2.0))) * ((1.0 - exp(d * T)) / (1.0 - g * exp(d * T)));

      // Build f
      complex<double> f;
      f = complex<double>(0.0, 1.0);
      f = exp(C + D * v + f * phi * x);


      // Build realArgument
      complex<double> realArgument;
      realArgument = std::complex<double>(0.0, 1.0);
      realArgument = (exp(0.0 - realArgument * phi * log(K)) * f) / (realArgument * phi);

      return real(realArgument);
    }

    void HestonSolver::EvaluateHestonExactSurface(DataVector& alpha, double maturity) {
      if (!this->bGridConstructed)
        throw new application_exception("HestonSolver::EvaluateHestonPriceExact : The grid wasn't initialized before!");

      if (this->numAssets != 1 || this->payoffType != "std_euro_call")
        throw new application_exception("HestonSolver::EvaluateHestonPriceExact : Can only solve in closed form for a European call option with one asset!");

      double tmp;

      for (size_t i = 0; i < this->myGrid->getStorage()->size(); i++) {
        std::string coords = this->myGridStorage->get(i)->getCoordsStringBB(*this->myBoundingBox);
        std::stringstream coordsStream(coords);
        double* dblFuncValues = new double[dim];

        for (size_t j = 0; j < this->dim; j++) {
          coordsStream >> tmp;

          if (this->useLogTransform && j == 0)
            dblFuncValues[j] = exp(tmp);
          else
            dblFuncValues[j] = tmp;
        }

        alpha[i] = EvaluateHestonPriceExact(dblFuncValues[0], dblFuncValues[1], this->volvols->get(0), this->thetas->get(0), this->kappas->get(0), this->hMatrix->get(0, 1), this->r, maturity, this->dStrike) ;
        delete dblFuncValues;
      }

      OperationHierarchisation* myHierarchisation = sg::op_factory::createOperationHierarchisation(*this->myGrid);
      myHierarchisation->doHierarchisation(alpha);
      delete myHierarchisation;
    }

    void HestonSolver::EvaluateHestonExactSurfacePut(DataVector& alpha, double maturity) {
      if (!this->bGridConstructed)
        throw new application_exception("HestonSolver::EvaluateHestonPriceExact : The grid wasn't initialized before!");

      if (this->numAssets != 1 || this->payoffType != "std_euro_put")
        throw new application_exception("HestonSolver::EvaluateHestonPriceExact : Can only solve in closed form for a European put option with one asset!");

      double tmp;

      for (size_t i = 0; i < this->myGrid->getStorage()->size(); i++) {
        std::string coords = this->myGridStorage->get(i)->getCoordsStringBB(*this->myBoundingBox);
        std::stringstream coordsStream(coords);
        double* dblFuncValues = new double[dim];

        for (size_t j = 0; j < this->dim; j++) {
          coordsStream >> tmp;

          if (this->useLogTransform && j == 0)
            dblFuncValues[j] = exp(tmp);
          else
            dblFuncValues[j] = tmp;
        }

        alpha[i] = EvaluateHestonPriceExactPut(dblFuncValues[0], dblFuncValues[1], this->volvols->get(0), this->thetas->get(0), this->kappas->get(0), this->hMatrix->get(0, 1), this->r, maturity, this->dStrike) ;
        delete dblFuncValues;
      }

      OperationHierarchisation* myHierarchisation = sg::op_factory::createOperationHierarchisation(*this->myGrid);
      myHierarchisation->doHierarchisation(alpha);
      delete myHierarchisation;
    }

    void HestonSolver::CompareHestonBs1d(double maturity, double v) {
      if (this->numAssets != 1 || this->payoffType != "std_euro_call")
        throw new application_exception("HestonSolver::EvaluateHestonPriceExact : Can only solve in closed form for a European call option with one asset!");

      size_t dim1d = 1;
      int levels1d = this->levels;

      // Build a new 1d grid for stock price
      DimensionBoundary* boundaries1d = new sg::base::DimensionBoundary[dim1d];

      boundaries1d[0].leftBoundary = this->myBoundingBox->getBoundary(0).leftBoundary;
      boundaries1d[0].rightBoundary = this->myBoundingBox->getBoundary(0).rightBoundary;
      boundaries1d[0].bDirichletLeft = this->myBoundingBox->getBoundary(0).bDirichletLeft;
      boundaries1d[0].bDirichletRight = this->myBoundingBox->getBoundary(0).bDirichletRight;

      BoundingBox* boundingBox1d = new BoundingBox(dim1d, boundaries1d);

      Grid* grid1d = new LinearTrapezoidBoundaryGrid(*boundingBox1d);

      GridGenerator* myGenerator = grid1d->createGridGenerator();
      myGenerator->regular(levels1d);

      sg::base::DataVector* alphaHeston = new sg::base::DataVector(grid1d->getSize());
      sg::base::DataVector* alphaBS = new sg::base::DataVector(grid1d->getSize());

      EvaluateHestonExact1d(*alphaHeston, grid1d, boundingBox1d, maturity, v);
      EvaluateBsExact1d(*alphaBS, grid1d, boundingBox1d, maturity, 0.259);

      // Switch over the grids so we can print the 1d ones
      Grid* gridTemp = this->myGrid;
      this->myGrid = grid1d;

      printGrid(*alphaHeston, 50, "hestonExact1d.gnuplot");
      printGrid(*alphaBS, 50, "bsExact1d.gnuplot");

      alphaHeston->sub(*alphaBS);

      printGrid(*alphaHeston, 50, "hestonMinusBsExact1d.gnuplot");

      this->myGrid = gridTemp;

      delete boundaries1d;
      delete grid1d;
      delete myGenerator;
      delete alphaHeston;
      delete alphaBS;
    }

    void HestonSolver::EvaluateHestonExact1d(DataVector& alpha, Grid* grid1d, BoundingBox* boundingBox1d, double maturity, double v) {
      double tmp;

      for (size_t i = 0; i < grid1d->getStorage()->size(); i++) {
        std::string coords = grid1d->getStorage()->get(i)->getCoordsStringBB(*boundingBox1d);
        std::stringstream coordsStream(coords);
        double* dblFuncValues = new double[1];
        coordsStream >> tmp;
        dblFuncValues[0] = tmp;

        alpha[i] = EvaluateHestonPriceExact(dblFuncValues[0], v, this->volvols->get(0), this->thetas->get(0), this->kappas->get(0), this->hMatrix->get(0, 1), this->r, maturity, this->dStrike) ;
        delete dblFuncValues;
      }

      OperationHierarchisation* myHierarchisation = sg::op_factory::createOperationHierarchisation(*grid1d);
      myHierarchisation->doHierarchisation(alpha);
      delete myHierarchisation;
    }

    void HestonSolver::EvaluateBsExact1d(DataVector& alpha, Grid* grid1d, BoundingBox* boundingBox1d, double maturity, double sigma) {
      double tmp;

      sg::finance::BlackScholesSolver* myBSSolver = new sg::finance::BlackScholesSolver(false);

      for (size_t i = 0; i < grid1d->getStorage()->size(); i++) {
        std::string coords = grid1d->getStorage()->get(i)->getCoordsStringBB(*boundingBox1d);
        std::stringstream coordsStream(coords);
        double* dblFuncValues = new double[1];
        coordsStream >> tmp;
        dblFuncValues[0] = tmp;

        alpha[i] = myBSSolver->getAnalyticSolution1D(dblFuncValues[0], true, maturity, sigma, this->r, this->dStrike);
        delete dblFuncValues;
      }

      OperationHierarchisation* myHierarchisation = sg::op_factory::createOperationHierarchisation(*grid1d);
      myHierarchisation->doHierarchisation(alpha);
      delete myHierarchisation;
      delete myBSSolver;
    }

    double HestonSolver::EvaluateHestonPriceExact(double S, double v, double xi, double theta, double kappa, double rho, double r, double T, double K) {
      double int1 = 0.5 + (1.0 / M_PI) * GaussLobattoInt(0.001, 1000.0, 1e-10, 100000, xi, theta, kappa, rho, r, T, K, S, v, 1);
      double int2 = 0.5 + (1.0 / M_PI) * GaussLobattoInt(0.001, 1000.0, 1e-10, 100000, xi, theta, kappa, rho, r, T, K, S, v, 2);
      return S * int1 - K * exp((-1.0) * r * T) * int2;
    }

    double HestonSolver::EvaluateHestonPriceExact(double S, double v, double maturity) {
      return EvaluateHestonPriceExact(S, v, this->volvols->get(0), this->thetas->get(0), this->kappas->get(0), this->hMatrix->get(0, 1), this->r, maturity, this->dStrike);
    }

    double HestonSolver::EvaluateHestonPriceExactPut(double S, double v, double xi, double theta, double kappa, double rho, double r, double T, double K) {
      double int1 = 0.5 + (1.0 / M_PI) * GaussLobattoInt(0.001, 1000.0, 1e-10, 100000, xi, theta, kappa, rho, r, T, K, S, v, 1);
      double int2 = 0.5 + (1.0 / M_PI) * GaussLobattoInt(0.001, 1000.0, 1e-10, 100000, xi, theta, kappa, rho, r, T, K, S, v, 2);
      return S * int1 - K * exp((-1.0) * r * T) * int2 + K * exp((-1.0) * r * T) - S;
    }

    double HestonSolver::EvaluateHestonPriceExactPut(double S, double v, double maturity) {
      return EvaluateHestonPriceExactPut(S, v, this->volvols->get(0), this->thetas->get(0), this->kappas->get(0), this->hMatrix->get(0, 1), this->r, maturity, this->dStrike);
    }

    void HestonSolver::GetBsExactSolution(sg::base::DataVector& alphaBS, double maturity) {
      sg::finance::BlackScholesSolver* myBSSolver = new sg::finance::BlackScholesSolver(false);

      double S, v;

      for (size_t i = 0; i < this->myGridStorage->size(); i++) {
        std::string coords = this->myGridStorage->get(i)->getCoordsStringBB(*this->myBoundingBox);
        std::stringstream coordsStream(coords);
        coordsStream >> S;
        coordsStream >> v;

        if (this->useLogTransform)
          alphaBS[i] = myBSSolver->getAnalyticSolution1D(exp(S), true, maturity, sqrt(v), this->r, this->dStrike);
        else
          alphaBS[i] = myBSSolver->getAnalyticSolution1D(S, true, maturity, sqrt(v), this->r, this->dStrike);
      }

      OperationHierarchisation* myHierarchisation = sg::op_factory::createOperationHierarchisation(*this->myGrid);
      myHierarchisation->doHierarchisation(alphaBS);
      delete myHierarchisation;
    }

    void HestonSolver::CompareHestonBsExact(sg::base::DataVector& alpha, double maturity) {
      DataVector* alphaHeston = new DataVector(getNumberGridPoints());
      DataVector* alphaBS = new DataVector(getNumberGridPoints());

      EvaluateHestonExactSurface(*alphaHeston, maturity);
      GetBsExactSolution(*alphaBS, maturity);

      // Find the difference (heston - BS)
      for (size_t i = 0; i < this->myGridStorage->size(); i++) {
        alpha[i] = (*alphaHeston)[i] - (*alphaBS)[i];
      }
    }

    void HestonSolver::CompareHestonNumericToBsExact(sg::base::DataVector& alphaHestonNumeric, sg::base::DataVector& alphaBS, sg::base::DataVector& error, double maturity) {
      GetBsExactSolution(alphaBS, maturity);

      // Find the difference (heston - BS)
      for (size_t i = 0; i < this->myGridStorage->size(); i++) {
        error[i] = alphaHestonNumeric[i] - alphaBS[i];
      }
    }


    double HestonSolver::GaussLobattoIntStep(
      double a, double b,
      double fa, double fb,
      size_t& neval,
      size_t maxeval,
      double acc
      , double xi, double theta, double kappa, double rho, double r, double T, double K, double S, double v, int type) {

      // Constants used in the algorithm
      const double alpha = std::sqrt(2.0 / 3.0);
      const double beta  = 1.0 / std::sqrt(5.0);

      if (neval >= maxeval) {
        throw new application_exception("HestonSolver::Gauss-Lobatto : Maximum number of evaluations reached in GaussLobatto.");
      }

      // Here the abcissa points and function values for both the 4-point
      // and the 7-point rule are calculated (the points at the end of
      // interval come from the function call, i.e., fa and fb. Also note
      // the 7-point rule re-uses all the points of the 4-point rule.)
      const double h = (b - a) / 2;
      const double m = (a + b) / 2;

      const double mll = m - alpha * h;
      const double ml = m - beta * h;
      const double mr = m + beta * h;
      const double mrr = m + alpha * h;

      const double fmll = EvaluateHestonClosedFormIntegralFunction(mll, xi, theta, kappa, rho, r, T, K, S, v, type);
      const double fml = EvaluateHestonClosedFormIntegralFunction(ml, xi, theta, kappa, rho, r, T, K, S, v, type);
      const double fm  = EvaluateHestonClosedFormIntegralFunction(m, xi, theta, kappa, rho, r, T, K, S, v, type);
      const double fmr = EvaluateHestonClosedFormIntegralFunction(mr, xi, theta, kappa, rho, r, T, K, S, v, type);
      const double fmrr = EvaluateHestonClosedFormIntegralFunction(mrr, xi, theta, kappa, rho, r, T, K, S, v, type);
      neval += 5;

      // Both the 4-point and 7-point rule integrals are evaluted
      const double integral2 = (h / 6) * (fa + fb + 5 * (fml + fmr));
      const double integral1 = (h / 1470) * (77 * (fa + fb)
                                             + 432 * (fmll + fmrr) + 625 * (fml + fmr) + 672 * fm);

      // The difference betwen the 4-point and 7-point integrals is the
      // estimate of the accuracy
      const double estacc = (integral1 - integral2);

      // The volatile keyword should prevent the floating point
      // destination value from being stored in extended precision
      // registers which actually have a very different
      // std::numeric_limits<double>::epsilon().
      volatile double dist = acc + estacc;

      if (dist == acc || mll <= a || b <= mrr) {
        if (!(m > a && b > m)) {
          throw new application_exception("HestonSolver::Gauss-Lobatto : Integration reached an interval with no more machine numbers!");
        }

        return integral1;
      } else {
        return  GaussLobattoIntStep(a, mll, fa, fmll, neval, maxeval, acc, xi, theta, kappa, rho, r, T, K, S, v, type)
                + GaussLobattoIntStep(mll, ml, fmll, fml, neval, maxeval, acc, xi, theta, kappa, rho, r, T, K, S, v, type)
                + GaussLobattoIntStep(ml, m, fml, fm, neval, maxeval, acc, xi, theta, kappa, rho, r, T, K, S, v, type)
                + GaussLobattoIntStep(m, mr, fm, fmr, neval, maxeval, acc, xi, theta, kappa, rho, r, T, K, S, v, type)
                + GaussLobattoIntStep(mr, mrr, fmr, fmrr, neval, maxeval, acc, xi, theta, kappa, rho, r, T, K, S, v, type)
                + GaussLobattoIntStep(mrr, b, fmrr, fb, neval, maxeval, acc, xi, theta, kappa, rho, r, T, K, S, v, type);

      }
    }

    double HestonSolver::GaussLobattoInt(double a, double b,
                                         double abstol,
                                         size_t maxeval
                                         , double xi, double theta, double kappa, double rho, double r, double T, double K, double S, double v, int type) {
      const double tol_epsunit = abstol / std::numeric_limits<double>::epsilon();
      size_t neval = 0;
      return GaussLobattoIntStep(a, b,
                                 EvaluateHestonClosedFormIntegralFunction(a, xi, theta, kappa, rho, r, T, K, S, v, type), EvaluateHestonClosedFormIntegralFunction(b, xi, theta, kappa, rho, r, T, K, S, v, type),
                                 neval,
                                 maxeval,
                                 tol_epsunit, xi, theta, kappa, rho, r, T, K, S, v, type);
    }

    void HestonSolver::CompareHestonSolutionToExact(sg::base::DataVector* solution, sg::base::DataVector* exact, std::string filename, size_t PointsPerDimension) {
      DimensionBoundary dimOne;
      DimensionBoundary dimTwo;
      std::ofstream fileout;

      fileout.open(filename.c_str());
      OperationEval* myEval = sg::op_factory::createOperationEval(*myGrid);


      if (myGrid->getStorage()->dim() == 2) {
        dimOne = myGrid->getBoundingBox()->getBoundary(0);
        dimTwo = myGrid->getBoundingBox()->getBoundary(1);

        double offset_x = dimOne.leftBoundary;
        double offset_y = dimTwo.leftBoundary;
        double inc_x = ((dimOne.rightBoundary - dimOne.leftBoundary) / double(PointsPerDimension - 1));
        double inc_y = ((dimTwo.rightBoundary - dimTwo.leftBoundary) / double(PointsPerDimension - 1));

        size_t points = (size_t)PointsPerDimension;

        for (size_t i = 0; i < points; i++) {
          for (size_t j = 0; j < points; j++) {
            std::vector<double> point;
            point.push_back(offset_x + (((double)(i))*inc_x));
            point.push_back(offset_y + (((double)(j))*inc_y));
            fileout << (offset_x + ((double)(i))*inc_x) << " " << (offset_y + ((double)(j))*inc_y) << " " << (myEval->eval(*solution, point) - myEval->eval(*exact, point)) << std::endl;
          }

          fileout << std::endl;
        }
      }

      delete myEval;
      // close filehandle
      fileout.close();
    }

    double HestonSolver::EvalSinglePoint1Asset(double s, double v, DataVector& alphaVec) {
      OperationEval* myEval = sg::op_factory::createOperationEval(*myGrid);
      std::vector<double> point;
      point.push_back(s);
      point.push_back(v);
      return myEval->eval(alphaVec, point);
    }

  }
}

