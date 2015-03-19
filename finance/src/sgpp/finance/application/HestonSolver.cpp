// Copyright (C) 2008-today The SG++ Project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/finance/algorithm/HestonParabolicPDESolverSystemEuroAmer.hpp>
#include <sgpp/finance/application/HestonSolver.hpp>
#include <sgpp/finance/application/BlackScholesSolver.hpp>
#include <sgpp/solver/ode/Euler.hpp>
#include <sgpp/solver/ode/CrankNicolson.hpp>
#include <sgpp/solver/ode/AdamsBashforth.hpp>
#include <sgpp/solver/ode/VarTimestep.hpp>
#include <sgpp/solver/ode/StepsizeControlH.hpp>
#include <sgpp/solver/ode/StepsizeControlBDF.hpp>
#include <sgpp/solver/ode/StepsizeControlEJ.hpp>
#include <sgpp/solver/sle/BiCGStab.hpp>
#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/exception/application_exception.hpp>

#include <sgpp/solver/sle/ConjugateGradients.hpp>

#include <sgpp/base/operation/BaseOpFactory.hpp>
#include <sgpp/pde/operation/PdeOpFactory.hpp>

#include <cstdlib>
#include <sstream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <complex>
#include <limits>

#include <sgpp/globaldef.hpp>


using namespace SGPP::pde;
using namespace SGPP::solver;
using namespace SGPP::base;
using namespace std;

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif


namespace SGPP {
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

      this->myGrid = new LinearTruncatedBoundaryGrid(BoundingBox);

      GridGenerator* myGenerator = this->myGrid->createGridGenerator();
      myGenerator->regular(this->levels);
      delete myGenerator;

      this->myBoundingBox = this->myGrid->getBoundingBox();
      this->myGridStorage = this->myGrid->getStorage();

      this->bGridConstructed = true;
    }

    void HestonSolver::refineInitialGridWithPayoff(DataVector& alpha, float_t strike, std::string payoffType, float_t dStrikeDistance) {
      size_t nRefinements = 0;

      this->dStrike = strike;
      this->payoffType = payoffType;

      if (this->useLogTransform == false) {
        if (this->bGridConstructed) {

          DataVector refineVector(alpha.getSize());

          if (payoffType == "std_euro_call" || payoffType == "std_euro_put" || payoffType == "std_amer_put") {
            this->tBoundaryType = "Dirichlet";
            float_t tmp;
            float_t* dblFuncValues = new float_t[dim];
            float_t dDistance = 0.0;

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
                dDistance = fabs(((tmp / static_cast<float_t>(this->dim)) - strike));
              }

              if (payoffType == "std_euro_put" || payoffType == "std_amer_put") {
                dDistance = fabs((strike - (tmp / static_cast<float_t>(this->dim))));
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

    void HestonSolver::refineInitialGridWithPayoffToMaxLevel(DataVector& alpha, float_t strike, std::string payoffType, float_t dStrikeDistance, SGPP::base::GridIndex::level_type maxLevel) {
      size_t nRefinements = 0;

      this->dStrike = strike;
      this->payoffType = payoffType;

      if (this->useLogTransform == false) {
        if (this->bGridConstructed) {

          DataVector refineVector(alpha.getSize());

          if (payoffType == "std_euro_call" || payoffType == "std_euro_put" || payoffType == "std_amer_put") {
            float_t tmp;
            float_t* dblFuncValues = new float_t[dim];
            float_t dDistance = 0.0;

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
                dDistance = fabs(((tmp / static_cast<float_t>(this->dim)) - strike));
              }

              if (payoffType == "std_euro_put" || payoffType == "std_amer_put") {
                dDistance = fabs((strike - (tmp / static_cast<float_t>(this->dim))));
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

    void HestonSolver::setStochasticData(DataVector& thetas_arg, DataVector& kappas_arg, DataVector& volvols_arg, DataMatrix& hMatrix_arg, float_t r) {
      this->thetas = new SGPP::base::DataVector(thetas_arg);
      this->kappas = new SGPP::base::DataVector(kappas_arg);
      this->volvols = new SGPP::base::DataVector(volvols_arg);
      this->hMatrix = new SGPP::base::DataMatrix(hMatrix_arg);
      this->r = r;

      bStochasticDataAlloc = true;
    }

    void HestonSolver::solveCrankNicolson(size_t numTimesteps, float_t timestepsize, size_t maxCGIterations, float_t epsilonCG, DataVector& alpha, size_t NumImEul) {
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

        std::cout << "Average Grid size: " << static_cast<float_t>(myHestonSystem->getSumGridPointsComplete()) / static_cast<float_t>(numTimesteps) << std::endl;
        std::cout << "Average Grid size (Inner): " << static_cast<float_t>(myHestonSystem->getSumGridPointsInner()) / static_cast<float_t>(numTimesteps) << std::endl << std::endl << std::endl;

        if (this->myScreen != NULL) {
          std::cout << "Time to solve: " << this->dNeededTime << " seconds" << std::endl;
          this->myScreen->writeEmptyLines(2);
        }

        this->finInnerGridSize = getNumberInnerGridPoints();
        this->avgInnerGridSize = static_cast<size_t>((static_cast<float_t>(myHestonSystem->getSumGridPointsInner()) / static_cast<float_t>(numTimesteps)) + 0.5);
        this->nNeededIterations = myEuler->getNumberIterations() + myCN->getNumberIterations();

        delete myHestonSystem;
        delete myCG;
        delete myCN;
        delete myEuler;
        delete myStopwatch;

        this->current_time += (static_cast<float_t>(numTimesteps) * timestepsize);
      } else {
        throw new application_exception("HestonSolver::solveCrankNicolson : A grid wasn't constructed before or stochastic parameters weren't set!");
      }
    }


    void HestonSolver::initGridWithPayoff(DataVector& alpha, float_t strike, std::string payoffType) {
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

    float_t HestonSolver::get1DEuroCallPayoffValue(float_t assetValue, float_t strike) {
      if (assetValue <= strike) {
        return 0.0;
      } else {
        return assetValue - strike;
      }
    }

    void HestonSolver::solveImplicitEuler(size_t numTimesteps, float_t timestepsize, size_t maxCGIterations, float_t epsilonCG, SGPP::base::DataVector& alpha, bool verbose, bool generateAnimation, size_t numEvalsAnimation) {
      throw new application_exception("This scheme is not implemented for the Heston solver!");
    }

    void HestonSolver::solveExplicitEuler(size_t numTimesteps, float_t timestepsize, size_t maxCGIterations, float_t epsilonCG, SGPP::base::DataVector& alpha, bool verbose, bool generateAnimation, size_t numEvalsAnimation) {
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
      this->myScreen->writeTitle("SGpp - Heston Solver, 1.0.0", "The SG++ Project (C) 2009-2010, by Sam Maurus (CSE Master's Thesis)");
      this->myScreen->writeStartSolve("Multidimensional Heston Solver");
    }

    void HestonSolver::setEnableCoarseningData(std::string adaptSolveMode, std::string refineMode, SGPP::base::GridIndex::level_type refineMaxLevel, int numCoarsenPoints, float_t coarsenThreshold, float_t refineThreshold) {
      this->useCoarsen = true;
      this->coarsenThreshold = coarsenThreshold;
      this->refineThreshold = refineThreshold;
      this->refineMaxLevel = refineMaxLevel;
      this->adaptSolveMode = adaptSolveMode;
      this->refineMode = refineMode;
      this->numCoarsenPoints = numCoarsenPoints;
    }

    size_t HestonSolver::getGridPointsAtMoney(std::string payoffType, float_t strike, float_t eps) {
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
              float_t stockPriceSum = 0.0;

              for (size_t j = 0; j < this->dim; j = j + 2) {
                stockPriceSum += coords[j];
              }

              if ( ((stockPriceSum / static_cast<float_t>(this->dim)) < (strike - eps)) || ((stockPriceSum / static_cast<float_t>(this->dim)) > (strike + eps)) ) {
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

    void HestonSolver::initCartesianGridWithPayoff(DataVector& alpha, float_t strike, std::string payoffType) {
      float_t tmp;

      //BlackScholesSolver* bsSolver = new BlackScholesSolver();

      if (this->bGridConstructed) {
        std::ofstream fileout;

        for (size_t i = 0; i < this->myGrid->getStorage()->size(); i++) {
          std::string coords = this->myGridStorage->get(i)->getCoordsStringBB(*this->myBoundingBox);

          GridIndex* curPoint = (*myGridStorage)[i];

          std::stringstream coordsStream(coords);
          float_t* dblFuncValues = new float_t[dim];

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
              float_t normalPayoff = std::max<float_t>(((tmp / static_cast<float_t>(numAssets)) - strike), 0.0);
              float_t vRange = this->myBoundingBox->getBoundary(1).rightBoundary - this->myBoundingBox->getBoundary(1).leftBoundary;
              float_t sPayoffDiff = dblFuncValues[0] - normalPayoff;
              alpha[i] = normalPayoff + ((dblFuncValues[1] - this->myBoundingBox->getBoundary(1).leftBoundary) / vRange) * sPayoffDiff;
            } else {
              // Payoff function
              alpha[i] = std::max<float_t>(((tmp / static_cast<float_t>(numAssets)) - strike), 0.0);
            }
          } else if (payoffType == "std_euro_put") {
            if (!curPoint->isInnerPoint()) {
              if (numAssets == 1 && dblFuncValues[1] == this->myBoundingBox->getBoundary(1).rightBoundary) {
                // Vmax boundary for a single asset. Strike price
                alpha[i] = strike;
              } else if (numAssets == 1 && dblFuncValues[0] == this->myBoundingBox->getBoundary(0).rightBoundary) {
                // Smax boundary for a single asset. exponential function
                //float_t constantC = 20.0;
                float_t constantB = (strike * pow(exp(this->myBoundingBox->getBoundary(1).leftBoundary - this->myBoundingBox->getBoundary(1).rightBoundary), 20.0)) / (1 - exp(this->myBoundingBox->getBoundary(1).leftBoundary - this->myBoundingBox->getBoundary(1).rightBoundary));
                float_t constantA = strike + constantB;
                alpha[i] = constantA * pow(exp(dblFuncValues[1] - this->myBoundingBox->getBoundary(1).rightBoundary), 20.0) - constantB;
              } else {

                // Get the payoff function value for this point
                float_t payoffFuncVal = 0.0;

                for (size_t k = 0; k < numAssets; k++) {
                  payoffFuncVal += dblFuncValues[2 * k];
                }

                payoffFuncVal = std::max<float_t>(strike - payoffFuncVal / (static_cast<float_t>(numAssets)), 0.0);

                // Get the max-volatility function value for this point
                float_t maxVolFuncVal = 0.0;

                for (size_t k = 0; k < numAssets; k++) {
                  float_t stockPrice = dblFuncValues[2 * k];
                  float_t sMaxValue = this->myBoundingBox->getBoundary(2 * k).rightBoundary;
                  maxVolFuncVal +=  (strike / (2.0 * sMaxValue)) * stockPrice;
                }

                maxVolFuncVal = strike - maxVolFuncVal / (static_cast<float_t>(numAssets));

                // Get the fraction that we are away from the v=vMin value in the direction of vMax
                float_t volFraction = 0.0;

                for (size_t k = 0; k < numAssets; k++) {
                  float_t vValue = dblFuncValues[2 * k + 1];
                  float_t vRightValue = this->myBoundingBox->getBoundary(2 * k + 1).rightBoundary;
                  float_t vLeftValue = this->myBoundingBox->getBoundary(2 * k + 1).leftBoundary;
                  float_t vRange = vRightValue - vLeftValue;
                  volFraction += (vValue - vLeftValue) / vRange;
                }

                volFraction = volFraction / (static_cast<float_t>(numAssets));

                // Now set the value of the point as the payoff function + fraction*maxVol function
                alpha[i] = std::max<float_t>(payoffFuncVal + volFraction * (maxVolFuncVal - payoffFuncVal), 0.0);
              }
            } else {
              // Non-boundary point. Just set the payoff value.
              tmp = 0.0;

              for (size_t j = 0; j < numAssets; j++) {
                tmp += dblFuncValues[2 * j];
              }

              alpha[i] = std::max<float_t>(strike - ((tmp / static_cast<float_t>(numAssets))), 0.0);
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
          SGPP::base::GridIndex* curPoint = (*myGridStorage)[i];

          if (curPoint->isInnerPoint() && alpha.get(i) != 0)
            numNonZeroInner++;
        }

        int k = 0;
        k++;

        OperationHierarchisation* myHierarchisation = SGPP::op_factory::createOperationHierarchisation(*this->myGrid);
        myHierarchisation->doHierarchisation(alpha);
        delete myHierarchisation;
      } else {
        throw new application_exception("HestonSolver::initCartesianGridWithPayoff : A grid wasn't constructed before!");
      }
    }

    void HestonSolver::initLogTransformedGridWithPayoff(DataVector& alpha, float_t strike, std::string payoffType) {
      float_t tmp;

      //BlackScholesSolver* bsSolver = new BlackScholesSolver();

      if (this->bGridConstructed) {
        for (size_t i = 0; i < this->myGrid->getStorage()->size(); i++) {
          std::string coords = this->myGridStorage->get(i)->getCoordsStringBB(*this->myBoundingBox);
          std::stringstream coordsStream(coords);
          float_t* dblFuncValues = new float_t[dim];

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
              float_t sumStockPrices = 0.0;
              float_t kMultiplier = 0;

              for (size_t k = 0; k < numAssets; k++) {
                // Get the sum of all stock prices
                sumStockPrices += exp(dblFuncValues[2 * k]);

                // At the same time, figure out the fraction of the way we are away from the max-v value for this stock
                float_t vValue = dblFuncValues[2 * k + 1];
                float_t vRightValue = this->myBoundingBox->getBoundary(2 * k + 1).rightBoundary;
                float_t vLeftValue = this->myBoundingBox->getBoundary(2 * k + 1).leftBoundary;
                float_t vRange = vRightValue - vLeftValue;
                kMultiplier += (vRightValue - vValue) / vRange;
              }

              alpha[i] = std::max<float_t>((sumStockPrices - strike * kMultiplier) / static_cast<float_t>(numAssets), 0.0);
            } else {
              // Non-boundary point. Just set the payoff value.
              tmp = 0.0;

              for (size_t j = 0; j < numAssets; j++) {
                tmp += exp(dblFuncValues[2 * j]);
              }

              alpha[i] = std::max<float_t>(((tmp / static_cast<float_t>(numAssets)) - strike), 0.0);
            }
          } else if (payoffType == "std_euro_put") {
            if (!curPoint->isInnerPoint()) {
              if (numAssets == 1 && dblFuncValues[1] == this->myBoundingBox->getBoundary(1).rightBoundary) {
                // Vmax boundary for a single asset. Strike price
                alpha[i] = strike;
              } else if (numAssets == 1 && dblFuncValues[0] == this->myBoundingBox->getBoundary(0).rightBoundary) {
                // Smax boundary for a single asset. exponential function
                //float_t constantC = 20.0;
                float_t constantB = (strike * pow(exp(this->myBoundingBox->getBoundary(1).leftBoundary - this->myBoundingBox->getBoundary(1).rightBoundary), 20.0)) / (1 - exp(this->myBoundingBox->getBoundary(1).leftBoundary - this->myBoundingBox->getBoundary(1).rightBoundary));
                float_t constantA = strike + constantB;
                alpha[i] = constantA * pow(exp(dblFuncValues[1] - this->myBoundingBox->getBoundary(1).rightBoundary), 20.0) - constantB;
              } else {
                // Get the payoff function value for this point
                float_t payoffFuncVal = 0.0;

                for (size_t k = 0; k < numAssets; k++) {
                  payoffFuncVal += exp(dblFuncValues[2 * k]);
                }

                payoffFuncVal = std::max<float_t>(strike - payoffFuncVal / (static_cast<float_t>(numAssets)), 0.0);

                // Get the max-volatility function value for this point
                float_t maxVolFuncVal = 0.0;

                for (size_t k = 0; k < numAssets; k++) {
                  float_t stockPrice = exp(dblFuncValues[2 * k]);
                  float_t sMinValue = exp(this->myBoundingBox->getBoundary(2 * k).leftBoundary);
                  float_t sMaxValue = exp(this->myBoundingBox->getBoundary(2 * k).rightBoundary);
                  maxVolFuncVal +=  (strike / (2.0 * (sMaxValue - sMinValue))) * (stockPrice - sMinValue);
                }

                maxVolFuncVal = strike - maxVolFuncVal / (static_cast<float_t>(numAssets));

                // Get the fraction that we are away from the v=vMin value in the direction of vMax
                float_t volFraction = 0.0;

                for (size_t k = 0; k < numAssets; k++) {
                  float_t vValue = dblFuncValues[2 * k + 1];
                  float_t vRightValue = this->myBoundingBox->getBoundary(2 * k + 1).rightBoundary;
                  float_t vLeftValue = this->myBoundingBox->getBoundary(2 * k + 1).leftBoundary;
                  float_t vRange = vRightValue - vLeftValue;
                  volFraction += (vValue - vLeftValue) / vRange;
                }

                volFraction = volFraction / (static_cast<float_t>(numAssets));

                // Now set the value of the point as the payoff function + fraction*maxVol function
                alpha[i] = std::max<float_t>(payoffFuncVal + volFraction * (maxVolFuncVal - payoffFuncVal), 0.0);
              }
            } else {
              // Non-boundary point. Just set the payoff value.
              tmp = 0.0;

              for (size_t j = 0; j < numAssets; j++) {
                tmp += exp(dblFuncValues[2 * j]);
              }

              alpha[i] = std::max<float_t>(strike - ((tmp / static_cast<float_t>(numAssets))), 0.0);
            }
          } else {
            throw new application_exception("HestonSolver::initLogTransformedGridWithPayoff : An unknown payoff-type was specified!");
          }

          delete[] dblFuncValues;
        }

        OperationHierarchisation* myHierarchisation = SGPP::op_factory::createOperationHierarchisation(*this->myGrid);
        myHierarchisation->doHierarchisation(alpha);
        delete myHierarchisation;
      } else {
        throw new application_exception("HestonSolver::initLogTransformedGridWithPayoff : A grid wasn't constructed before!");
      }
    }

    float_t HestonSolver::evalOption(std::vector<float_t>& eval_point, SGPP::base::DataVector& alpha) {
      std::vector<float_t> trans_eval = eval_point;

      // apply needed coordinate transformations
      if (this->useLogTransform) {
        for (size_t i = 0; i < eval_point.size(); i = i + 2) { // Here we know that the variance dimension is not log-transformed, so we shouldn't log its value
          trans_eval[i] = log(trans_eval[i]);
        }
      }

      SGPP::base::OperationEval* myEval = SGPP::op_factory::createOperationEval(*this->myGrid);
      float_t result = myEval->eval(alpha, trans_eval);
      delete myEval;

      return result;
    }

    void HestonSolver::transformPoint(SGPP::base::DataVector& point) {
      SGPP::base::DataVector tmp_point(point);

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

    float_t HestonSolver::getNeededTimeToSolve() {
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

    float_t EvaluateHestonClosedFormIntegralFunction(float_t phi, float_t xi, float_t theta, float_t kappa, float_t rho, float_t r, float_t T, float_t K, float_t S, float_t v, int type) {
      float_t a, b, u, x;

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
      complex<float_t> dBuilder1;
      complex<float_t> dBuilder2;
      complex<float_t> d;
      dBuilder1 = std::complex<float_t>(0.0, 1.0);
      dBuilder1 = pow(dBuilder1 * phi * rho * xi - b, 2.0);
      dBuilder2 = std::complex<float_t>(0.0, 1.0);
      dBuilder2 = dBuilder2 * phi;
      dBuilder2 = dBuilder2 * float_t(2.0) * u - pow(phi, float_t(2.0));
      dBuilder2 = dBuilder2 * float_t(pow(xi, 2.0));
      d = sqrt(dBuilder1 - dBuilder2);

      // Build g
      complex<float_t> g;
      g = std::complex<float_t>(0.0, 1.0);
      g = (b - rho * xi * phi * g + d) / (b - rho * xi * phi * g - d);


      // Build C
      complex<float_t> C;
      C = std::complex<float_t>(0.0, 1.0);
      C = (r * phi * C * T) + (a / float_t(pow(xi, 2.0))) * ((b - rho * xi * phi * C + d) * T - float_t(2.0) * log( (float_t(1.0) - g * exp(d * T)) / (float_t(1.0) - g)  ));

      // Build D
      complex<float_t> D;
      D = std::complex<float_t>(0.0, 1.0);
      D = ((b - rho * xi * phi * D + d) / float_t(pow(xi, 2.0))) * ((float_t(1.0) - exp(d * T)) / (float_t(1.0) - g * exp(d * T)));

      // Build f
      complex<float_t> f;
      f = complex<float_t>(0.0, 1.0);
      f = exp(C + D * v + f * phi * x);


      // Build realArgument
      complex<float_t> realArgument;
      realArgument = std::complex<float_t>(0.0, 1.0);
      realArgument = (exp(float_t(0.0) - realArgument * phi * float_t(log(K))) * f) / (realArgument * phi);

      return real(realArgument);
    }

    void HestonSolver::EvaluateHestonExactSurface(DataVector& alpha, float_t maturity) {
      if (!this->bGridConstructed)
        throw new application_exception("HestonSolver::EvaluateHestonPriceExact : The grid wasn't initialized before!");

      if (this->numAssets != 1 || this->payoffType != "std_euro_call")
        throw new application_exception("HestonSolver::EvaluateHestonPriceExact : Can only solve in closed form for a European call option with one asset!");

      float_t tmp;

      for (size_t i = 0; i < this->myGrid->getStorage()->size(); i++) {
        std::string coords = this->myGridStorage->get(i)->getCoordsStringBB(*this->myBoundingBox);
        std::stringstream coordsStream(coords);
        float_t* dblFuncValues = new float_t[dim];

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

      OperationHierarchisation* myHierarchisation = SGPP::op_factory::createOperationHierarchisation(*this->myGrid);
      myHierarchisation->doHierarchisation(alpha);
      delete myHierarchisation;
    }

    void HestonSolver::EvaluateHestonExactSurfacePut(DataVector& alpha, float_t maturity) {
      if (!this->bGridConstructed)
        throw new application_exception("HestonSolver::EvaluateHestonPriceExact : The grid wasn't initialized before!");

      if (this->numAssets != 1 || this->payoffType != "std_euro_put")
        throw new application_exception("HestonSolver::EvaluateHestonPriceExact : Can only solve in closed form for a European put option with one asset!");

      float_t tmp;

      for (size_t i = 0; i < this->myGrid->getStorage()->size(); i++) {
        std::string coords = this->myGridStorage->get(i)->getCoordsStringBB(*this->myBoundingBox);
        std::stringstream coordsStream(coords);
        float_t* dblFuncValues = new float_t[dim];

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

      OperationHierarchisation* myHierarchisation = SGPP::op_factory::createOperationHierarchisation(*this->myGrid);
      myHierarchisation->doHierarchisation(alpha);
      delete myHierarchisation;
    }

    void HestonSolver::CompareHestonBs1d(float_t maturity, float_t v) {
      if (this->numAssets != 1 || this->payoffType != "std_euro_call")
        throw new application_exception("HestonSolver::EvaluateHestonPriceExact : Can only solve in closed form for a European call option with one asset!");

      size_t dim1d = 1;
      int levels1d = this->levels;

      // Build a new 1d grid for stock price
      DimensionBoundary* boundaries1d = new SGPP::base::DimensionBoundary[dim1d];

      boundaries1d[0].leftBoundary = this->myBoundingBox->getBoundary(0).leftBoundary;
      boundaries1d[0].rightBoundary = this->myBoundingBox->getBoundary(0).rightBoundary;
      boundaries1d[0].bDirichletLeft = this->myBoundingBox->getBoundary(0).bDirichletLeft;
      boundaries1d[0].bDirichletRight = this->myBoundingBox->getBoundary(0).bDirichletRight;

      BoundingBox* boundingBox1d = new BoundingBox(dim1d, boundaries1d);

      Grid* grid1d = new LinearTruncatedBoundaryGrid(*boundingBox1d);

      GridGenerator* myGenerator = grid1d->createGridGenerator();
      myGenerator->regular(levels1d);

      SGPP::base::DataVector* alphaHeston = new SGPP::base::DataVector(grid1d->getSize());
      SGPP::base::DataVector* alphaBS = new SGPP::base::DataVector(grid1d->getSize());

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

    void HestonSolver::EvaluateHestonExact1d(DataVector& alpha, Grid* grid1d, BoundingBox* boundingBox1d, float_t maturity, float_t v) {
      float_t tmp;

      for (size_t i = 0; i < grid1d->getStorage()->size(); i++) {
        std::string coords = grid1d->getStorage()->get(i)->getCoordsStringBB(*boundingBox1d);
        std::stringstream coordsStream(coords);
        float_t* dblFuncValues = new float_t[1];
        coordsStream >> tmp;
        dblFuncValues[0] = tmp;

        alpha[i] = EvaluateHestonPriceExact(dblFuncValues[0], v, this->volvols->get(0), this->thetas->get(0), this->kappas->get(0), this->hMatrix->get(0, 1), this->r, maturity, this->dStrike) ;
        delete dblFuncValues;
      }

      OperationHierarchisation* myHierarchisation = SGPP::op_factory::createOperationHierarchisation(*grid1d);
      myHierarchisation->doHierarchisation(alpha);
      delete myHierarchisation;
    }

    void HestonSolver::EvaluateBsExact1d(DataVector& alpha, Grid* grid1d, BoundingBox* boundingBox1d, float_t maturity, float_t sigma) {
      float_t tmp;

      SGPP::finance::BlackScholesSolver* myBSSolver = new SGPP::finance::BlackScholesSolver(false);

      for (size_t i = 0; i < grid1d->getStorage()->size(); i++) {
        std::string coords = grid1d->getStorage()->get(i)->getCoordsStringBB(*boundingBox1d);
        std::stringstream coordsStream(coords);
        float_t* dblFuncValues = new float_t[1];
        coordsStream >> tmp;
        dblFuncValues[0] = tmp;

        alpha[i] = myBSSolver->getAnalyticSolution1D(dblFuncValues[0], true, maturity, sigma, this->r, this->dStrike);
        delete dblFuncValues;
      }

      OperationHierarchisation* myHierarchisation = SGPP::op_factory::createOperationHierarchisation(*grid1d);
      myHierarchisation->doHierarchisation(alpha);
      delete myHierarchisation;
      delete myBSSolver;
    }

    float_t HestonSolver::EvaluateHestonPriceExact(float_t S, float_t v, float_t xi, float_t theta, float_t kappa, float_t rho, float_t r, float_t T, float_t K) {
      float_t int1 = 0.5 + (1.0 / M_PI) * GaussLobattoInt(0.001, 1000.0, 1e-10, 100000, xi, theta, kappa, rho, r, T, K, S, v, 1);
      float_t int2 = 0.5 + (1.0 / M_PI) * GaussLobattoInt(0.001, 1000.0, 1e-10, 100000, xi, theta, kappa, rho, r, T, K, S, v, 2);
      return S * int1 - K * exp((-1.0) * r * T) * int2;
    }

    float_t HestonSolver::EvaluateHestonPriceExact(float_t S, float_t v, float_t maturity) {
      return EvaluateHestonPriceExact(S, v, this->volvols->get(0), this->thetas->get(0), this->kappas->get(0), this->hMatrix->get(0, 1), this->r, maturity, this->dStrike);
    }

    float_t HestonSolver::EvaluateHestonPriceExactPut(float_t S, float_t v, float_t xi, float_t theta, float_t kappa, float_t rho, float_t r, float_t T, float_t K) {
      float_t int1 = 0.5 + (1.0 / M_PI) * GaussLobattoInt(0.001, 1000.0, 1e-10, 100000, xi, theta, kappa, rho, r, T, K, S, v, 1);
      float_t int2 = 0.5 + (1.0 / M_PI) * GaussLobattoInt(0.001, 1000.0, 1e-10, 100000, xi, theta, kappa, rho, r, T, K, S, v, 2);
      return S * int1 - K * exp((-1.0) * r * T) * int2 + K * exp((-1.0) * r * T) - S;
    }

    float_t HestonSolver::EvaluateHestonPriceExactPut(float_t S, float_t v, float_t maturity) {
      return EvaluateHestonPriceExactPut(S, v, this->volvols->get(0), this->thetas->get(0), this->kappas->get(0), this->hMatrix->get(0, 1), this->r, maturity, this->dStrike);
    }

    void HestonSolver::GetBsExactSolution(SGPP::base::DataVector& alphaBS, float_t maturity) {
      SGPP::finance::BlackScholesSolver* myBSSolver = new SGPP::finance::BlackScholesSolver(false);

      float_t S, v;

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

      OperationHierarchisation* myHierarchisation = SGPP::op_factory::createOperationHierarchisation(*this->myGrid);
      myHierarchisation->doHierarchisation(alphaBS);
      delete myHierarchisation;
    }

    void HestonSolver::CompareHestonBsExact(SGPP::base::DataVector& alpha, float_t maturity) {
      DataVector* alphaHeston = new DataVector(getNumberGridPoints());
      DataVector* alphaBS = new DataVector(getNumberGridPoints());

      EvaluateHestonExactSurface(*alphaHeston, maturity);
      GetBsExactSolution(*alphaBS, maturity);

      // Find the difference (heston - BS)
      for (size_t i = 0; i < this->myGridStorage->size(); i++) {
        alpha[i] = (*alphaHeston)[i] - (*alphaBS)[i];
      }
    }

    void HestonSolver::CompareHestonNumericToBsExact(SGPP::base::DataVector& alphaHestonNumeric, SGPP::base::DataVector& alphaBS, SGPP::base::DataVector& error, float_t maturity) {
      GetBsExactSolution(alphaBS, maturity);

      // Find the difference (heston - BS)
      for (size_t i = 0; i < this->myGridStorage->size(); i++) {
        error[i] = alphaHestonNumeric[i] - alphaBS[i];
      }
    }


    float_t HestonSolver::GaussLobattoIntStep(
      float_t a, float_t b,
      float_t fa, float_t fb,
      size_t& neval,
      size_t maxeval,
      float_t acc
      , float_t xi, float_t theta, float_t kappa, float_t rho, float_t r, float_t T, float_t K, float_t S, float_t v, int type) {

      // Constants used in the algorithm
      const float_t alpha = std::sqrt(2.0 / 3.0);
      const float_t beta  = 1.0 / std::sqrt(5.0);

      if (neval >= maxeval) {
        throw new application_exception("HestonSolver::Gauss-Lobatto : Maximum number of evaluations reached in GaussLobatto.");
      }

      // Here the abcissa points and function values for both the 4-point
      // and the 7-point rule are calculated (the points at the end of
      // interval come from the function call, i.e., fa and fb. Also note
      // the 7-point rule re-uses all the points of the 4-point rule.)
      const float_t h = (b - a) / 2;
      const float_t m = (a + b) / 2;

      const float_t mll = m - alpha * h;
      const float_t ml = m - beta * h;
      const float_t mr = m + beta * h;
      const float_t mrr = m + alpha * h;

      const float_t fmll = EvaluateHestonClosedFormIntegralFunction(mll, xi, theta, kappa, rho, r, T, K, S, v, type);
      const float_t fml = EvaluateHestonClosedFormIntegralFunction(ml, xi, theta, kappa, rho, r, T, K, S, v, type);
      const float_t fm  = EvaluateHestonClosedFormIntegralFunction(m, xi, theta, kappa, rho, r, T, K, S, v, type);
      const float_t fmr = EvaluateHestonClosedFormIntegralFunction(mr, xi, theta, kappa, rho, r, T, K, S, v, type);
      const float_t fmrr = EvaluateHestonClosedFormIntegralFunction(mrr, xi, theta, kappa, rho, r, T, K, S, v, type);
      neval += 5;

      // Both the 4-point and 7-point rule integrals are evaluted
      const float_t integral2 = (h / 6) * (fa + fb + 5 * (fml + fmr));
      const float_t integral1 = (h / 1470) * (77 * (fa + fb)
                                              + 432 * (fmll + fmrr) + 625 * (fml + fmr) + 672 * fm);

      // The difference betwen the 4-point and 7-point integrals is the
      // estimate of the accuracy
      const float_t estacc = (integral1 - integral2);

      // The volatile keyword should prevent the floating point
      // destination value from being stored in extended precision
      // registers which actually have a very different
      // std::numeric_limits<float_t>::epsilon().
      volatile float_t dist = acc + estacc;

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

    float_t HestonSolver::GaussLobattoInt(float_t a, float_t b,
                                          float_t abstol,
                                          size_t maxeval
                                          , float_t xi, float_t theta, float_t kappa, float_t rho, float_t r, float_t T, float_t K, float_t S, float_t v, int type) {
      const float_t tol_epsunit = abstol / std::numeric_limits<float_t>::epsilon();
      size_t neval = 0;
      return GaussLobattoIntStep(a, b,
                                 EvaluateHestonClosedFormIntegralFunction(a, xi, theta, kappa, rho, r, T, K, S, v, type), EvaluateHestonClosedFormIntegralFunction(b, xi, theta, kappa, rho, r, T, K, S, v, type),
                                 neval,
                                 maxeval,
                                 tol_epsunit, xi, theta, kappa, rho, r, T, K, S, v, type);
    }

    void HestonSolver::CompareHestonSolutionToExact(SGPP::base::DataVector* solution, SGPP::base::DataVector* exact, std::string filename, size_t PointsPerDimension) {
      DimensionBoundary dimOne;
      DimensionBoundary dimTwo;
      std::ofstream fileout;

      fileout.open(filename.c_str());
      OperationEval* myEval = SGPP::op_factory::createOperationEval(*myGrid);


      if (myGrid->getStorage()->dim() == 2) {
        dimOne = myGrid->getBoundingBox()->getBoundary(0);
        dimTwo = myGrid->getBoundingBox()->getBoundary(1);

        float_t offset_x = dimOne.leftBoundary;
        float_t offset_y = dimTwo.leftBoundary;
        float_t inc_x = ((dimOne.rightBoundary - dimOne.leftBoundary) / float_t(PointsPerDimension - 1));
        float_t inc_y = ((dimTwo.rightBoundary - dimTwo.leftBoundary) / float_t(PointsPerDimension - 1));

        size_t points = (size_t)PointsPerDimension;

        for (size_t i = 0; i < points; i++) {
          for (size_t j = 0; j < points; j++) {
            std::vector<float_t> point;
            point.push_back(offset_x + (((float_t)(i))*inc_x));
            point.push_back(offset_y + (((float_t)(j))*inc_y));
            fileout << (offset_x + ((float_t)(i))*inc_x) << " " << (offset_y + ((float_t)(j))*inc_y) << " " << (myEval->eval(*solution, point) - myEval->eval(*exact, point)) << std::endl;
          }

          fileout << std::endl;
        }
      }

      delete myEval;
      // close filehandle
      fileout.close();
    }

    float_t HestonSolver::EvalSinglePoint1Asset(float_t s, float_t v, DataVector& alphaVec) {
      OperationEval* myEval = SGPP::op_factory::createOperationEval(*myGrid);
      std::vector<float_t> point;
      point.push_back(s);
      point.push_back(v);
      return myEval->eval(alphaVec, point);
    }

  }
}
