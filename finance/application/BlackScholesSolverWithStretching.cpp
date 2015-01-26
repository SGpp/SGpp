/* ****************************************************************************
* Copyright (C) 2011 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Alexander Heinecke (Alexander.Heinecke@mytum.de), Sarpkan Selcuk (Sarpkan.Selcuk@mytum.de)

#include "finance/algorithm/BlackScholesParabolicPDESolverSystem.hpp"
#include "finance/algorithm/BlackScholesParabolicPDESolverSystemEuroAmer.hpp"
#include "finance/algorithm/BlackScholesParabolicPDESolverSystemEuroAmerParallelOMP.hpp"
#include "finance/application/BlackScholesSolverWithStretching.hpp"
#include "solver/ode/Euler.hpp"
#include "solver/ode/CrankNicolson.hpp"
#include "solver/ode/StepsizeControlH.hpp"
#include "solver/ode/StepsizeControlBDF.hpp"
#include "solver/ode/StepsizeControlEJ.hpp"
#include "solver/sle/BiCGStab.hpp"
#include "base/grid/Grid.hpp"
#include "base/exception/application_exception.hpp"
#include "base/operation/BaseOpFactory.hpp"
#include "base/datatypes/DataVector.hpp"
#include <cstdlib>
#include <sstream>
#include <cmath>
#include <fstream>
#include <iomanip>

namespace sg {
  namespace finance {

    BlackScholesSolverWithStretching::BlackScholesSolverWithStretching(bool useLogTransform, std::string OptionType) : BlackScholesSolver(useLogTransform) {
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
      //  std::cout<<"BSSolverWithStretching\n";
      this->tBoundaryType = "freeBoundaries";

      // @todo set to random value to remove compiler warnings due to uninitialized members
      this->myStretching = NULL;
    }

    BlackScholesSolverWithStretching::~BlackScholesSolverWithStretching() {
    }

    void BlackScholesSolverWithStretching::getGridNormalDistribution(sg::base::DataVector& alpha, std::vector<double>& norm_mu, std::vector<double>& norm_sigma) {
      if (this->bGridConstructed) {
        double tmp;
        double value;
        sg::base::StdNormalDistribution myNormDistr;

        for (size_t i = 0; i < this->myGrid->getStorage()->size(); i++) {
          std::string coords = this->myGridStorage->get(i)->getCoordsStringStretching(*(this->myStretching));
          std::stringstream coordsStream(coords);

          value = 1.0;

          for (size_t j = 0; j < this->dim; j++) {
            coordsStream >> tmp;

            if (this->useLogTransform == false) {
              value *= myNormDistr.getDensity(tmp, norm_mu[j], norm_sigma[j]);
            } else {
              value *= myNormDistr.getDensity(exp(tmp), norm_mu[j], norm_sigma[j]);
            }
          }

          alpha[i] = value;
        }
      } else {
        throw new sg::base::application_exception("BlackScholesSolverWithStretching::getGridNormalDistribution : The grid wasn't initialized before!");
      }
    }


    void BlackScholesSolverWithStretching::constructGridStretching(sg::base::Stretching& stretching, int level) {
      this->dim = stretching.getDimensions();
      this->levels = level;

      this->myGrid = new sg::base::LinearStretchedTrapezoidBoundaryGrid(stretching);

      sg::base::GridGenerator* myGenerator = this->myGrid->createGridGenerator();
      myGenerator->regular(this->levels);
      delete myGenerator;

      this->myStretching = this->myGrid->getStretching();
      this->myGridStorage = this->myGrid->getStorage();

      //std::string serGrid;
      //myGrid->serialize(serGrid);
      //std::cout << serGrid << std::endl;

      this->bGridConstructed = true;
    }

    void BlackScholesSolverWithStretching::constructGrid(sg::base::BoundingBox& myBoundingBox, size_t level) {
      throw new sg::base::application_exception("BlackScholesSolverWithStretching::constructGrid : This solver does not support sg::base::BoundingBox, use constructGridStretching instead!");
    }

    void BlackScholesSolverWithStretching::refineInitialGridWithPayoff(sg::base::DataVector& alpha, double strike, std::string payoffType, double dStrikeDistance) {
      size_t nRefinements = 0;

      this->dStrike = strike;
      this->payoffType = payoffType;

      if (this->useLogTransform == false) {
        if (this->bGridConstructed) {

          sg::base::DataVector refineVector(alpha.getSize());

          if (payoffType == "std_euro_call" || payoffType == "std_euro_put") {
            this->tBoundaryType = "Dirichlet";

            double tmp;
            double* dblFuncValues = new double[dim];
            double dDistance = 0.0;

            for (size_t i = 0; i < this->myGrid->getStorage()->size(); i++) {
              std::string coords = this->myGridStorage->get(i)->getCoordsStringStretching(*(this->myStretching));
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

              if (payoffType == "std_euro_put") {
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

            sg::base::SurplusRefinementFunctor* myRefineFunc = new sg::base::SurplusRefinementFunctor(&refineVector, nRefinements, 0.0);

            this->myGrid->createGridGenerator()->refine(myRefineFunc);

            delete myRefineFunc;

            alpha.resize(this->myGridStorage->size());

            // reinit the grid with the payoff function
            initGridWithPayoff(alpha, strike, payoffType);
          } else {
            throw new sg::base::application_exception("BlackScholesSolverWithStretching::refineInitialGridWithPayoff : An unsupported payoffType was specified!");
          }
        } else {
          throw new sg::base::application_exception("BlackScholesSolverWithStretching::refineInitialGridWithPayoff : The grid wasn't initialized before!");
        }
      }
    }

    void BlackScholesSolverWithStretching::refineInitialGridWithPayoffToMaxLevel(sg::base::DataVector& alpha, double strike, std::string payoffType, double dStrikeDistance, sg::base::GridIndex::level_type maxLevel) {
      size_t nRefinements = 0;

      this->dStrike = strike;
      this->payoffType = payoffType;

      if (this->useLogTransform == false) {
        if (this->bGridConstructed) {

          sg::base::DataVector refineVector(alpha.getSize());

          if (payoffType == "std_euro_call" || payoffType == "std_euro_put") {
            this->tBoundaryType = "Dirichlet";

            double tmp;
            double* dblFuncValues = new double[dim];
            double dDistance = 0.0;

            for (size_t i = 0; i < this->myGrid->getStorage()->size(); i++) {
              std::string coords = this->myGridStorage->get(i)->getCoordsStringStretching(*this->myStretching);
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

              if (payoffType == "std_euro_put") {
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

            sg::base::SurplusRefinementFunctor* myRefineFunc = new sg::base::SurplusRefinementFunctor(&refineVector, nRefinements, 0.0);

            this->myGrid->createGridGenerator()->refineMaxLevel(myRefineFunc, maxLevel);

            delete myRefineFunc;

            alpha.resize(this->myGridStorage->size());

            // reinit the grid with the payoff function
            initGridWithPayoff(alpha, strike, payoffType);
          } else {
            throw new sg::base::application_exception("BlackScholesSolverWithStretching::refineInitialGridWithPayoffToMaxLevel : An unsupported payoffType was specified!");
          }
        } else {
          throw new sg::base::application_exception("BlackScholesSolverWithStretching::refineInitialGridWithPayoffToMaxLevel : The grid wasn't initialized before!");
        }
      }
    }


    void BlackScholesSolverWithStretching::initGridWithPayoff(sg::base::DataVector& alpha, double strike, std::string payoffType) {
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

    void BlackScholesSolverWithStretching::initScreen() {
      this->myScreen = new sg::base::ScreenOutput();
      this->myScreen->writeTitle("SGpp - Black Scholes Solver with sg::base::Stretching, 2.1.0", "TUM (C) 2009-2010, by Alexander Heinecke and Sarpkan Selcuk");
      this->myScreen->writeStartSolve("Multidimensional Black Scholes Solver with sg::base::Stretching");
    }


    void BlackScholesSolverWithStretching::printPayoffInterpolationError2D(sg::base::DataVector& alpha, std::string tFilename, size_t numTestpoints, double strike) {
      if (this->useLogTransform == false) {
        if (this->bGridConstructed) {
          if (this->myGrid->getStorage()->getStretching()->getDimensions() == 2) {
            if (numTestpoints < 2)
              numTestpoints = 2;

            double dInc = (2.0 * strike) / static_cast<double>(numTestpoints - 1);

            double dX = 0.0;
            double dY = 2 * strike;

            std::ofstream file;
            file.open(tFilename.c_str());

            sg::base::OperationEval* myEval = sg::op_factory::createOperationEval(*this->myGrid);

            for (size_t i = 0; i < numTestpoints; i++) {
              std::vector<double> point;

              point.push_back(dX);
              point.push_back(dY);

              double result = myEval->eval(alpha, point);

              file << std::scientific << std::setprecision( 16 ) << dX << " " << dY << " " << result << std::endl;

              dX += dInc;
              dY -= dInc;
            }

            delete myEval;

            file.close();
          }
        } else {
          throw new sg::base::application_exception("BlackScholesSolverWithStretching::getPayoffInterpolationError : A grid wasn't constructed before!");
        }
      }
    }

    size_t BlackScholesSolverWithStretching::getGridPointsAtMoney(std::string payoffType, double strike, double eps) {
      size_t nPoints = 0;

      if (this->useLogTransform == false) {
        if (this->bGridConstructed) {
          for (size_t i = 0; i < this->myGrid->getStorage()->size(); i++) {
            bool isAtMoney = true;
            sg::base::DataVector coords(this->dim);
            this->myGridStorage->get(i)->getCoordsStretching(coords, *this->myStretching);

            if (payoffType == "std_euro_call" || payoffType == "std_euro_put") {
              for (size_t d = 0; d < this->dim; d++) {
                if ( ((coords.sum() / static_cast<double>(this->dim)) < (strike - eps)) || ((coords.sum() / static_cast<double>(this->dim)) > (strike + eps)) ) {
                  isAtMoney = false;
                }

              }
            } else {
              throw new sg::base::application_exception("BlackScholesSolverWithStretching::getGridPointsAtMoney : An unknown payoff-type was specified!");
            }

            if (isAtMoney == true) {
              nPoints++;
            }
          }
        } else {
          throw new sg::base::application_exception("BlackScholesSolverWithStretching::getGridPointsAtMoney : A grid wasn't constructed before!");
        }
      }

      return nPoints;
    }

    void BlackScholesSolverWithStretching::initCartesianGridWithPayoff(sg::base::DataVector& alpha, double strike, std::string payoffType) {
      double tmp;

      if (this->bGridConstructed) {
        for (size_t i = 0; i < this->myGrid->getStorage()->size(); i++) {
          std::string coords = this->myGridStorage->get(i)->getCoordsStringStretching(*this->myStretching);
          std::stringstream coordsStream(coords);
          double* dblFuncValues = new double[dim];

          for (size_t j = 0; j < this->dim; j++) {
            coordsStream >> tmp;

            dblFuncValues[j] = tmp;
          }

          if (payoffType == "std_euro_call") {
            tmp = 0.0;

            for (size_t j = 0; j < dim; j++) {
              tmp += dblFuncValues[j];
            }

            alpha[i] = std::max<double>(((tmp / static_cast<double>(dim)) - strike), 0.0);
          } else if (payoffType == "std_euro_put") {
            tmp = 0.0;

            for (size_t j = 0; j < dim; j++) {
              tmp += dblFuncValues[j];
            }

            alpha[i] = std::max<double>(strike - ((tmp / static_cast<double>(dim))), 0.0);
          } else {
            throw new sg::base::application_exception("BlackScholesSolverWithStretching::initCartesianGridWithPayoff : An unknown payoff-type was specified!");
          }

          delete[] dblFuncValues;
        }

        sg::base::OperationHierarchisation* myHierarchisation = sg::op_factory::createOperationHierarchisation(*this->myGrid);
        myHierarchisation->doHierarchisation(alpha);
        delete myHierarchisation;
      } else {
        throw new sg::base::application_exception("BlackScholesSolverWithStretching::initCartesianGridWithPayoff : A grid wasn't constructed before!");
      }
    }

    void BlackScholesSolverWithStretching::initLogTransformedGridWithPayoff(sg::base::DataVector& alpha, double strike, std::string payoffType) {
      double tmp;

      if (this->bGridConstructed) {
        for (size_t i = 0; i < this->myGrid->getStorage()->size(); i++) {
          std::string coords = this->myGridStorage->get(i)->getCoordsStringStretching(*this->myStretching);
          std::stringstream coordsStream(coords);
          double* dblFuncValues = new double[dim];

          for (size_t j = 0; j < this->dim; j++) {
            coordsStream >> tmp;

            dblFuncValues[j] = tmp;
          }

          if (payoffType == "std_euro_call") {
            tmp = 0.0;

            for (size_t j = 0; j < dim; j++) {
              tmp += exp(dblFuncValues[j]);
            }

            alpha[i] = std::max<double>(((tmp / static_cast<double>(dim)) - strike), 0.0);
          } else if (payoffType == "std_euro_put") {
            tmp = 0.0;

            for (size_t j = 0; j < dim; j++) {
              tmp += exp(dblFuncValues[j]);
            }

            alpha[i] = std::max<double>(strike - ((tmp / static_cast<double>(dim))), 0.0);
          } else {
            throw new sg::base::application_exception("BlackScholesSolverWithStretching::initLogTransformedGridWithPayoff : An unknown payoff-type was specified!");
          }

          delete[] dblFuncValues;
        }

        sg::base::OperationHierarchisation* myHierarchisation = sg::op_factory::createOperationHierarchisation(*this->myGrid);
        myHierarchisation->doHierarchisation(alpha);
        delete myHierarchisation;
      } else {
        throw new sg::base::application_exception("BlackScholesSolverWithStretching::initLogTransformedGridWithPayoff : A grid wasn't constructed before!");
      }
    }

    void BlackScholesSolverWithStretching::getAnalyticAlpha1D(base::DataVector& alpha_analytic, double strike, double t, std::string payoffType, bool hierarchized) {
      double coord;

      if (dim != 1) {
        throw new base::application_exception("BlackScholesSolver::getAnalyticAlpha1D : A grid wasn't constructed before!");
      }

      if (!this->bGridConstructed) {
        throw new base::application_exception("BlackScholesSolver::getAnalyticAlpha1D : function only available for dim = 1!");
      }

      // compute values of analytic solution on given grid
      for (size_t i = 0; i < this->myGridStorage->size(); i++) {
        std::string coords = this->myGridStorage->get(i)->getCoordsStringStretching(*this->myStretching);
        std::stringstream coordsStream(coords);
        coordsStream >> coord;

        if (useLogTransform) {
          coord = exp(coord);
        }

        if (payoffType == "std_euro_call") {
          alpha_analytic[i] = this->getAnalyticSolution1D(coord, true, t, this->sigmas->get(0), this->r, strike);
        } else if (payoffType == "std_euro_put") {
          alpha_analytic[i] = this->getAnalyticSolution1D(coord, false, t, this->sigmas->get(0), this->r, strike);
        }
      }

      if (hierarchized) {
        // hierarchize computed values
        base::OperationHierarchisation* myHier = sg::op_factory::createOperationHierarchisation(*this->myGrid);
        myHier->doHierarchisation(alpha_analytic);

        delete myHier;
      }
    }

    void BlackScholesSolverWithStretching::printGrid(sg::base::DataVector& alpha, size_t PointesPerDimension, std::string tfilename) const {
      sg::base::GridPrinterForStretching myPrinter(*this->myGrid);
      myPrinter.printGrid(alpha, tfilename, PointesPerDimension);
    }

    void BlackScholesSolverWithStretching::printGridDomainStretching(sg::base::DataVector& alpha, size_t PointesPerDimension, sg::base::Stretching& GridArea, std::string tfilename) const {
      sg::base::GridPrinterForStretching myPrinter(*this->myGrid);
      myPrinter.printGridDomainStretching(alpha, tfilename, GridArea, PointesPerDimension);
    }

    void BlackScholesSolverWithStretching::printGridDomain(sg::base::DataVector& alpha, size_t PointesPerDimension, sg::base::BoundingBox& GridArea, std::string tfilename)const {
      throw new sg::base::application_exception("BlackScholesSolverWithStretching::printGridDomain: sg::base::BoundingBox not supported, use printGridDomainStretching instead!");
    }

    void BlackScholesSolverWithStretching::printSparseGrid(sg::base::DataVector& alpha, std::string tfilename, bool bSurplus) const {
      sg::base::GridPrinterForStretching myPrinter(*this->myGrid);
      myPrinter.printSparseGrid(alpha, tfilename, bSurplus);
    }

    void BlackScholesSolverWithStretching::printSparseGridExpTransform(sg::base::DataVector& alpha, std::string tfilename, bool bSurplus) const {
      sg::base::GridPrinterForStretching myPrinter(*this->myGrid);
      myPrinter.printSparseGridExpTransform(alpha, tfilename, bSurplus);
    }

  }
}
