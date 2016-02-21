// Copyright (C) 2008-today The SG++ Project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/finance/algorithm/BlackScholesParabolicPDESolverSystem.hpp>
#include <sgpp/finance/algorithm/BlackScholesParabolicPDESolverSystemEuroAmer.hpp>
#include <sgpp/finance/algorithm/BlackScholesParabolicPDESolverSystemEuroAmerParallelOMP.hpp>
#include <sgpp/finance/application/BlackScholesSolverWithStretching.hpp>
#include <sgpp/solver/ode/Euler.hpp>
#include <sgpp/solver/ode/CrankNicolson.hpp>
#include <sgpp/solver/ode/StepsizeControlH.hpp>
#include <sgpp/solver/ode/StepsizeControlBDF.hpp>
#include <sgpp/solver/ode/StepsizeControlEJ.hpp>
#include <sgpp/solver/sle/BiCGStab.hpp>
#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/exception/application_exception.hpp>
#include <sgpp/base/operation/BaseOpFactory.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>

#include <sgpp/globaldef.hpp>

#include <cstdlib>
#include <sstream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <string>
#include <vector>

namespace SGPP {
namespace finance {

BlackScholesSolverWithStretching::BlackScholesSolverWithStretching(bool useLogTransform,
                                                                   std::string OptionType)
    : BlackScholesSolver(useLogTransform) {
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

  this->myStretching = NULL;
}

BlackScholesSolverWithStretching::~BlackScholesSolverWithStretching() {}

void BlackScholesSolverWithStretching::getGridNormalDistribution(SGPP::base::DataVector& alpha,
                                                                 std::vector<float_t>& norm_mu,
                                                                 std::vector<float_t>& norm_sigma) {
  if (this->bGridConstructed) {
    float_t tmp;
    float_t value;
    SGPP::base::StdNormalDistribution myNormDistr;

    for (size_t i = 0; i < this->myGrid->getStorage().size(); i++) {
      std::string coords =
          this->myGridStorage->get(i)->getCoordsStringStretching(*(this->myStretching));
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
    throw SGPP::base::application_exception(
        "BlackScholesSolverWithStretching::getGridNormalDistribution : The grid wasn't initialized "
        "before!");
  }
}

void BlackScholesSolverWithStretching::constructGridStretching(SGPP::base::Stretching& stretching,
                                                               int level) {
  this->dim = stretching.getDimensions();
  this->levels = level;

  this->myGrid = new SGPP::base::LinearStretchedBoundaryGrid(stretching);

  this->myGrid->getGenerator().regular(this->levels);

  this->myStretching = &this->myGrid->getStretching();
  this->myGridStorage = &this->myGrid->getStorage();

  // std::string serGrid;
  // myGrid->serialize(serGrid);
  // std::cout << serGrid << std::endl;

  this->bGridConstructed = true;
}

void BlackScholesSolverWithStretching::constructGrid(SGPP::base::BoundingBox& myBoundingBox,
                                                     size_t level) {
  throw SGPP::base::application_exception(
      "BlackScholesSolverWithStretching::constructGrid : This solver does not support "
      "SGPP::base::BoundingBox, use constructGridStretching instead!");
}

void BlackScholesSolverWithStretching::refineInitialGridWithPayoff(SGPP::base::DataVector& alpha,
                                                                   float_t strike,
                                                                   std::string payoffType,
                                                                   float_t dStrikeDistance) {
  size_t nRefinements = 0;

  this->dStrike = strike;
  this->payoffType = payoffType;

  if (this->useLogTransform == false) {
    if (this->bGridConstructed) {
      SGPP::base::DataVector refineVector(alpha.getSize());

      if (payoffType == "std_euro_call" || payoffType == "std_euro_put") {
        this->tBoundaryType = "Dirichlet";

        float_t tmp;
        float_t* dblFuncValues = new float_t[dim];
        float_t dDistance = 0.0;

        for (size_t i = 0; i < this->myGrid->getStorage().size(); i++) {
          std::string coords =
              this->myGridStorage->get(i)->getCoordsStringStretching(*(this->myStretching));
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

          if (payoffType == "std_euro_put") {
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

        SGPP::base::SurplusRefinementFunctor myRefineFunc(refineVector, nRefinements, 0.0);
        this->myGrid->getGenerator().refine(myRefineFunc);

        alpha.resize(this->myGridStorage->size());

        // reinit the grid with the payoff function
        initGridWithPayoff(alpha, strike, payoffType);
      } else {
        throw SGPP::base::application_exception(
            "BlackScholesSolverWithStretching::refineInitialGridWithPayoff : An unsupported "
            "payoffType was specified!");
      }
    } else {
      throw SGPP::base::application_exception(
          "BlackScholesSolverWithStretching::refineInitialGridWithPayoff : The grid wasn't "
          "initialized before!");
    }
  }
}

void BlackScholesSolverWithStretching::refineInitialGridWithPayoffToMaxLevel(
    SGPP::base::DataVector& alpha, float_t strike, std::string payoffType, float_t dStrikeDistance,
    SGPP::base::GridIndex::level_type maxLevel) {
  size_t nRefinements = 0;

  this->dStrike = strike;
  this->payoffType = payoffType;

  if (this->useLogTransform == false) {
    if (this->bGridConstructed) {
      SGPP::base::DataVector refineVector(alpha.getSize());

      if (payoffType == "std_euro_call" || payoffType == "std_euro_put") {
        this->tBoundaryType = "Dirichlet";

        float_t tmp;
        float_t* dblFuncValues = new float_t[dim];
        float_t dDistance = 0.0;

        for (size_t i = 0; i < this->myGrid->getStorage().size(); i++) {
          std::string coords =
              this->myGridStorage->get(i)->getCoordsStringStretching(*this->myStretching);
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

          if (payoffType == "std_euro_put") {
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

        SGPP::base::SurplusRefinementFunctor myRefineFunc(refineVector, nRefinements, 0.0);
        this->myGrid->getGenerator().refineMaxLevel(myRefineFunc, maxLevel);

        alpha.resize(this->myGridStorage->size());

        // reinit the grid with the payoff function
        initGridWithPayoff(alpha, strike, payoffType);
      } else {
        throw SGPP::base::application_exception(
            "BlackScholesSolverWithStretching::refineInitialGridWithPayoffToMaxLevel : An "
            "unsupported payoffType was specified!");
      }
    } else {
      throw SGPP::base::application_exception(
          "BlackScholesSolverWithStretching::refineInitialGridWithPayoffToMaxLevel : The grid "
          "wasn't initialized before!");
    }
  }
}

void BlackScholesSolverWithStretching::initGridWithPayoff(SGPP::base::DataVector& alpha,
                                                          float_t strike, std::string payoffType) {
  this->dStrike = strike;
  this->payoffType = payoffType;

  if (payoffType == "std_euro_call" || payoffType == "std_euro_put" ||
      payoffType == "std_amer_put") {
    this->tBoundaryType = "Dirichlet";
  }

  if (this->useLogTransform) {
    initLogTransformedGridWithPayoff(alpha, strike, payoffType);
  } else {
    initCartesianGridWithPayoff(alpha, strike, payoffType);
  }
}

void BlackScholesSolverWithStretching::initScreen() {
  this->myScreen = new SGPP::base::ScreenOutput();
  this->myScreen->writeTitle(
      "SGpp - Black Scholes Solver with SGPP::base::Stretching, 2.1.0",
      "The SG++ Project (C) 2009-2010, by Alexander Heinecke and Sarpkan Selcuk");
  this->myScreen->writeStartSolve(
      "Multidimensional Black Scholes Solver with SGPP::base::Stretching");
}

void BlackScholesSolverWithStretching::printPayoffInterpolationError2D(
    SGPP::base::DataVector& alpha, std::string tFilename, size_t numTestpoints, float_t strike) {
  if (this->useLogTransform == false) {
    if (this->bGridConstructed) {
      if (this->myGrid->getStorage().getStretching()->getDimensions() == 2) {
        if (numTestpoints < 2) numTestpoints = 2;

        float_t dInc = (2.0 * strike) / static_cast<float_t>(numTestpoints - 1);

        float_t dX = 0.0;
        float_t dY = 2 * strike;

        std::ofstream file;
        file.open(tFilename.c_str());

        std::unique_ptr<SGPP::base::OperationEval> myEval(
            SGPP::op_factory::createOperationEval(*this->myGrid));

        for (size_t i = 0; i < numTestpoints; i++) {
          std::vector<float_t> point;

          point.push_back(dX);
          point.push_back(dY);

          float_t result = myEval->eval(alpha, point);

          file << std::scientific << std::setprecision(16) << dX << " " << dY << " " << result
               << std::endl;

          dX += dInc;
          dY -= dInc;
        }

        file.close();
      }
    } else {
      throw SGPP::base::application_exception(
          "BlackScholesSolverWithStretching::getPayoffInterpolationError : A grid wasn't "
          "constructed before!");
    }
  }
}

size_t BlackScholesSolverWithStretching::getGridPointsAtMoney(std::string payoffType,
                                                              float_t strike, float_t eps) {
  size_t nPoints = 0;

  if (this->useLogTransform == false) {
    if (this->bGridConstructed) {
      for (size_t i = 0; i < this->myGrid->getStorage().size(); i++) {
        bool isAtMoney = true;
        SGPP::base::DataVector coords(this->dim);
        this->myGridStorage->get(i)->getCoordsStretching(coords, *this->myStretching);

        if (payoffType == "std_euro_call" || payoffType == "std_euro_put") {
          for (size_t d = 0; d < this->dim; d++) {
            if (((coords.sum() / static_cast<float_t>(this->dim)) < (strike - eps)) ||
                ((coords.sum() / static_cast<float_t>(this->dim)) > (strike + eps))) {
              isAtMoney = false;
            }
          }
        } else {
          throw SGPP::base::application_exception(
              "BlackScholesSolverWithStretching::getGridPointsAtMoney : An unknown payoff-type was "
              "specified!");
        }

        if (isAtMoney == true) {
          nPoints++;
        }
      }
    } else {
      throw SGPP::base::application_exception(
          "BlackScholesSolverWithStretching::getGridPointsAtMoney : A grid wasn't constructed "
          "before!");
    }
  }

  return nPoints;
}

void BlackScholesSolverWithStretching::initCartesianGridWithPayoff(SGPP::base::DataVector& alpha,
                                                                   float_t strike,
                                                                   std::string payoffType) {
  float_t tmp;

  if (this->bGridConstructed) {
    for (size_t i = 0; i < this->myGrid->getStorage().size(); i++) {
      std::string coords =
          this->myGridStorage->get(i)->getCoordsStringStretching(*this->myStretching);
      std::stringstream coordsStream(coords);
      float_t* dblFuncValues = new float_t[dim];

      for (size_t j = 0; j < this->dim; j++) {
        coordsStream >> tmp;

        dblFuncValues[j] = tmp;
      }

      if (payoffType == "std_euro_call") {
        tmp = 0.0;

        for (size_t j = 0; j < dim; j++) {
          tmp += dblFuncValues[j];
        }

        alpha[i] = std::max<float_t>(((tmp / static_cast<float_t>(dim)) - strike), 0.0);
      } else if (payoffType == "std_euro_put") {
        tmp = 0.0;

        for (size_t j = 0; j < dim; j++) {
          tmp += dblFuncValues[j];
        }

        alpha[i] = std::max<float_t>(strike - ((tmp / static_cast<float_t>(dim))), 0.0);
      } else {
        throw SGPP::base::application_exception(
            "BlackScholesSolverWithStretching::initCartesianGridWithPayoff : An unknown "
            "payoff-type was specified!");
      }

      delete[] dblFuncValues;
    }

    SGPP::op_factory::createOperationHierarchisation(*this->myGrid)->doHierarchisation(alpha);
  } else {
    throw SGPP::base::application_exception(
        "BlackScholesSolverWithStretching::initCartesianGridWithPayoff : A grid wasn't constructed "
        "before!");
  }
}

void BlackScholesSolverWithStretching::initLogTransformedGridWithPayoff(
    SGPP::base::DataVector& alpha, float_t strike, std::string payoffType) {
  float_t tmp;

  if (this->bGridConstructed) {
    for (size_t i = 0; i < this->myGrid->getStorage().size(); i++) {
      std::string coords =
          this->myGridStorage->get(i)->getCoordsStringStretching(*this->myStretching);
      std::stringstream coordsStream(coords);
      float_t* dblFuncValues = new float_t[dim];

      for (size_t j = 0; j < this->dim; j++) {
        coordsStream >> tmp;

        dblFuncValues[j] = tmp;
      }

      if (payoffType == "std_euro_call") {
        tmp = 0.0;

        for (size_t j = 0; j < dim; j++) {
          tmp += exp(dblFuncValues[j]);
        }

        alpha[i] = std::max<float_t>(((tmp / static_cast<float_t>(dim)) - strike), 0.0);
      } else if (payoffType == "std_euro_put") {
        tmp = 0.0;

        for (size_t j = 0; j < dim; j++) {
          tmp += exp(dblFuncValues[j]);
        }

        alpha[i] = std::max<float_t>(strike - ((tmp / static_cast<float_t>(dim))), 0.0);
      } else {
        throw SGPP::base::application_exception(
            "BlackScholesSolverWithStretching::initLogTransformedGridWithPayoff : An unknown "
            "payoff-type was specified!");
      }

      delete[] dblFuncValues;
    }

    SGPP::op_factory::createOperationHierarchisation(*this->myGrid)->doHierarchisation(alpha);
  } else {
    throw SGPP::base::application_exception(
        "BlackScholesSolverWithStretching::initLogTransformedGridWithPayoff : A grid wasn't "
        "constructed before!");
  }
}

void BlackScholesSolverWithStretching::getAnalyticAlpha1D(base::DataVector& alpha_analytic,
                                                          float_t strike, float_t t,
                                                          std::string payoffType,
                                                          bool hierarchized) {
  float_t coord;

  if (dim != 1) {
    throw base::application_exception(
        "BlackScholesSolver::getAnalyticAlpha1D : A grid wasn't constructed before!");
  }

  if (!this->bGridConstructed) {
    throw base::application_exception(
        "BlackScholesSolver::getAnalyticAlpha1D : function only available for dim = 1!");
  }

  // compute values of analytic solution on given grid
  for (size_t i = 0; i < this->myGridStorage->size(); i++) {
    std::string coords =
        this->myGridStorage->get(i)->getCoordsStringStretching(*this->myStretching);
    std::stringstream coordsStream(coords);
    coordsStream >> coord;

    if (useLogTransform) {
      coord = exp(coord);
    }

    if (payoffType == "std_euro_call") {
      alpha_analytic[i] =
          this->getAnalyticSolution1D(coord, true, t, this->sigmas->get(0), this->r, strike);
    } else if (payoffType == "std_euro_put") {
      alpha_analytic[i] =
          this->getAnalyticSolution1D(coord, false, t, this->sigmas->get(0), this->r, strike);
    }
  }

  if (hierarchized) {
    // hierarchize computed values
    SGPP::op_factory::createOperationHierarchisation(*this->myGrid)->
        doHierarchisation(alpha_analytic);
  }
}

void BlackScholesSolverWithStretching::printGrid(SGPP::base::DataVector& alpha,
                                                 size_t PointesPerDimension,
                                                 std::string tfilename) const {
  SGPP::base::GridPrinterForStretching myPrinter(*this->myGrid);
  myPrinter.printGrid(alpha, tfilename, PointesPerDimension);
}

void BlackScholesSolverWithStretching::printGridDomainStretching(SGPP::base::DataVector& alpha,
                                                                 size_t PointesPerDimension,
                                                                 SGPP::base::Stretching& GridArea,
                                                                 std::string tfilename) const {
  SGPP::base::GridPrinterForStretching myPrinter(*this->myGrid);
  myPrinter.printGridDomainStretching(alpha, tfilename, GridArea, PointesPerDimension);
}

void BlackScholesSolverWithStretching::printGridDomain(SGPP::base::DataVector& alpha,
                                                       size_t PointesPerDimension,
                                                       SGPP::base::BoundingBox& GridArea,
                                                       std::string tfilename) const {
  throw SGPP::base::application_exception(
      "BlackScholesSolverWithStretching::printGridDomain: SGPP::base::BoundingBox not supported, "
      "use printGridDomainStretching instead!");
}

void BlackScholesSolverWithStretching::printSparseGrid(SGPP::base::DataVector& alpha,
                                                       std::string tfilename, bool bSurplus) const {
  SGPP::base::GridPrinterForStretching myPrinter(*this->myGrid);
  myPrinter.printSparseGrid(alpha, tfilename, bSurplus);
}

void BlackScholesSolverWithStretching::printSparseGridExpTransform(SGPP::base::DataVector& alpha,
                                                                   std::string tfilename,
                                                                   bool bSurplus) const {
  SGPP::base::GridPrinterForStretching myPrinter(*this->myGrid);
  myPrinter.printSparseGridExpTransform(alpha, tfilename, bSurplus);
}
}  // namespace finance
}  // namespace SGPP
