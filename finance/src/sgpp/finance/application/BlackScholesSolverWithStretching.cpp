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

namespace sgpp {
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

void BlackScholesSolverWithStretching::getGridNormalDistribution(sgpp::base::DataVector& alpha,
                                                                 std::vector<double>& norm_mu,
                                                                 std::vector<double>& norm_sigma) {
  if (this->bGridConstructed) {
    double tmp;
    double value;
    sgpp::base::StdNormalDistribution myNormDistr;

    for (size_t i = 0; i < this->myGrid->getSize(); i++) {
      std::string coords =
          this->myGridStorage->getGridIndex(i).getCoordsStringStretching(*(this->myStretching));
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
    throw sgpp::base::application_exception(
        "BlackScholesSolverWithStretching::getGridNormalDistribution : The grid wasn't initialized "
        "before!");
  }
}

void BlackScholesSolverWithStretching::constructGridStretching(sgpp::base::Stretching& stretching,
                                                               int level) {
  this->dim = stretching.getDimensions();
  this->levels = level;

  this->myGrid = new sgpp::base::LinearStretchedBoundaryGrid(stretching);

  this->myGrid->getGenerator().regular(this->levels);

  this->myStretching = &this->myGrid->getStretching();
  this->myGridStorage = &this->myGrid->getStorage();

  // std::string serGrid;
  // myGrid->serialize(serGrid);
  // std::cout << serGrid << std::endl;

  this->bGridConstructed = true;
}

void BlackScholesSolverWithStretching::constructGrid(sgpp::base::BoundingBox& myBoundingBox,
                                                     size_t level) {
  throw sgpp::base::application_exception(
      "BlackScholesSolverWithStretching::constructGrid : This solver does not support "
      "sgpp::base::BoundingBox, use constructGridStretching instead!");
}

void BlackScholesSolverWithStretching::refineInitialGridWithPayoff(sgpp::base::DataVector& alpha,
                                                                   double strike,
                                                                   std::string payoffType,
                                                                   double dStrikeDistance) {
  size_t nRefinements = 0;

  this->dStrike = strike;
  this->payoffType = payoffType;

  if (this->useLogTransform == false) {
    if (this->bGridConstructed) {
      sgpp::base::DataVector refineVector(alpha.getSize());

      if (payoffType == "std_euro_call" || payoffType == "std_euro_put") {
        this->tBoundaryType = "Dirichlet";

        double tmp;
        double* dblFuncValues = new double[dim];
        double dDistance = 0.0;

        for (size_t i = 0; i < this->myGrid->getSize(); i++) {
          std::string coords =
              this->myGridStorage->getGridIndex(i).getCoordsStringStretching(*(this->myStretching));
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

        sgpp::base::SurplusRefinementFunctor myRefineFunc(refineVector, nRefinements, 0.0);
        this->myGrid->getGenerator().refine(myRefineFunc);

        alpha.resize(this->myGridStorage->getSize());

        // reinit the grid with the payoff function
        initGridWithPayoff(alpha, strike, payoffType);
      } else {
        throw sgpp::base::application_exception(
            "BlackScholesSolverWithStretching::refineInitialGridWithPayoff : An unsupported "
            "payoffType was specified!");
      }
    } else {
      throw sgpp::base::application_exception(
          "BlackScholesSolverWithStretching::refineInitialGridWithPayoff : The grid wasn't "
          "initialized before!");
    }
  }
}

void BlackScholesSolverWithStretching::refineInitialGridWithPayoffToMaxLevel(
    sgpp::base::DataVector& alpha, double strike, std::string payoffType, double dStrikeDistance,
    sgpp::base::GridIndex::level_type maxLevel) {
  size_t nRefinements = 0;

  this->dStrike = strike;
  this->payoffType = payoffType;

  if (this->useLogTransform == false) {
    if (this->bGridConstructed) {
      sgpp::base::DataVector refineVector(alpha.getSize());

      if (payoffType == "std_euro_call" || payoffType == "std_euro_put") {
        this->tBoundaryType = "Dirichlet";

        double tmp;
        double* dblFuncValues = new double[dim];
        double dDistance = 0.0;

        for (size_t i = 0; i < this->myGrid->getSize(); i++) {
          std::string coords =
              this->myGridStorage->getGridIndex(i).getCoordsStringStretching(*this->myStretching);
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

        sgpp::base::SurplusRefinementFunctor myRefineFunc(refineVector, nRefinements, 0.0);
        this->myGrid->getGenerator().refineMaxLevel(myRefineFunc, maxLevel);

        alpha.resize(this->myGridStorage->getSize());

        // reinit the grid with the payoff function
        initGridWithPayoff(alpha, strike, payoffType);
      } else {
        throw sgpp::base::application_exception(
            "BlackScholesSolverWithStretching::refineInitialGridWithPayoffToMaxLevel : An "
            "unsupported payoffType was specified!");
      }
    } else {
      throw sgpp::base::application_exception(
          "BlackScholesSolverWithStretching::refineInitialGridWithPayoffToMaxLevel : The grid "
          "wasn't initialized before!");
    }
  }
}

void BlackScholesSolverWithStretching::initGridWithPayoff(sgpp::base::DataVector& alpha,
                                                          double strike, std::string payoffType) {
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
  this->myScreen = new sgpp::base::ScreenOutput();
  this->myScreen->writeTitle(
      "SGpp - Black Scholes Solver with sgpp::base::Stretching, 2.1.0",
      "The SG++ Project (C) 2009-2010, by Alexander Heinecke and Sarpkan Selcuk");
  this->myScreen->writeStartSolve(
      "Multidimensional Black Scholes Solver with sgpp::base::Stretching");
}

void BlackScholesSolverWithStretching::printPayoffInterpolationError2D(
    sgpp::base::DataVector& alpha, std::string tFilename, size_t numTestpoints, double strike) {
  if (this->useLogTransform == false) {
    if (this->bGridConstructed) {
      if (this->myGrid->getStorage().getStretching()->getDimensions() == 2) {
        if (numTestpoints < 2) numTestpoints = 2;

        double dInc = (2.0 * strike) / static_cast<double>(numTestpoints - 1);

        double dX = 0.0;
        double dY = 2 * strike;

        std::ofstream file;
        file.open(tFilename.c_str());

        std::unique_ptr<sgpp::base::OperationEval> myEval(
            sgpp::op_factory::createOperationEval(*this->myGrid));

        for (size_t i = 0; i < numTestpoints; i++) {
          std::vector<double> point;

          point.push_back(dX);
          point.push_back(dY);

          double result = myEval->eval(alpha, point);

          file << std::scientific << std::setprecision(16) << dX << " " << dY << " " << result
               << std::endl;

          dX += dInc;
          dY -= dInc;
        }

        file.close();
      }
    } else {
      throw sgpp::base::application_exception(
          "BlackScholesSolverWithStretching::getPayoffInterpolationError : A grid wasn't "
          "constructed before!");
    }
  }
}

size_t BlackScholesSolverWithStretching::getGridPointsAtMoney(std::string payoffType,
                                                              double strike, double eps) {
  size_t nPoints = 0;

  if (this->useLogTransform == false) {
    if (this->bGridConstructed) {
      for (size_t i = 0; i < this->myGrid->getSize(); i++) {
        bool isAtMoney = true;
        sgpp::base::DataVector coords(this->dim);
        this->myGridStorage->getGridIndex(i).getCoordsStretching(coords, *this->myStretching);

        if (payoffType == "std_euro_call" || payoffType == "std_euro_put") {
          for (size_t d = 0; d < this->dim; d++) {
            if (((coords.sum() / static_cast<double>(this->dim)) < (strike - eps)) ||
                ((coords.sum() / static_cast<double>(this->dim)) > (strike + eps))) {
              isAtMoney = false;
            }
          }
        } else {
          throw sgpp::base::application_exception(
              "BlackScholesSolverWithStretching::getGridPointsAtMoney : An unknown payoff-type was "
              "specified!");
        }

        if (isAtMoney == true) {
          nPoints++;
        }
      }
    } else {
      throw sgpp::base::application_exception(
          "BlackScholesSolverWithStretching::getGridPointsAtMoney : A grid wasn't constructed "
          "before!");
    }
  }

  return nPoints;
}

void BlackScholesSolverWithStretching::initCartesianGridWithPayoff(sgpp::base::DataVector& alpha,
                                                                   double strike,
                                                                   std::string payoffType) {
  double tmp;

  if (this->bGridConstructed) {
    for (size_t i = 0; i < this->myGrid->getSize(); i++) {
      std::string coords =
          this->myGridStorage->getGridIndex(i).getCoordsStringStretching(*this->myStretching);
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
        throw sgpp::base::application_exception(
            "BlackScholesSolverWithStretching::initCartesianGridWithPayoff : An unknown "
            "payoff-type was specified!");
      }

      delete[] dblFuncValues;
    }

    sgpp::op_factory::createOperationHierarchisation(*this->myGrid)->doHierarchisation(alpha);
  } else {
    throw sgpp::base::application_exception(
        "BlackScholesSolverWithStretching::initCartesianGridWithPayoff : A grid wasn't constructed "
        "before!");
  }
}

void BlackScholesSolverWithStretching::initLogTransformedGridWithPayoff(
    sgpp::base::DataVector& alpha, double strike, std::string payoffType) {
  double tmp;

  if (this->bGridConstructed) {
    for (size_t i = 0; i < this->myGrid->getSize(); i++) {
      std::string coords =
          this->myGridStorage->getGridIndex(i).getCoordsStringStretching(*this->myStretching);
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
        throw sgpp::base::application_exception(
            "BlackScholesSolverWithStretching::initLogTransformedGridWithPayoff : An unknown "
            "payoff-type was specified!");
      }

      delete[] dblFuncValues;
    }

    sgpp::op_factory::createOperationHierarchisation(*this->myGrid)->doHierarchisation(alpha);
  } else {
    throw sgpp::base::application_exception(
        "BlackScholesSolverWithStretching::initLogTransformedGridWithPayoff : A grid wasn't "
        "constructed before!");
  }
}

void BlackScholesSolverWithStretching::getAnalyticAlpha1D(base::DataVector& alpha_analytic,
                                                          double strike, double t,
                                                          std::string payoffType,
                                                          bool hierarchized) {
  double coord;

  if (dim != 1) {
    throw base::application_exception(
        "BlackScholesSolver::getAnalyticAlpha1D : A grid wasn't constructed before!");
  }

  if (!this->bGridConstructed) {
    throw base::application_exception(
        "BlackScholesSolver::getAnalyticAlpha1D : function only available for dim = 1!");
  }

  // compute values of analytic solution on given grid
  for (size_t i = 0; i < this->myGridStorage->getSize(); i++) {
    std::string coords =
        this->myGridStorage->getGridIndex(i).getCoordsStringStretching(*this->myStretching);
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
    sgpp::op_factory::createOperationHierarchisation(*this->myGrid)->
        doHierarchisation(alpha_analytic);
  }
}

void BlackScholesSolverWithStretching::printGrid(sgpp::base::DataVector& alpha,
                                                 size_t PointesPerDimension,
                                                 std::string tfilename) const {
  sgpp::base::GridPrinterForStretching myPrinter(*this->myGrid);
  myPrinter.printGrid(alpha, tfilename, PointesPerDimension);
}

void BlackScholesSolverWithStretching::printGridDomainStretching(sgpp::base::DataVector& alpha,
                                                                 size_t PointesPerDimension,
                                                                 sgpp::base::Stretching& GridArea,
                                                                 std::string tfilename) const {
  sgpp::base::GridPrinterForStretching myPrinter(*this->myGrid);
  myPrinter.printGridDomainStretching(alpha, tfilename, GridArea, PointesPerDimension);
}

void BlackScholesSolverWithStretching::printGridDomain(sgpp::base::DataVector& alpha,
                                                       size_t PointesPerDimension,
                                                       sgpp::base::BoundingBox& GridArea,
                                                       std::string tfilename) const {
  throw sgpp::base::application_exception(
      "BlackScholesSolverWithStretching::printGridDomain: sgpp::base::BoundingBox not supported, "
      "use printGridDomainStretching instead!");
}

void BlackScholesSolverWithStretching::printSparseGrid(sgpp::base::DataVector& alpha,
                                                       std::string tfilename, bool bSurplus) const {
  sgpp::base::GridPrinterForStretching myPrinter(*this->myGrid);
  myPrinter.printSparseGrid(alpha, tfilename, bSurplus);
}

void BlackScholesSolverWithStretching::printSparseGridExpTransform(sgpp::base::DataVector& alpha,
                                                                   std::string tfilename,
                                                                   bool bSurplus) const {
  sgpp::base::GridPrinterForStretching myPrinter(*this->myGrid);
  myPrinter.printSparseGridExpTransform(alpha, tfilename, bSurplus);
}
}  // namespace finance
}  // namespace sgpp
