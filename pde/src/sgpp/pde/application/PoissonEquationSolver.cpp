// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/pde/application/PoissonEquationSolver.hpp>
#include <sgpp/pde/algorithm/PoissonEquationEllipticPDESolverSystemDirichlet.hpp>
#include <sgpp/solver/sle/ConjugateGradients.hpp>
#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/exception/application_exception.hpp>
#include <sgpp/base/tools/SGppStopwatch.hpp>
#include <sgpp/base/operation/BaseOpFactory.hpp>

#include <sgpp/globaldef.hpp>

#include <stdlib.h>
#include <sstream>
#include <fstream>
#include <string>

namespace SGPP {
namespace pde {

PoissonEquationSolver::PoissonEquationSolver() : EllipticPDESolver() {
  this->bGridConstructed = false;
  this->myScreen = NULL;
}

PoissonEquationSolver::~PoissonEquationSolver() {
  if (this->myScreen != NULL) {
    delete this->myScreen;
  }
}

void PoissonEquationSolver::constructGrid(base::BoundingBox& BoundingBox, int level) {
  this->dim = BoundingBox.getDimensions();
  this->levels = level;

  this->myGrid = new base::LinearBoundaryGrid(BoundingBox);

  base::GridGenerator* myGenerator = this->myGrid->createGridGenerator();
  myGenerator->regular(this->levels);
  delete myGenerator;

  this->myBoundingBox = this->myGrid->getBoundingBox();
  this->myGridStorage = this->myGrid->getStorage();

  this->bGridConstructed = true;
}

void PoissonEquationSolver::solvePDE(base::DataVector& alpha, base::DataVector& rhs,
                                     size_t maxCGIterations, float_t epsilonCG, bool verbose) {
  float_t dTimeAlpha = 0.0;
  float_t dTimeRHS = 0.0;
  float_t dTimeSolver = 0.0;

  base::SGppStopwatch* myStopwatch = new base::SGppStopwatch();
  solver::ConjugateGradients* myCG = new solver::ConjugateGradients(maxCGIterations, epsilonCG);
  PoissonEquationEllipticPDESolverSystemDirichlet* mySystem =
      new PoissonEquationEllipticPDESolverSystemDirichlet(*(this->myGrid), rhs);

  std::cout << "Gridpoints (complete grid): " << mySystem->getNumGridPointsComplete() << std::endl;
  std::cout << "Gridpoints (inner grid): " << mySystem->getNumGridPointsInner() << std::endl
            << std::endl
            << std::endl;

  myStopwatch->start();
  base::DataVector* alpha_solve = mySystem->getGridCoefficientsForCG();
  dTimeAlpha = myStopwatch->stop();
  std::cout << "coefficients has been initialized for solving!" << std::endl;
  myStopwatch->start();
  base::DataVector* rhs_solve = mySystem->generateRHS();
  dTimeRHS = myStopwatch->stop();
  std::cout << "right hand side has been initialized for solving!" << std::endl
            << std::endl
            << std::endl;

  myStopwatch->start();
  myCG->solve(*mySystem, *alpha_solve, *rhs_solve, true, verbose, 0.0);

  // Copy result into coefficient vector of the boundary grid
  mySystem->getSolutionBoundGrid(alpha, *alpha_solve);
  dTimeSolver = myStopwatch->stop();

  std::cout << std::endl << std::endl;
  std::cout << "Gridpoints (complete grid): " << mySystem->getNumGridPointsComplete() << std::endl;
  std::cout << "Gridpoints (inner grid): " << mySystem->getNumGridPointsInner() << std::endl
            << std::endl
            << std::endl;

  std::cout << "Timings for solving Poisson Equation" << std::endl;
  std::cout << "------------------------------------" << std::endl;
  std::cout << "Time for creating CG coeffs: " << dTimeAlpha << std::endl;
  std::cout << "Time for creating RHS: " << dTimeRHS << std::endl;
  std::cout << "Time for solving: " << dTimeSolver << std::endl << std::endl;
  std::cout << "Time: " << dTimeAlpha + dTimeRHS + dTimeSolver << std::endl
            << std::endl
            << std::endl;

  delete myCG;
  delete mySystem;  // alpha_solver and rhs_solve are allocated and freed here!!
  delete myStopwatch;
}

void PoissonEquationSolver::initGridWithSmoothHeat(base::DataVector& alpha, float_t mu,
                                                   float_t sigma, float_t factor) {
  if (this->bGridConstructed) {
    float_t tmp;
    float_t* dblFuncValues = new float_t[this->dim];

    for (size_t i = 0; i < this->myGrid->getStorage()->size(); i++) {
      std::string coords = this->myGridStorage->get(i)->getCoordsStringBB(*this->myBoundingBox);
      std::stringstream coordsStream(coords);
      bool isInner = true;

      for (size_t j = 0; j < this->dim; j++) {
        coordsStream >> tmp;

        // determine if a grid point is an inner grid point
        if ((tmp != this->myBoundingBox->getBoundary(j).leftBoundary &&
             tmp != this->myBoundingBox->getBoundary(j).rightBoundary)) {
          // Nothtin to do, test is that qay hence == for floating point values is unsave
        } else {
          isInner = false;
        }

        dblFuncValues[j] = tmp;
      }

      if (isInner == false) {
        tmp = 1.0;

        for (size_t j = 0; j < this->dim; j++) {
          tmp *= factor * factor *
                 ((1.0 / (sigma * 2.0 * 3.145)) * exp((-0.5) * ((dblFuncValues[j] - mu) / sigma) *
                                                      ((dblFuncValues[j] - mu) / sigma)));
        }
      } else {
        tmp = 0.0;
      }

      alpha[i] = tmp;
    }

    delete[] dblFuncValues;

    base::OperationHierarchisation* myHierarchisation =
        SGPP::op_factory::createOperationHierarchisation(*this->myGrid);
    myHierarchisation->doHierarchisation(alpha);
    delete myHierarchisation;
  } else {
    throw new base::application_exception(
        "HeatEquationSolver::initGridWithSmoothHeat : A grid wasn't constructed before!");
  }
}

void PoissonEquationSolver::initGridWithSmoothHeatFullDomain(base::DataVector& alpha, float_t mu,
                                                             float_t sigma, float_t factor) {
  if (this->bGridConstructed) {
    float_t tmp;
    float_t* dblFuncValues = new float_t[this->dim];

    for (size_t i = 0; i < this->myGrid->getStorage()->size(); i++) {
      std::string coords = this->myGridStorage->get(i)->getCoordsStringBB(*this->myBoundingBox);
      std::stringstream coordsStream(coords);

      for (size_t j = 0; j < this->dim; j++) {
        coordsStream >> tmp;

        dblFuncValues[j] = tmp;
      }

      tmp = 1.0;

      for (size_t j = 0; j < this->dim; j++) {
        tmp *= factor * factor *
               ((1.0 / (sigma * 2.0 * 3.145)) * exp((-0.5) * ((dblFuncValues[j] - mu) / sigma) *
                                                    ((dblFuncValues[j] - mu) / sigma)));
      }

      alpha[i] = tmp;
    }

    delete[] dblFuncValues;

    base::OperationHierarchisation* myHierarchisation =
        SGPP::op_factory::createOperationHierarchisation(*this->myGrid);
    myHierarchisation->doHierarchisation(alpha);
    delete myHierarchisation;
  } else {
    throw new base::application_exception(
        "HeatEquationSolver::initGridWithSmoothHeatFullDomain : A grid wasn't constructed before!");
  }
}

void PoissonEquationSolver::initGridWithExpHeat(base::DataVector& alpha, float_t factor) {
  if (this->bGridConstructed) {
    float_t tmp;
    float_t* dblFuncValues = new float_t[this->dim];
    float_t* rightBound = new float_t[this->dim];

    base::BoundingBox* tmpBB = this->myGrid->getBoundingBox();

    for (size_t j = 0; j < this->dim; j++) {
      rightBound[j] = (tmpBB->getBoundary(j)).rightBoundary;
    }

    for (size_t i = 0; i < this->myGrid->getStorage()->size(); i++) {
      std::string coords = this->myGridStorage->get(i)->getCoordsStringBB(*this->myBoundingBox);
      std::stringstream coordsStream(coords);
      bool isInner = true;
      tmp = 0.0;

      for (size_t j = 0; j < this->dim; j++) {
        coordsStream >> tmp;

        // determine if a grid point is an inner grid point
        if ((tmp != this->myBoundingBox->getBoundary(j).leftBoundary &&
             tmp != this->myBoundingBox->getBoundary(j).rightBoundary)) {
          // Nothtin to do, test is that qay hence == for floating point values is unsave
        } else {
          isInner = false;
        }

        dblFuncValues[j] = tmp;
      }

      if (isInner == false) {
        tmp = 1.0;

        for (size_t j = 0; j < this->dim; j++) {
          tmp *= exp((dblFuncValues[j] - rightBound[j]) * factor);
        }
      } else {
        tmp = 0.0;
      }

      alpha[i] = tmp;
    }

    delete[] dblFuncValues;

    base::OperationHierarchisation* myHierarchisation =
        SGPP::op_factory::createOperationHierarchisation(*this->myGrid);
    myHierarchisation->doHierarchisation(alpha);
    delete myHierarchisation;
  } else {
    throw new base::application_exception(
        "PoissonEquationSolver::initGridWithExpHeat : A grid wasn't constructed before!");
  }
}

void PoissonEquationSolver::initGridWithExpHeatFullDomain(base::DataVector& alpha, float_t factor) {
  if (this->bGridConstructed) {
    float_t tmp;
    float_t* dblFuncValues = new float_t[this->dim];
    float_t* rightBound = new float_t[this->dim];

    base::BoundingBox* tmpBB = this->myGrid->getBoundingBox();

    for (size_t j = 0; j < this->dim; j++) {
      rightBound[j] = (tmpBB->getBoundary(j)).rightBoundary;
    }

    for (size_t i = 0; i < this->myGrid->getStorage()->size(); i++) {
      std::string coords = this->myGridStorage->get(i)->getCoordsStringBB(*this->myBoundingBox);
      std::stringstream coordsStream(coords);
      tmp = 0.0;

      for (size_t j = 0; j < this->dim; j++) {
        coordsStream >> tmp;

        dblFuncValues[j] = tmp;
      }

      tmp = 1.0;

      for (size_t j = 0; j < this->dim; j++) {
        tmp *= exp((dblFuncValues[j] - rightBound[j]) * factor);
      }

      alpha[i] = tmp;
    }

    delete[] dblFuncValues;

    base::OperationHierarchisation* myHierarchisation =
        SGPP::op_factory::createOperationHierarchisation(*this->myGrid);
    myHierarchisation->doHierarchisation(alpha);
    delete myHierarchisation;
  } else {
    throw new base::application_exception(
        "PoissonEquationSolver::initGridWithExpHeat : A grid wasn't constructed before!");
  }
}

void PoissonEquationSolver::storeInnerRHS(base::DataVector& alpha, std::string tFilename) {
  base::SGppStopwatch* myStopwatch = new base::SGppStopwatch();
  PoissonEquationEllipticPDESolverSystemDirichlet* mySystem =
      new PoissonEquationEllipticPDESolverSystemDirichlet(*(this->myGrid), alpha);

  std::cout << "Exporting inner right-hand-side..." << std::endl;
  myStopwatch->start();
  base::DataVector* rhs_inner = mySystem->generateRHS();

  size_t nCoefs = rhs_inner->getSize();
  std::ofstream outfile(tFilename.c_str());

  for (size_t i = 0; i < nCoefs; i++) {
    outfile << std::scientific << rhs_inner->get(i) << std::endl;
  }

  outfile.close();
  std::cout << "Exporting inner right-hand-side... DONE! (" << myStopwatch->stop() << " s)"
            << std::endl
            << std::endl
            << std::endl;

  delete mySystem;  // rhs_inner are allocated and freed here!!
  delete myStopwatch;
}

void PoissonEquationSolver::storeInnerSolution(base::DataVector& alpha, size_t maxCGIterations,
                                               float_t epsilonCG, std::string tFilename) {
  solver::ConjugateGradients* myCG = new solver::ConjugateGradients(maxCGIterations, epsilonCG);
  PoissonEquationEllipticPDESolverSystemDirichlet* mySystem =
      new PoissonEquationEllipticPDESolverSystemDirichlet(*(this->myGrid), alpha);

  std::cout << "Exporting inner solution..." << std::endl;

  base::DataVector* alpha_solve = mySystem->getGridCoefficientsForCG();
  base::DataVector* rhs_solve = mySystem->generateRHS();

  myCG->solve(*mySystem, *alpha_solve, *rhs_solve, true, false, 0.0);

  size_t nCoefs = alpha_solve->getSize();
  std::ofstream outfile(tFilename.c_str());

  for (size_t i = 0; i < nCoefs; i++) {
    outfile << std::scientific << alpha_solve->get(i) << std::endl;
  }

  outfile.close();

  std::cout << "Exporting inner solution... DONE!" << std::endl;

  delete myCG;
  delete mySystem;  // alpha_solver and rhs_solve are allocated and freed here!!
}

void PoissonEquationSolver::initScreen() {
  this->myScreen = new base::ScreenOutput();
  this->myScreen->writeTitle("SGpp - Poisson Equation Solver, 1.0.0",
                             "Alexander Heinecke, (C) 2009-2011");
}
}  // namespace pde
}  // namespace SGPP
