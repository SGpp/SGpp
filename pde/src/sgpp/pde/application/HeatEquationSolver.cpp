// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/pde/algorithm/HeatEquationParabolicPDESolverSystem.hpp>
#include <sgpp/pde/algorithm/HeatEquationParabolicPDESolverSystemParallelOMP.hpp>
#include <sgpp/pde/application/HeatEquationSolver.hpp>
#include <sgpp/solver/ode/Euler.hpp>
#include <sgpp/solver/ode/CrankNicolson.hpp>
#include <sgpp/solver/sle/ConjugateGradients.hpp>
#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/exception/application_exception.hpp>
#include <sgpp/base/operation/BaseOpFactory.hpp>

#include <sgpp/globaldef.hpp>

#include <stdlib.h>
#include <sstream>
#include <fstream>
#include <string>

namespace sgpp {
namespace pde {

HeatEquationSolver::HeatEquationSolver() : ParabolicPDESolver() {
  this->bGridConstructed = false;
  this->myScreen = nullptr;
}

HeatEquationSolver::~HeatEquationSolver() {
  if (this->myScreen != nullptr) {
    delete this->myScreen;
  }
}

void HeatEquationSolver::constructGrid(base::BoundingBox& BoundingBox, size_t level) {
  this->dim = BoundingBox.getDimension();
  this->levels = static_cast<int>(level);

  this->myGrid = new base::LinearBoundaryGrid(BoundingBox);

  this->myGrid->getGenerator().regular(levels);

  this->myBoundingBox = &this->myGrid->getBoundingBox();
  this->myGridStorage = &this->myGrid->getStorage();

  this->bGridConstructed = true;
}

void HeatEquationSolver::setHeatCoefficient(double a) { this->a = a; }

void HeatEquationSolver::solveExplicitEuler(size_t numTimesteps, double timestepsize,
                                            size_t maxCGIterations, double epsilonCG,
                                            base::DataVector& alpha, bool verbose,
                                            bool generateAnimation) {
  if (this->bGridConstructed) {
    this->myScreen->writeStartSolve("Multidimensional Heat Equation Solver");
    double dNeededTime;
    solver::Euler* myEuler = new solver::Euler(
        "ExEul", numTimesteps, timestepsize, generateAnimation, this->myScreen);
    solver::ConjugateGradients* myCG = new solver::ConjugateGradients(maxCGIterations, epsilonCG);
#ifdef _OPENMP
    HeatEquationParabolicPDESolverSystemParallelOMP* myHESolver =
        new HeatEquationParabolicPDESolverSystemParallelOMP(*this->myGrid, alpha, this->a,
                                                            timestepsize, "ExEul");
#else
    HeatEquationParabolicPDESolverSystem* myHESolver = new HeatEquationParabolicPDESolverSystem(
        *this->myGrid, alpha, this->a, timestepsize, "ExEul");
#endif
    base::SGppStopwatch* myStopwatch = new base::SGppStopwatch();

    myStopwatch->start();
    myEuler->solve(*myCG, *myHESolver, verbose);
    dNeededTime = myStopwatch->stop();

    if (this->myScreen != nullptr) {
      std::cout << "Time to solve: " << dNeededTime << " seconds" << std::endl;
      this->myScreen->writeEmptyLines(2);
    }

    delete myStopwatch;
    delete myCG;
    delete myEuler;
  } else {
    throw base::application_exception(
        "HeatEquationSolver::solveExplicitEuler : A grid wasn't constructed before!");
  }
}

void HeatEquationSolver::solveImplicitEuler(size_t numTimesteps, double timestepsize,
                                            size_t maxCGIterations, double epsilonCG,
                                            base::DataVector& alpha, bool verbose,
                                            bool generateAnimation) {
  if (this->bGridConstructed) {
    this->myScreen->writeStartSolve("Multidimensional Heat Equation Solver");
    double dNeededTime;
    solver::Euler* myEuler = new solver::Euler(
        "ImEul", numTimesteps, timestepsize, generateAnimation, this->myScreen);
    solver::ConjugateGradients* myCG = new solver::ConjugateGradients(maxCGIterations, epsilonCG);
#ifdef _OPENMP
    HeatEquationParabolicPDESolverSystemParallelOMP* myHESolver =
        new HeatEquationParabolicPDESolverSystemParallelOMP(*this->myGrid, alpha, this->a,
                                                            timestepsize, "ImEul");
#else
    HeatEquationParabolicPDESolverSystem* myHESolver = new HeatEquationParabolicPDESolverSystem(
        *this->myGrid, alpha, this->a, timestepsize, "ImEul");
#endif
    base::SGppStopwatch* myStopwatch = new base::SGppStopwatch();

    myStopwatch->start();
    myEuler->solve(*myCG, *myHESolver, verbose);
    dNeededTime = myStopwatch->stop();

    if (this->myScreen != nullptr) {
      std::cout << "Time to solve: " << dNeededTime << " seconds" << std::endl;
      this->myScreen->writeEmptyLines(2);
    }

    delete myStopwatch;
    delete myHESolver;
    delete myCG;
    delete myEuler;
  } else {
    throw base::application_exception(
        "HeatEquationSolver::solveImplicitEuler : A grid wasn't constructed before!");
  }
}

void HeatEquationSolver::solveCrankNicolson(size_t numTimesteps, double timestepsize,
                                            size_t maxCGIterations, double epsilonCG,
                                            base::DataVector& alpha, size_t NumImEul) {
  if (this->bGridConstructed) {
    this->myScreen->writeStartSolve("Multidimensional Heat Equation Solver");
    double dNeededTime;
    solver::ConjugateGradients* myCG = new solver::ConjugateGradients(maxCGIterations, epsilonCG);
#ifdef _OPENMP
    HeatEquationParabolicPDESolverSystemParallelOMP* myHESolver =
        new HeatEquationParabolicPDESolverSystemParallelOMP(*this->myGrid, alpha, this->a,
                                                            timestepsize, "CrNic");
#else
    HeatEquationParabolicPDESolverSystem* myHESolver = new HeatEquationParabolicPDESolverSystem(
        *this->myGrid, alpha, this->a, timestepsize, "CrNic");
#endif
    base::SGppStopwatch* myStopwatch = new base::SGppStopwatch();

    size_t numCNSteps;
    size_t numIESteps;

    numCNSteps = numTimesteps;

    if (numTimesteps > NumImEul) {
      numCNSteps = numTimesteps - NumImEul;
    }

    numIESteps = NumImEul;

    solver::Euler* myEuler =
        new solver::Euler("ImEul", numIESteps, timestepsize, false, this->myScreen);
    solver::CrankNicolson* myCN = new solver::CrankNicolson(numCNSteps, timestepsize);

    myStopwatch->start();

    if (numIESteps > 0) {
      myEuler->solve(*myCG, *myHESolver, false);
    }

    myCN->solve(*myCG, *myHESolver, false);
    dNeededTime = myStopwatch->stop();

    if (this->myScreen != nullptr) {
      std::cout << "Time to solve: " << dNeededTime << " seconds" << std::endl;
      this->myScreen->writeEmptyLines(2);
    }

    delete myStopwatch;
    delete myHESolver;
    delete myCG;
    delete myCN;
    delete myEuler;
  } else {
    throw base::application_exception(
        "HeatEquationSolver::solveCrankNicolson : A grid wasn't constructed before!");
  }
}

void HeatEquationSolver::initGridWithSmoothHeat(base::DataVector& alpha, double mu, double sigma,
                                                double factor) {
  if (this->bGridConstructed) {
    double tmp;
    double* dblFuncValues = new double[this->dim];

    for (size_t i = 0; i < this->myGrid->getSize(); i++) {
      std::string coords = this->myGridStorage->getCoordinates(
          this->myGridStorage->getPoint(i)).toString();
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

    sgpp::op_factory::createOperationHierarchisation(*this->myGrid)->doHierarchisation(alpha);
  } else {
    throw base::application_exception(
        "HeatEquationSolver::initGridWithSmoothHeat : A grid wasn't constructed before!");
  }
}

void HeatEquationSolver::initScreen() {
  this->myScreen = new base::ScreenOutput();
  this->myScreen->writeTitle("SGpp - Heat Equation Solver, 1.0.0",
                             "Alexander Heinecke, (C) 2009-2011");
}

void HeatEquationSolver::storeInnerRHS(base::DataVector& alpha, std::string tFilename,
                                       double timestepsize) {
  if (this->bGridConstructed) {
    HeatEquationParabolicPDESolverSystem* myHESolver = new HeatEquationParabolicPDESolverSystem(
        *this->myGrid, alpha, this->a, timestepsize, "ImEul");
    base::SGppStopwatch* myStopwatch = new base::SGppStopwatch();

    myStopwatch->start();
    std::cout << "Exporting inner right-hand-side..." << std::endl;
    base::DataVector* rhs_inner = myHESolver->generateRHS();

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

    delete myHESolver;
  } else {
    throw base::application_exception(
        "HeatEquationSolver::storeInnerMatrix : A grid wasn't constructed before!");
  }
}

void HeatEquationSolver::storeInnerSolution(base::DataVector& alpha, size_t numTimesteps,
                                            double timestepsize, size_t maxCGIterations,
                                            double epsilonCG, std::string tFilename) {
  if (this->bGridConstructed) {
    solver::Euler* myEuler =
        new solver::Euler("ImEul", numTimesteps, timestepsize, false, this->myScreen);
    solver::ConjugateGradients* myCG = new solver::ConjugateGradients(maxCGIterations, epsilonCG);
    HeatEquationParabolicPDESolverSystem* myHESolver = new HeatEquationParabolicPDESolverSystem(
        *this->myGrid, alpha, this->a, timestepsize, "ImEul");
    base::SGppStopwatch* myStopwatch = new base::SGppStopwatch();

    myStopwatch->start();
    std::cout << "Exporting inner solution..." << std::endl;
    myEuler->solve(*myCG, *myHESolver, false);

    base::DataVector* alpha_solve = myHESolver->getGridCoefficientsForCG();
    size_t nCoefs = alpha_solve->getSize();
    std::ofstream outfile(tFilename.c_str());

    for (size_t i = 0; i < nCoefs; i++) {
      outfile << std::scientific << alpha_solve->get(i) << std::endl;
    }

    outfile.close();

    std::cout << "Exporting inner solution... DONE!" << std::endl;

    delete myHESolver;
    delete myCG;
    delete myEuler;
  } else {
    throw base::application_exception(
        "HeatEquationSolver::solveImplicitEuler : A grid wasn't constructed before!");
  }
}
}  // namespace pde
}  // namespace sgpp
