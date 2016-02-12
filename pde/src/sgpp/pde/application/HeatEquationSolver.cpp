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

namespace SGPP {
namespace pde {

HeatEquationSolver::HeatEquationSolver() : ParabolicPDESolver() {
  this->bGridConstructed = false;
  this->myScreen = NULL;
}

HeatEquationSolver::~HeatEquationSolver() {
  if (this->myScreen != NULL) {
    delete this->myScreen;
  }
}

void HeatEquationSolver::constructGrid(base::BoundingBox& BoundingBox, int level) {
  this->dim = BoundingBox.getDimensions();
  this->levels = level;

  this->myGrid = new base::LinearBoundaryGrid(BoundingBox);

  base::GridGenerator* myGenerator = this->myGrid->createGridGenerator();
  myGenerator->regular(levels);
  delete myGenerator;

  this->myBoundingBox = this->myGrid->getBoundingBox();
  this->myGridStorage = this->myGrid->getStorage();

  this->bGridConstructed = true;
}

void HeatEquationSolver::setHeatCoefficient(float_t a) { this->a = a; }

void HeatEquationSolver::solveExplicitEuler(size_t numTimesteps, float_t timestepsize,
                                            size_t maxCGIterations, float_t epsilonCG,
                                            base::DataVector& alpha, bool verbose,
                                            bool generateAnimation, size_t numEvalsAnimation) {
  if (this->bGridConstructed) {
    this->myScreen->writeStartSolve("Multidimensional Heat Equation Solver");
    float_t dNeededTime;
    solver::Euler* myEuler = new solver::Euler(
        "ExEul", numTimesteps, timestepsize, generateAnimation, numEvalsAnimation, this->myScreen);
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

    if (this->myScreen != NULL) {
      std::cout << "Time to solve: " << dNeededTime << " seconds" << std::endl;
      this->myScreen->writeEmptyLines(2);
    }

    delete myStopwatch;
    delete myCG;
    delete myEuler;
  } else {
    throw new base::application_exception(
        "HeatEquationSolver::solveExplicitEuler : A grid wasn't constructed before!");
  }
}

void HeatEquationSolver::solveImplicitEuler(size_t numTimesteps, float_t timestepsize,
                                            size_t maxCGIterations, float_t epsilonCG,
                                            base::DataVector& alpha, bool verbose,
                                            bool generateAnimation, size_t numEvalsAnimation) {
  if (this->bGridConstructed) {
    this->myScreen->writeStartSolve("Multidimensional Heat Equation Solver");
    float_t dNeededTime;
    solver::Euler* myEuler = new solver::Euler(
        "ImEul", numTimesteps, timestepsize, generateAnimation, numEvalsAnimation, this->myScreen);
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

    if (this->myScreen != NULL) {
      std::cout << "Time to solve: " << dNeededTime << " seconds" << std::endl;
      this->myScreen->writeEmptyLines(2);
    }

    delete myStopwatch;
    delete myHESolver;
    delete myCG;
    delete myEuler;
  } else {
    throw new base::application_exception(
        "HeatEquationSolver::solveImplicitEuler : A grid wasn't constructed before!");
  }
}

void HeatEquationSolver::solveCrankNicolson(size_t numTimesteps, float_t timestepsize,
                                            size_t maxCGIterations, float_t epsilonCG,
                                            base::DataVector& alpha, size_t NumImEul) {
  if (this->bGridConstructed) {
    this->myScreen->writeStartSolve("Multidimensional Heat Equation Solver");
    float_t dNeededTime;
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
        new solver::Euler("ImEul", numIESteps, timestepsize, false, 0, this->myScreen);
    solver::CrankNicolson* myCN = new solver::CrankNicolson(numCNSteps, timestepsize);

    myStopwatch->start();

    if (numIESteps > 0) {
      myEuler->solve(*myCG, *myHESolver, false);
    }

    myCN->solve(*myCG, *myHESolver, false);
    dNeededTime = myStopwatch->stop();

    if (this->myScreen != NULL) {
      std::cout << "Time to solve: " << dNeededTime << " seconds" << std::endl;
      this->myScreen->writeEmptyLines(2);
    }

    delete myStopwatch;
    delete myHESolver;
    delete myCG;
    delete myCN;
    delete myEuler;
  } else {
    throw new base::application_exception(
        "HeatEquationSolver::solveCrankNicolson : A grid wasn't constructed before!");
  }
}

void HeatEquationSolver::initGridWithSmoothHeat(base::DataVector& alpha, float_t mu, float_t sigma,
                                                float_t factor) {
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
        "HeatEquationSolver::initGridWithSmoothHeat : A grid wasn't constructed before!");
  }
}

void HeatEquationSolver::initScreen() {
  this->myScreen = new base::ScreenOutput();
  this->myScreen->writeTitle("SGpp - Heat Equation Solver, 1.0.0",
                             "Alexander Heinecke, (C) 2009-2011");
}

void HeatEquationSolver::storeInnerRHS(base::DataVector& alpha, std::string tFilename,
                                       float_t timestepsize) {
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
    throw new base::application_exception(
        "HeatEquationSolver::storeInnerMatrix : A grid wasn't constructed before!");
  }
}

void HeatEquationSolver::storeInnerSolution(base::DataVector& alpha, size_t numTimesteps,
                                            float_t timestepsize, size_t maxCGIterations,
                                            float_t epsilonCG, std::string tFilename) {
  if (this->bGridConstructed) {
    solver::Euler* myEuler =
        new solver::Euler("ImEul", numTimesteps, timestepsize, false, 0, this->myScreen);
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
    throw new base::application_exception(
        "HeatEquationSolver::solveImplicitEuler : A grid wasn't constructed before!");
  }
}
}  // namespace pde
}  // namespace SGPP
