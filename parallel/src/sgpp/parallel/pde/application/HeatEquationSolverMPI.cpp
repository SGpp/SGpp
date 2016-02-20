// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/parallel/tools/MPI/SGppMPITools.hpp>
#include <sgpp/parallel/solver/sle/ConjugateGradientsMPI.hpp>
#include <sgpp/parallel/pde/application/HeatEquationSolverMPI.hpp>
#include <sgpp/parallel/pde/algorithm/HeatEquationParabolicPDESolverSystemParallelMPI.hpp>
#include <sgpp/parallel/pde/algorithm/HeatEquationParabolicPDESolverSystemVectorizedMPI.hpp>

#include <sgpp/solver/ode/Euler.hpp>
#include <sgpp/solver/ode/CrankNicolson.hpp>
#include <sgpp/solver/sle/ConjugateGradients.hpp>

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/exception/application_exception.hpp>
#include <sgpp/base/operation/BaseOpFactory.hpp>

#include <sstream>
#include <cstdlib>
#include <cstring>

#include <sgpp/globaldef.hpp>


namespace SGPP {
namespace parallel {

HeatEquationSolverMPI::HeatEquationSolverMPI() :
  SGPP::pde::ParabolicPDESolver() {
  this->bGridConstructed = false;
  this->myScreen = NULL;
}

HeatEquationSolverMPI::~HeatEquationSolverMPI() {
  if (this->myScreen != NULL) {
    delete this->myScreen;
  }
}

void HeatEquationSolverMPI::constructGrid(SGPP::base::BoundingBox& BoundingBox,
    int level) {
  this->dim = BoundingBox.getDimensions();
  this->levels = level;

  this->myGrid = new SGPP::base::LinearBoundaryGrid(BoundingBox);

  std::unique_ptr<SGPP::base::GridGenerator> myGenerator = this->myGrid->createGridGenerator();
  myGenerator->regular(this->levels);
  delete myGenerator;

  this->myBoundingBox = this->myGrid->getBoundingBox();
  this->myGridStorage = this->myGrid->getStorage();

  this->bGridConstructed = true;
}

void HeatEquationSolverMPI::setHeatCoefficient(double a) {
  this->a = a;
}

void HeatEquationSolverMPI::solveExplicitEuler(size_t numTimesteps,
    double timestepsize, size_t maxCGIterations, double epsilonCG,
    SGPP::base::DataVector& alpha, bool verbose, bool generateAnimation,
    size_t numEvalsAnimation) {
  throw new SGPP::base::application_exception("HeatEquationSolver::solveExplicitEuler : Explicit Euler is not supported!");
}

void HeatEquationSolverMPI::solveImplicitEuler(size_t numTimesteps,
    double timestepsize, size_t maxCGIterations, double epsilonCG,
    SGPP::base::DataVector& alpha, bool verbose, bool generateAnimation,
    size_t numEvalsAnimation) {
  if (this->bGridConstructed) {
    if (this->myScreen != NULL) {
      this->myScreen->writeStartSolve("Multidimensional Heat Equation Solver");
    }

    SGPP::solver::SLESolver* myCG;
    SGPP::solver::OperationParabolicPDESolverSystem* myHESolver;

    double dNeededTime;
    SGPP::solver::Euler* myEuler = new SGPP::solver::Euler("ImEul", numTimesteps,
        timestepsize, generateAnimation, numEvalsAnimation, this->myScreen);

    // read env variable, which solver type should be selected
    char* alg_selector = getenv("SGPP_PDE_SOLVER_ALG");

    if (alg_selector != NULL) {
      if (! strcmp(alg_selector, "X86SIMD")) {
        myCG = new solver::ConjugateGradients(maxCGIterations, epsilonCG);
        myHESolver = new HeatEquationParabolicPDESolverSystemVectorizedMPI(
          *this->myGrid, alpha, this->a, timestepsize, "ImEul");
      } else if (! strcmp(alg_selector, "OCL")) {
        myCG = new solver::ConjugateGradients(maxCGIterations, epsilonCG);
        myHESolver = new HeatEquationParabolicPDESolverSystemVectorizedMPI(
          *this->myGrid, alpha, this->a, timestepsize, "ImEul");
      } else {
        throw new base::application_exception("HeatEquationSolverMPI::solveImplicitEuler : You have selected an unsupport vectorization method!");
      }
    } else {
      myCG = new ConjugateGradientsMPI(maxCGIterations, epsilonCG);
      myHESolver = new HeatEquationParabolicPDESolverSystemParallelMPI(*this->myGrid,
          alpha, this->a, timestepsize, "ImEul");
    }

    SGPP::base::SGppStopwatch* myStopwatch = new SGPP::base::SGppStopwatch();

    MPI_Barrier(MPI_COMM_WORLD);
    myStopwatch->start();
    myEuler->solve(*myCG, *myHESolver, verbose);
    MPI_Barrier(MPI_COMM_WORLD);
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
    throw new SGPP::base::application_exception("HeatEquationSolver::solveImplicitEuler : A grid wasn't constructed before!");
  }
}

void HeatEquationSolverMPI::solveCrankNicolson(size_t numTimesteps,
    double timestepsize, size_t maxCGIterations, double epsilonCG,
    SGPP::base::DataVector& alpha, size_t NumImEul) {
  if (this->bGridConstructed) {
    if (this->myScreen != NULL) {
      this->myScreen->writeStartSolve("Multidimensional Heat Equation Solver");
    }

    SGPP::solver::SLESolver* myCG;
    SGPP::solver::OperationParabolicPDESolverSystem* myHESolver;
    double dNeededTime;

    char* alg_selector = getenv("SGPP_PDE_SOLVER_ALG");

    if (alg_selector != NULL) {
      if (! strcmp(alg_selector, "X86SIMD")) {
        myCG = new solver::ConjugateGradients(maxCGIterations, epsilonCG);
        myHESolver = new HeatEquationParabolicPDESolverSystemVectorizedMPI(
          *this->myGrid, alpha, this->a, timestepsize, "CrNic");
      } else if (! strcmp(alg_selector, "OCL")) {
        myCG = new solver::ConjugateGradients(maxCGIterations, epsilonCG);
        myHESolver = new HeatEquationParabolicPDESolverSystemVectorizedMPI(
          *this->myGrid, alpha, this->a, timestepsize, "CrNic");
      } else {
        throw new base::application_exception("HeatEquationSolverMPI::solveCrankNicolson : You have selected an unsupport vectorization method!");
      }
    } else {
      myCG = new ConjugateGradientsMPI(maxCGIterations, epsilonCG);
      myHESolver = new HeatEquationParabolicPDESolverSystemParallelMPI(*this->myGrid,
          alpha, this->a, timestepsize, "CrNic");
    }

    SGPP::base::SGppStopwatch* myStopwatch = new SGPP::base::SGppStopwatch();

    size_t numCNSteps;
    size_t numIESteps;

    numCNSteps = numTimesteps;

    if (numTimesteps > NumImEul) {
      numCNSteps = numTimesteps - NumImEul;
    }

    numIESteps = NumImEul;

    SGPP::solver::Euler* myEuler = new SGPP::solver::Euler("ImEul", numIESteps,
        timestepsize, false, 0, this->myScreen);
    SGPP::solver::CrankNicolson* myCN = new SGPP::solver::CrankNicolson(numCNSteps,
        timestepsize);

    MPI_Barrier(MPI_COMM_WORLD);
    myStopwatch->start();

    if (numIESteps > 0) {
      myHESolver->setODESolver("ImEul");
      myEuler->solve(*myCG, *myHESolver, false);
    }

    myHESolver->setODESolver("CrNic");
    myCN->solve(*myCG, *myHESolver, false);
    MPI_Barrier(MPI_COMM_WORLD);
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
    throw new SGPP::base::application_exception("HeatEquationSolver::solveCrankNicolson : A grid wasn't constructed before!");
  }
}

void HeatEquationSolverMPI::initGridWithSmoothHeat(SGPP::base::DataVector&
    alpha, double mu, double sigma, double factor) {
  if (this->bGridConstructed) {
    double tmp;
    double* dblFuncValues = new double[this->dim];

    for (size_t i = 0; i < this->myGrid->getStorage()->size(); i++) {
      std::string coords = this->myGridStorage->get(i)->getCoordsStringBB(
                             *this->myBoundingBox);
      std::stringstream coordsStream(coords);

      for (size_t j = 0; j < this->dim; j++) {
        coordsStream >> tmp;

        dblFuncValues[j] = tmp;
      }

      tmp = 1.0;

      for (size_t j = 0; j < this->dim; j++) {
        tmp *=  factor * factor * ((1.0 / (sigma * 2.0 * 3.145)) * exp((-0.5) * ((
                                     dblFuncValues[j] - mu) / sigma) * ((dblFuncValues[j] - mu) / sigma)));
      }

      alpha[i] = tmp;
    }

    delete[] dblFuncValues;

    SGPP::base::OperationHierarchisation* myHierarchisation =
      SGPP::op_factory::createOperationHierarchisation(*this->myGrid);
    myHierarchisation->doHierarchisation(alpha);
    delete myHierarchisation;
  } else {
    throw new SGPP::base::application_exception("HeatEquationSolver::initGridWithSmoothHeat : A grid wasn't constructed before!");
  }

}

void HeatEquationSolverMPI::initScreen() {
  this->myScreen = new SGPP::base::ScreenOutput();
  this->myScreen->writeTitle("SGpp - Heat Equation Solver, 1.0.0",
                             "Alexander Heinecke, (C) 2009-2011");
}

}
}
