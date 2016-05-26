// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/parallel/tools/MPI/SGppMPITools.hpp>
#include <sgpp/parallel/solver/sle/ConjugateGradientsMPI.hpp>
#include <sgpp/parallel/pde/application/PoissonEquationSolverMPI.hpp>
#include <sgpp/parallel/pde/algorithm/PoissonEquationEllipticPDESolverSystemDirichletParallelMPI.hpp>
#include <sgpp/parallel/pde/algorithm/PoissonEquationEllipticPDESolverSystemDirichletVectorizedMPI.hpp>

#include <sgpp/solver/sle/ConjugateGradients.hpp>

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/exception/application_exception.hpp>
#include <sgpp/base/tools/SGppStopwatch.hpp>
#include <sgpp/base/operation/BaseOpFactory.hpp>
#include <sgpp/globaldef.hpp>

#include <cstdlib>
#include <sstream>
#include <cstring>
#include <string>

namespace sgpp {
namespace parallel {

PoissonEquationSolverMPI::PoissonEquationSolverMPI() : sgpp::pde::EllipticPDESolver() {
  this->bGridConstructed = false;
  this->myScreen = NULL;
}

PoissonEquationSolverMPI::~PoissonEquationSolverMPI() {
  if (this->myScreen != NULL) {
    delete this->myScreen;
  }
}

void PoissonEquationSolverMPI::constructGrid(sgpp::base::BoundingBox& BoundingBox, size_t level) {
  this->dim = BoundingBox.getDimensions();
  this->levels = static_cast<int>(level);

  this->myGrid = new sgpp::base::LinearBoundaryGrid(BoundingBox);

  this->myGrid->getGenerator().regular(this->levels);

  this->myBoundingBox = &this->myGrid->getBoundingBox();
  this->myGridStorage = &this->myGrid->getStorage();

  this->bGridConstructed = true;
}

void PoissonEquationSolverMPI::solvePDE(sgpp::base::DataVector& alpha, sgpp::base::DataVector& rhs,
                                        size_t maxCGIterations, double epsilonCG, bool verbose) {
  double dTimeAlpha = 0.0;
  double dTimeRHS = 0.0;
  double dTimeSolver = 0.0;
  sgpp::solver::SLESolver* myCG;
  sgpp::pde::OperationEllipticPDESolverSystemDirichlet* mySystem;

  sgpp::base::SGppStopwatch* myStopwatch = new sgpp::base::SGppStopwatch();

  // read env variable, which solver type should be selected
  char* alg_selector = getenv("SGPP_PDE_SOLVER_ALG");

  if (alg_selector != NULL) {
    if (!strcmp(alg_selector, "X86SIMD")) {
      myCG = new solver::ConjugateGradients(maxCGIterations, epsilonCG);
      mySystem =
          new PoissonEquationEllipticPDESolverSystemDirichletVectorizedMPI(*(this->myGrid), rhs);
    } else if (!strcmp(alg_selector, "OCL")) {
      myCG = new solver::ConjugateGradients(maxCGIterations, epsilonCG);
      mySystem =
          new PoissonEquationEllipticPDESolverSystemDirichletVectorizedMPI(*(this->myGrid), rhs);
    } else {
      throw base::application_exception(
          "BlackScholesSolverMPI::solveImplicitEuler : You have selected an unsupport "
          "vectorization method!");
    }
  } else {
    myCG = new ConjugateGradientsMPI(maxCGIterations, epsilonCG);
    mySystem = new PoissonEquationEllipticPDESolverSystemDirichletParallelMPI(*(this->myGrid), rhs);
  }

  if (myGlobalMPIComm->getMyRank() == 0) {
    std::cout << "Gridpoints (complete grid): " << mySystem->getNumGridPointsComplete()
              << std::endl;
    std::cout << "Gridpoints (inner grid): " << mySystem->getNumGridPointsInner() << std::endl
              << std::endl
              << std::endl;
  }

  MPI_Barrier(MPI_COMM_WORLD);
  myStopwatch->start();
  sgpp::base::DataVector* alpha_solve = mySystem->getGridCoefficientsForCG();
  MPI_Barrier(MPI_COMM_WORLD);
  dTimeAlpha = myStopwatch->stop();

  // generate RHS
  sgpp::base::DataVector* rhs_solve = NULL;

  if (myGlobalMPIComm->getMyRank() == 0)
    std::cout << "coefficients has been initialized for solving!" << std::endl;

  MPI_Barrier(MPI_COMM_WORLD);
  myStopwatch->start();
  rhs_solve = mySystem->generateRHS();
  MPI_Barrier(MPI_COMM_WORLD);
  dTimeRHS = myStopwatch->stop();

  if (myGlobalMPIComm->getMyRank() == 0)
    std::cout << "right hand side has been initialized for solving!" << std::endl
              << std::endl
              << std::endl;

  MPI_Barrier(MPI_COMM_WORLD);
  myStopwatch->start();
  myCG->solve(*mySystem, *alpha_solve, *rhs_solve, true, verbose, 0.0);
  // Copy result into coefficient vector of the boundary grid
  mySystem->getSolutionBoundGrid(alpha, *alpha_solve);
  MPI_Barrier(MPI_COMM_WORLD);
  dTimeSolver = myStopwatch->stop();

  if (myGlobalMPIComm->getMyRank() == 0) {
    std::cout << std::endl << std::endl;
    std::cout << "Gridpoints (complete grid): " << mySystem->getNumGridPointsComplete()
              << std::endl;
    std::cout << "Gridpoints (inner grid): " << mySystem->getNumGridPointsInner() << std::endl
              << std::endl
              << std::endl;

    std::cout << "Timings for solving Poisson Equation" << std::endl;
    std::cout << "------------------------------------" << std::endl;
    std::cout << "Time for creating CG coeffs: " << dTimeAlpha << std::endl;
    std::cout << "Time for creating RHS: " << dTimeRHS << std::endl;
    std::cout << "Time for solving: " << dTimeSolver << std::endl << std::endl;
    std::cout << "Time: " << dTimeAlpha + dTimeRHS + dTimeSolver << std::endl;
    size_t entries = mySystem->getNumGridPointsInner() * mySystem->getNumGridPointsInner();
    size_t flop = (26 * this->dim + this->dim * this->dim) * entries;
    std::cout << "GFLOPS: "
              << 1.0 * static_cast<double>(myCG->getNumberIterations() * flop) / dTimeSolver /
                     1000000000.0
              << std::endl;
    std::cout << "GB/s: "
              << (static_cast<double>(myCG->getNumberIterations() * entries * sizeof(double))) /
                     dTimeSolver / 1000000000.0
              << std::endl
              << std::endl
              << std::endl;
  }

  delete myCG;
  delete mySystem;  // alpha_solver and rhs_solve are allocated and freed here!!
  delete myStopwatch;
}

void PoissonEquationSolverMPI::initGridWithSmoothHeat(sgpp::base::DataVector& alpha, double mu,
                                                      double sigma, double factor) {
  if (this->bGridConstructed) {
    double tmp;
    double* dblFuncValues = new double[this->dim];

    for (size_t i = 0; i < this->myGrid->getStorage().getSize(); i++) {
      std::string coords = this->myGrid->getStorage().getCoordinates(
          this->myGrid->getStorage().getPoint(i)).toString();
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

    std::unique_ptr<sgpp::base::OperationHierarchisation> myHierarchisation =
        sgpp::op_factory::createOperationHierarchisation(*this->myGrid);
    myHierarchisation->doHierarchisation(alpha);
  } else {
    throw sgpp::base::application_exception(
        "PoissonEquationSolverMPI::initGridWithSmoothHeat : A grid wasn't constructed before!");
  }
}

void PoissonEquationSolverMPI::initGridWithSmoothHeatFullDomain(sgpp::base::DataVector& alpha,
                                                                double mu, double sigma,
                                                                double factor) {
  if (this->bGridConstructed) {
    double tmp;
    double* dblFuncValues = new double[this->dim];

    for (size_t i = 0; i < this->myGrid->getStorage().getSize(); i++) {
      std::string coords = this->myGrid->getStorage().getCoordinates(
          this->myGrid->getStorage().getPoint(i)).toString();
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

    std::unique_ptr<sgpp::base::OperationHierarchisation> myHierarchisation =
        sgpp::op_factory::createOperationHierarchisation(*this->myGrid);
    myHierarchisation->doHierarchisation(alpha);
  } else {
    throw sgpp::base::application_exception(
        "HeatEquationSolverMPI::initGridWithSmoothHeatFullDomain : A grid wasn't constructed "
        "before!");
  }
}

void PoissonEquationSolverMPI::initGridWithExpHeat(sgpp::base::DataVector& alpha, double factor) {
  if (this->bGridConstructed) {
    double tmp;
    double* dblFuncValues = new double[this->dim];
    double* rightBound = new double[this->dim];

    sgpp::base::BoundingBox* tmpBB = &this->myGrid->getBoundingBox();

    for (size_t j = 0; j < this->dim; j++) {
      rightBound[j] = (tmpBB->getBoundary(j)).rightBoundary;
    }

    for (size_t i = 0; i < this->myGrid->getStorage().getSize(); i++) {
      std::string coords = this->myGrid->getStorage().getCoordinates(
          this->myGrid->getStorage().getPoint(i)).toString();
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

    std::unique_ptr<sgpp::base::OperationHierarchisation> myHierarchisation =
        sgpp::op_factory::createOperationHierarchisation(*this->myGrid);
    myHierarchisation->doHierarchisation(alpha);
  } else {
    throw sgpp::base::application_exception(
        "PoissonEquationSolverMPI::initGridWithExpHeat : A grid wasn't constructed before!");
  }
}

void PoissonEquationSolverMPI::initGridWithExpHeatFullDomain(sgpp::base::DataVector& alpha,
                                                             double factor) {
  if (this->bGridConstructed) {
    double tmp;
    double* dblFuncValues = new double[this->dim];
    double* rightBound = new double[this->dim];

    sgpp::base::BoundingBox* tmpBB = &this->myGrid->getBoundingBox();

    for (size_t j = 0; j < this->dim; j++) {
      rightBound[j] = (tmpBB->getBoundary(j)).rightBoundary;
    }

    for (size_t i = 0; i < this->myGrid->getStorage().getSize(); i++) {
      std::string coords = this->myGrid->getStorage().getCoordinates(
          this->myGrid->getStorage().getPoint(i)).toString();
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

    std::unique_ptr<sgpp::base::OperationHierarchisation> myHierarchisation =
        sgpp::op_factory::createOperationHierarchisation(*this->myGrid);
    myHierarchisation->doHierarchisation(alpha);
  } else {
    throw sgpp::base::application_exception(
        "PoissonEquationSolverMPI::initGridWithExpHeat : A grid wasn't constructed before!");
  }
}

void PoissonEquationSolverMPI::initScreen() {
  this->myScreen = new sgpp::base::ScreenOutput();
  this->myScreen->writeTitle("SGpp - Poisson Equation Solver, 1.0.0",
                             "Alexander Heinecke, (C) 2009-2011");
}
}  // namespace parallel
}  // namespace sgpp
