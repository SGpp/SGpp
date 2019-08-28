// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/pde/algorithm/HeatEquationParabolicPDESolverSystem.hpp>
#include <sgpp/pde/application/HeatEquationSolverWithStretching.hpp>
#include <sgpp/pde/algorithm/HeatEquationParabolicPDESolverSystemParallelOMP.hpp>
#include <sgpp/solver/ode/Euler.hpp>
#include <sgpp/solver/ode/CrankNicolson.hpp>
#include <sgpp/solver/sle/ConjugateGradients.hpp>
#include <sgpp/base/exception/application_exception.hpp>
#include <sgpp/base/operation/BaseOpFactory.hpp>
#include <sgpp/base/tools/SGppStopwatch.hpp>

#include <sgpp/globaldef.hpp>

#include <string>

namespace sgpp {
namespace pde {

HeatEquationSolverWithStretching::HeatEquationSolverWithStretching() : ParabolicPDESolver() {
  this->bGridConstructed = false;
  this->myScreen = nullptr;
}

HeatEquationSolverWithStretching::~HeatEquationSolverWithStretching() {
  if (this->myScreen != nullptr) {
    delete this->myScreen;
  }
}

void HeatEquationSolverWithStretching::constructGrid(base::Stretching& stretching, size_t level) {
  this->dim = stretching.getDimension();
  this->levels = static_cast<int>(level);

  this->myGrid = new base::LinearStretchedBoundaryGrid(stretching);

  this->myGrid->getGenerator().regular(this->levels);

  this->myStretching = &this->myGrid->getStretching();
  this->myGridStorage = &this->myGrid->getStorage();

  this->bGridConstructed = true;
}

void HeatEquationSolverWithStretching::constructGrid(base::BoundingBox& BoundingBox, size_t level) {
  std::cout << "I'm not supposed to be here, me is constructGrid\n";
}

void HeatEquationSolverWithStretching::setHeatCoefficient(double a) { this->a = a; }

void HeatEquationSolverWithStretching::solveExplicitEuler(size_t numTimesteps, double timestepsize,
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
        "HeatEquationSolverWithStretching::solveExplicitEuler : A grid wasn't constructed before!");
  }
}

void HeatEquationSolverWithStretching::solveImplicitEuler(size_t numTimesteps, double timestepsize,
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
        "HeatEquationSolverWithStretching::solveImplicitEuler : A grid wasn't constructed before!");
  }
}

void HeatEquationSolverWithStretching::solveCrankNicolson(size_t numTimesteps, double timestepsize,
                                                          size_t maxCGIterations, double epsilonCG,
                                                          base::DataVector& alpha,
                                                          size_t NumImEul) {
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
        "HeatEquationSolverWithStretching::solveCrankNicolson : A grid wasn't constructed before!");
  }
}

void HeatEquationSolverWithStretching::initGridWithSmoothHeat(base::DataVector& alpha, double mu,
                                                              double sigma, double factor) {
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
        "HeatEquationSolverWithStretching::initGridWithSmoothHeat : A grid wasn't constructed "
        "before!");
  }
}

void HeatEquationSolverWithStretching::initScreen() {
  this->myScreen = new base::ScreenOutput();
  this->myScreen->writeTitle("SGpp - Heat Equation Solver With Stretching, 1.0.1",
                             "Alexander Heinecke, Sarpkan Selcuk (C) 2009-2011");
}

void HeatEquationSolverWithStretching::printGrid(base::DataVector& alpha,
                                                 size_t PointesPerDimension,
                                                 std::string tfilename) const {
  base::GridPrinterForStretching myPrinter(*this->myGrid);
  myPrinter.printGrid(alpha, tfilename, PointesPerDimension);
}

void HeatEquationSolverWithStretching::printGridDomain(base::DataVector& alpha,
                                                       size_t PointesPerDimension,
                                                       base::BoundingBox& GridArea,
                                                       std::string tfilename) const {
  throw base::application_exception(
      "HeatEquationSolverWithStretching::printGridDomain : BoundingBox not supported with this "
      "solver, use printGridDomainStretching instead ");
}

void HeatEquationSolverWithStretching::printGridDomainStretching(base::DataVector& alpha,
                                                                 size_t PointesPerDimension,
                                                                 base::Stretching& GridArea,
                                                                 std::string tfilename) const {
  base::GridPrinterForStretching myPrinter(*this->myGrid);
  myPrinter.printGridDomainStretching(alpha, tfilename, GridArea, PointesPerDimension);
}

void HeatEquationSolverWithStretching::printSparseGrid(base::DataVector& alpha,
                                                       std::string tfilename, bool bSurplus) const {
  base::GridPrinterForStretching myPrinter(*this->myGrid);
  myPrinter.printSparseGrid(alpha, tfilename, bSurplus);
}

void HeatEquationSolverWithStretching::printSparseGridExpTransform(base::DataVector& alpha,
                                                                   std::string tfilename,
                                                                   bool bSurplus) const {
  base::GridPrinterForStretching myPrinter(*this->myGrid);
  myPrinter.printSparseGridExpTransform(alpha, tfilename, bSurplus);
}
}  // namespace pde
}  // namespace sgpp
