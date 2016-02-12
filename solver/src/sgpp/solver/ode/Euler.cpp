// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/grid/common/DirichletUpdateVector.hpp>
#include <sgpp/solver/ode/Euler.hpp>
#include <sgpp/base/operation/hash/OperationEval.hpp>
#include <sgpp/base/tools/GridPrinter.hpp>
#include <sgpp/base/exception/solver_exception.hpp>

#include <sgpp/globaldef.hpp>

#include <iostream>
#include <string>
#include <sstream>

namespace SGPP {
namespace solver {

Euler::Euler(std::string Mode, size_t imax, float_t timestepSize, bool generateAnimation,
             size_t numEvalsAnimation, SGPP::base::ScreenOutput* screen)
    : ODESolver(imax, timestepSize),
      bAnimation(generateAnimation),
      evalsAnimation(numEvalsAnimation),
      ExMode(Mode),
      myScreen(screen) {
  this->residuum = 0.0;

  if (Mode != "ExEul" && Mode != "ImEul") {
    throw new SGPP::base::solver_exception("Euler::Euler : An unknown Euler-Mode was specified!");
  }
}

Euler::~Euler() {}

void Euler::solve(SLESolver& LinearSystemSolver,
                  SGPP::solver::OperationParabolicPDESolverSystem& System, bool bIdentifyLastStep,
                  bool verbose) {
  size_t allIter = 0;
  SGPP::base::DataVector* rhs;

  // Do some animation creation exception handling
  size_t animationStep = this->nMaxIterations / 1500;

  if (animationStep == 0) {
    animationStep = 1;
  }

  // Create pictures of the animation, if specified
  if ((this->bAnimation == true)) {
    // Build filename
    std::string tFilename = "00000000000000000000000000000000";
    std::stringstream number;
    number << 0;
    tFilename.append(number.str());
    tFilename = tFilename.substr(tFilename.length() - 14, 14);
    tFilename.append(".gnuplot");

    // Print grid to file
    SGPP::base::GridPrinter myPrinter(*System.getGrid());
    myPrinter.printSparseGrid(*System.getGridCoefficients(), tFilename, false);
  }

  for (size_t i = 0; i < this->nMaxIterations; i++) {
    // generate right hand side
    rhs = System.generateRHS();

    // solve the system of the current timestep
    LinearSystemSolver.solve(System, *System.getGridCoefficientsForCG(), *rhs, true, false, -1.0);

    allIter += LinearSystemSolver.getNumberIterations();

    if (verbose == true) {
      if (myScreen == NULL) {
        std::cout << "Final residuum " << LinearSystemSolver.getResiduum() << "; with "
                  << LinearSystemSolver.getNumberIterations()
                  << " Iterations (Total Iter.: " << allIter << ")" << std::endl;
      }
    }

    if (myScreen != NULL) {
      std::stringstream soutput;
      soutput << "Final residuum " << LinearSystemSolver.getResiduum() << "; with "
              << LinearSystemSolver.getNumberIterations() << " Iterations (Total Iter.: " << allIter
              << ")";

      if (i < this->nMaxIterations - 1) {
        myScreen->update((size_t)(((float_t)(i + 1) * 100.0) / ((float_t) this->nMaxIterations)),
                         soutput.str());
      } else {
        myScreen->update(100, soutput.str());
      }
    }

    System.finishTimestep();

    if (bIdentifyLastStep == false) {
      System.coarsenAndRefine(false);
    } else {
      if (i < (this->nMaxIterations - 1)) {
        System.coarsenAndRefine(false);
      } else {
        System.coarsenAndRefine(true);
      }
    }

    // Create pictures of the animation, if specified
    if ((this->bAnimation == true) && (i % animationStep == 0)) {
      // Build filename
      std::string tFilename = "00000000000000000000000000000000";
      std::stringstream number;
      number << (i + 1);
      tFilename.append(number.str());
      tFilename = tFilename.substr(tFilename.length() - 14, 14);
      tFilename.append(".gnuplot");

      // Print grid to file
      SGPP::base::GridPrinter myPrinter(*System.getGrid());
      myPrinter.printSparseGrid(*System.getGridCoefficients(), tFilename, false);
    }
  }

  // write some empty lines to console
  if (myScreen != NULL) {
    myScreen->writeEmptyLines(2);
  }

  this->nIterations = allIter;
}

}  // namespace solver
}  // namespace SGPP
