/* ****************************************************************************
* Copyright (C) 2011 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Peter Hoffmann (peter.hoffmann@mytum.de)

#include "base/grid/common/DirichletUpdateVector.hpp"
#include "solver/ode/AdamsBashforth.hpp"
#include "base/operation/OperationEval.hpp"
#include "base/tools/GridPrinter.hpp"
#include "base/exception/solver_exception.hpp"

#include <iostream>
#include <string>
#include <sstream>

namespace sg {
  namespace solver {

    AdamsBashforth::AdamsBashforth(size_t imax, double timestepSize, sg::base::ScreenOutput* screen) : ODESolver(imax, timestepSize), myScreen(screen) {
      this->residuum = 0.0;

    }

    AdamsBashforth::~AdamsBashforth() {
    }

    void AdamsBashforth::solve(SLESolver& LinearSystemSolver, sg::pde::OperationParabolicPDESolverSystem& System, bool bIdentifyLastStep, bool verbose) {
      size_t allIter = 0;
      sg::base::DataVector* rhs;

      for (size_t i = 0; i < this->nMaxIterations; i++) {
        if (i > 0)
          System.setODESolver("AdBas");
        else
          System.setODESolver("ExEul");

        // generate right hand side
        rhs = System.generateRHS();


        // solve the system of the current timestep
        LinearSystemSolver.solve(System, *System.getGridCoefficientsForCG(), *rhs, true, false, -1.0);

        allIter += LinearSystemSolver.getNumberIterations();

        if (verbose == true) {
          if (myScreen == NULL) {
            std::cout << "Final residuum " << LinearSystemSolver.getResiduum() << "; with " << LinearSystemSolver.getNumberIterations() << " Iterations (Total Iter.: " << allIter << ")" << std::endl;
          }
        }

        if (myScreen != NULL) {
          std::stringstream soutput;
          soutput << "Final residuum " << LinearSystemSolver.getResiduum() << "; with " << LinearSystemSolver.getNumberIterations() << " Iterations (Total Iter.: " << allIter << ")";

          if (i < this->nMaxIterations - 1) {
            myScreen->update((size_t)(((double)(i + 1) * 100.0) / ((double)this->nMaxIterations)), soutput.str());
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

        System.saveAlpha();

      }

      // write some empty lines to console
      if (myScreen != NULL) {
        myScreen->writeEmptyLines(2);
      }

      this->nIterations = allIter;
    }

  }
}
