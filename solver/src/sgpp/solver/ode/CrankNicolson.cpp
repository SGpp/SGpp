// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/grid/common/DirichletUpdateVector.hpp>
#include <sgpp/solver/ode/CrankNicolson.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {
  namespace solver {

    CrankNicolson::CrankNicolson(size_t nTimesteps, float_t timestepSize, SGPP::base::ScreenOutput* screen) : ODESolver(nTimesteps, timestepSize), myScreen(screen) {
      this->residuum = 0.0;
    }

    CrankNicolson::~CrankNicolson() {
    }

    void CrankNicolson::solve(SLESolver& LinearSystemSolver, SGPP::solver::OperationParabolicPDESolverSystem& System, bool bIdentifyLastStep, bool verbose) {
      size_t allIter = 0;
      SGPP::base::DataVector* rhs = NULL;

      for (size_t i = 0; i < this->nMaxIterations; i++) {
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
            myScreen->update((size_t)(((float_t)(i + 1) * 100.0) / ((float_t)this->nMaxIterations)), soutput.str());
          } else {
            myScreen->update(100, soutput.str());
          }
        }

        System.finishTimestep();

        // Do some adjustments on the boundaries if needed, copy values back
        if (bIdentifyLastStep == false) {
          System.coarsenAndRefine(false);
        } else {
          if (i < (this->nMaxIterations - 1)) {
            System.coarsenAndRefine(false);
          } else {
            System.coarsenAndRefine(true);
          }
        }
      }

      // write some empty lines to console
      if (myScreen != NULL) {
        myScreen->writeEmptyLines(2);
      }

      this->nIterations = allIter;
    }

  }
}
