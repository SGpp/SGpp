// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/globaldef.hpp>

#include <sgpp/optimization/sle/solver/Auto.hpp>
#include <sgpp/optimization/sle/solver/Armadillo.hpp>
#include <sgpp/optimization/sle/solver/BiCGStab.hpp>
#include <sgpp/optimization/sle/solver/Eigen.hpp>
#include <sgpp/optimization/sle/solver/GaussianElimination.hpp>
#include <sgpp/optimization/sle/solver/Gmmpp.hpp>
#include <sgpp/optimization/sle/solver/UMFPACK.hpp>
#include <sgpp/optimization/tools/Printer.hpp>

#include <cstddef>
#include <algorithm>
#include <map>

namespace SGPP {
  namespace optimization {
    namespace sle_solver {

      const float_t Auto::MAX_NNZ_RATIO_FOR_SPARSE = 0.1;
      const float_t Auto::MAX_NNZ_RATIO_FOR_GMMPP = 0.05;
      const float_t Auto::ESTIMATE_NNZ_ROWS_SAMPLE_SIZE = 0.05;

      /**
       * Add a solver to the vector of solvers, if the solver is supported
       * (e.g. SG++ configured and compiled for use with the solve) and
       * it's not already in the vector.
       *
       * @param solver    linear solver
       * @param solvers   vector of solvers
       * @param supports  map indicating which solvers are supported
       */
      void addSLESolver(SLESolver* solver, std::vector<SLESolver*>& solvers,
                        const std::map<SLESolver*, bool>& supports) {
        // add solver if it's supported and not already in the vector
        if ((supports.at(solver)) &&
            (std::find(solvers.begin(), solvers.end(), solver) ==
                solvers.end())) {
          solvers.push_back(solver);
        }
      }

      bool Auto::solve(SLE& system, const std::vector<float_t>& b,
                       std::vector<float_t>& x) const {
        std::vector<std::vector<float_t>> B;
        std::vector<std::vector<float_t>> X;
        B.push_back(b);
        X.push_back(x);

        if (solve(system, B, X)) {
          x = X[0];
          return true;
        } else {
          return false;
        }
      }

      bool Auto::solve(SLE& system, const std::vector<std::vector<float_t>>& B,
                       std::vector<std::vector<float_t>>& X) const {
        printer.printStatusBegin(
            "Solving linear system (automatic method)...");

        Armadillo solverArmadillo;
        Eigen solverEigen;
        UMFPACK solverUmfpack;
        Gmmpp solverGmmpp;
        BiCGStab solverBicgstab;
        GaussianElimination solverGaussianElimination;

        std::map<SLESolver*, bool> supports;

        // by default, only BiCGStab is supported
        supports[&solverArmadillo] = false;
        supports[&solverEigen] = false;
        supports[&solverUmfpack] = false;
        supports[&solverGmmpp] = false;
        supports[&solverBicgstab] = true;
        supports[&solverGaussianElimination] = true;

#ifdef USEARMADILLO
        supports[&solverArmadillo] = true;
#endif /* USEARMADILLO */

#ifdef USEEIGEN
        supports[&solverEigen] = true;
#endif /* USEEIGEN */

#ifdef USEUMFPACK
        supports[&solverUmfpack] = true;
#endif /* USEUMFPACK */

#ifdef USEGMMPP
        supports[&solverGmmpp] = true;
#endif /* USEGMMPP */

        // solvers to be used, the solver which should be tried first
        // should be the first element
        std::vector<SLESolver*> solvers;
        const size_t n = system.getDimension();

        if (supports[&solverUmfpack] || supports[&solverGmmpp]) {
          // if at least one of the sparse solvers is supported
          // ==> estimate sparsity ratio of matrix by considering
          // every inc-th row
          size_t nrows = 0;
          size_t nnz = 0;
          size_t inc = static_cast<size_t>(
                         ESTIMATE_NNZ_ROWS_SAMPLE_SIZE *
                         static_cast<float_t>(n)) + 1;

          printer.printStatusUpdate("estimating sparsity pattern");

          for (size_t i = 0; i < n; i += inc) {
            nrows++;

            for (size_t j = 0; j < n; j++) {
              if (system.isMatrixEntryNonZero(i, j)) {
                nnz++;
              }
            }
          }

          // calculate estimate ratio nonzero entries
          float_t nnzRatio = static_cast<float_t>(nnz) /
                             (static_cast<float_t>(nrows) *
                                 static_cast<float_t>(n));

          // print ratio
          {
            char str[10];
            snprintf(str, 10, "%.1f%%", nnzRatio * 100.0);
            printer.printStatusUpdate("estimated nnz ratio: " +
                                      std::string(str));
            printer.printStatusNewLine();
          }

          if ((n > MAX_DIM_FOR_FULL) ||
              (nnzRatio <= MAX_NNZ_RATIO_FOR_GMMPP)) {
            if (B.size() == 1) {
              // prefer Gmm++ over UMFPACK
              addSLESolver(&solverGmmpp, solvers, supports);
              addSLESolver(&solverUmfpack, solvers, supports);
            } else {
              // prefer UMFPACK over Gmm++
              // (UMFPACK can solve multiple systems simultaneously)
              addSLESolver(&solverUmfpack, solvers, supports);
              addSLESolver(&solverGmmpp, solvers, supports);
            }
          } else if (nnzRatio <= MAX_NNZ_RATIO_FOR_SPARSE) {
            // prefer UMFPACK over Gmm++
            addSLESolver(&solverUmfpack, solvers, supports);
            addSLESolver(&solverGmmpp, solvers, supports);
          }
        }

        // add all remaining solvers (prefer Armadillo over Eigen)
        addSLESolver(&solverArmadillo, solvers, supports);
        addSLESolver(&solverEigen, solvers, supports);
        addSLESolver(&solverGmmpp, solvers, supports);
        addSLESolver(&solverUmfpack, solvers, supports);
        addSLESolver(&solverBicgstab, solvers, supports);
        addSLESolver(&solverGaussianElimination, solvers, supports);

        for (size_t i = 0; i < solvers.size(); i++) {
          // try solver
          bool result = solvers[i]->solve(system, B, X);

          if (result) {
            printer.printStatusEnd();
            return true;
          } else if ((solvers[i] == &solverGmmpp) &&
              (n > MAX_DIM_FOR_FULL)) {
            // don't use full solvers and return approximative solution
            printer.printStatusEnd(
                "warning: using non-converged solution of iterative "
                "solver, residual can be large "
                "(matrix too large to try other solvers)");
            return true;
          }
        }

        printer.printStatusEnd("error: could not solve linear system!");
        return false;
      }

    }
  }
}
