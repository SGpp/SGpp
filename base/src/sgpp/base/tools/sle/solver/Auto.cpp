// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/tools/Printer.hpp>
#include <sgpp/base/tools/sle/solver/Armadillo.hpp>
#include <sgpp/base/tools/sle/solver/Auto.hpp>
#include <sgpp/base/tools/sle/solver/BiCGStab.hpp>
#include <sgpp/base/tools/sle/solver/Eigen.hpp>
#include <sgpp/base/tools/sle/solver/GaussianElimination.hpp>
#include <sgpp/base/tools/sle/solver/Gmmpp.hpp>
#include <sgpp/base/tools/sle/solver/UMFPACK.hpp>
#include <sgpp/globaldef.hpp>

#include <algorithm>
#include <cstddef>
#include <map>
#include <string>
#include <vector>

namespace sgpp {
namespace base {
namespace sle_solver {

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
      (std::find(solvers.begin(), solvers.end(), solver) == solvers.end())) {
    solvers.push_back(solver);
  }
}

Auto::~Auto() {}

bool Auto::solve(SLE& system, base::DataVector& b, base::DataVector& x) const {
  base::DataMatrix B(b.getPointer(), b.getSize(), 1);
  base::DataMatrix X(B.getNrows(), B.getNcols());

  // call version for multiple RHSs
  if (solve(system, B, X)) {
    x.resize(X.getNrows());
    X.getColumn(0, x);
    return true;
  } else {
    return false;
  }
}

bool Auto::solve(SLE& system, base::DataMatrix& B, base::DataMatrix& X) const {
  Printer::getInstance().printStatusBegin("Solving linear system (automatic method)...");

  Armadillo solverArmadillo;
  Eigen solverEigen;
  UMFPACK solverUMFPACK;
  Gmmpp solverGmmpp;
  BiCGStab solverBiCGStab;
  GaussianElimination solverGaussianElimination;

  std::map<SLESolver*, bool> supports;

  // by default, only BiCGStab and GaussianElimination supported
  supports[&solverArmadillo] = false;
  supports[&solverEigen] = false;
  supports[&solverUMFPACK] = false;
  supports[&solverGmmpp] = false;
  supports[&solverBiCGStab] = true;
  supports[&solverGaussianElimination] = true;

#ifdef USE_ARMADILLO
  supports[&solverArmadillo] = true;
#endif /* USE_ARMADILLO */

#ifdef USE_EIGEN
  supports[&solverEigen] = true;
#endif /* USE_EIGEN */

#ifdef USE_UMFPACK
  supports[&solverUMFPACK] = true;
#endif /* USE_UMFPACK */

#ifdef USE_GMMPP
  supports[&solverGmmpp] = true;
#endif /* USE_GMMPP */

  // solvers to be used, the solver which should be tried first
  // should be the first element
  std::vector<SLESolver*> solvers;
  const size_t n = system.getDimension();

  if (supports[&solverUMFPACK] || supports[&solverGmmpp]) {
    // if at least one of the sparse solvers is supported
    // ==> estimate sparsity ratio of matrix by considering
    // every inc-th row
    size_t nrows = 0;
    size_t nnz = 0;
    size_t inc = static_cast<size_t>(ESTIMATE_NNZ_ROWS_SAMPLE_SIZE * static_cast<double>(n)) + 1;

    Printer::getInstance().printStatusUpdate("estimating sparsity pattern");

    for (size_t i = 0; i < n; i += inc) {
      nrows++;

      for (size_t j = 0; j < n; j++) {
        if (system.isMatrixEntryNonZero(i, j)) {
          nnz++;
        }
      }
    }

    // calculate estimate ratio nonzero entries
    double nnzRatio =
        static_cast<double>(nnz) / (static_cast<double>(nrows) * static_cast<double>(n));

    // print ratio
    {
      char str[10];
      snprintf(str, sizeof(str), "%.1f%%", nnzRatio * 100.0);
      Printer::getInstance().printStatusUpdate("estimated nnz ratio: " + std::string(str));
      Printer::getInstance().printStatusNewLine();
    }

    if (nnzRatio <= MAX_NNZ_RATIO_FOR_SPARSE) {
      // prefer UMFPACK over Gmm++
      addSLESolver(&solverUMFPACK, solvers, supports);
      addSLESolver(&solverGmmpp, solvers, supports);
    }
  }

  // add all remaining solvers (prefer Armadillo over Eigen)
  addSLESolver(&solverArmadillo, solvers, supports);
  addSLESolver(&solverEigen, solvers, supports);
  addSLESolver(&solverUMFPACK, solvers, supports);
  addSLESolver(&solverGmmpp, solvers, supports);

  // fallback solvers
  if (n <= MAX_DIM_FOR_GAUSSIAN) {
    addSLESolver(&solverGaussianElimination, solvers, supports);
  }

  addSLESolver(&solverBiCGStab, solvers, supports);
  addSLESolver(&solverGaussianElimination, solvers, supports);

  for (size_t i = 0; i < solvers.size(); i++) {
    // try solver
    bool result = solvers[i]->solve(system, B, X);

    if (result) {
      Printer::getInstance().printStatusEnd();
      return true;
    } else if ((solvers[i] == &solverGmmpp) && (n > MAX_DIM_FOR_FULL)) {
      // don't use full solvers and return approximative solution
      Printer::getInstance().printStatusEnd(
          "warning: using non-converged solution of iterative "
          "solver, residual can be large "
          "(matrix too large to try other solvers)");
      return true;
    }
  }

  Printer::getInstance().printStatusEnd("error: Could not solve linear system!");
  return false;
}
}  // namespace sle_solver
}  // namespace base
}  // namespace sgpp
