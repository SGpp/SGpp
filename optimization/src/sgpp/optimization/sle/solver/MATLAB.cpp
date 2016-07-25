// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/globaldef.hpp>

#include <sgpp/optimization/sle/solver/MATLAB.hpp>
#include <sgpp/optimization/sle/system/CloneableSLE.hpp>
#include <sgpp/optimization/tools/Printer.hpp>

#ifdef USE_MATLAB
#include <mex.h>
#endif /* USE_MATLAB */

#include <cstddef>
#include <cstring>
#include <iostream>
#include <algorithm>
#include <string>
#include <vector>

namespace sgpp {
namespace optimization {
namespace sle_solver {

MATLAB::~MATLAB() {}

bool MATLAB::solve(SLE& system, base::DataVector& b, base::DataVector& x) const {
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

bool MATLAB::solve(SLE& system, base::DataMatrix& B, base::DataMatrix& X) const {
#ifdef USE_MATLAB
  Printer::getInstance().printStatusBegin("Solving linear system (MATLAB)...");

  const size_t n = system.getDimension();
  const size_t m = B.getNcols();

  mxArray* plhs[1];
  mxArray* prhs[2];

  prhs[0] = mxCreateDoubleMatrix(n, n, mxREAL);
  double* arrayA = mxGetPr(prhs[0]);
  size_t nnz = 0;

// parallelize only if the system is cloneable
#pragma omp parallel if (system.isCloneable()) shared(system, arrayA, nnz) default(none)
  {
    SLE* system2 = &system;
#ifdef _OPENMP
    std::unique_ptr<CloneableSLE> clonedSLE;

    if (system.isCloneable() && (omp_get_max_threads() > 1)) {
      dynamic_cast<CloneableSLE&>(system).clone(clonedSLE);
      system2 = clonedSLE.get();
    }

#endif /* _OPENMP */

// copy system matrix to Armadillo matrix object
#pragma omp for ordered schedule(dynamic)

    for (size_t j = 0; j < n; j++) {
      for (size_t i = 0; i < n; i++) {
        arrayA[i + j * n] = system2->getMatrixEntry(i, j);

        // count nonzero entries
        // (not necessary, you can also remove that if you like)
        if (arrayA[i + j * n] != 0) {
#pragma omp atomic
          nnz++;
        }
      }

      // status message
      if (j % 100 == 0) {
#pragma omp ordered
        {
          char str[10];
          snprintf(str, sizeof(str), "%.1f%%",
                   static_cast<double>(j) / static_cast<double>(n) * 100.0);
          Printer::getInstance().printStatusUpdate("constructing matrix (" + std::string(str) +
                                                   ")");
        }
      }
    }
  }

  Printer::getInstance().printStatusUpdate("constructing matrix (100.0%)");
  Printer::getInstance().printStatusNewLine();

  // print ratio of nonzero entries
  {
    char str[10];
    const double nnzRatio = static_cast<double>(nnz) / static_cast<double>(n * n);
    snprintf(str, sizeof(str), "%.1f%%", nnzRatio * 100.0);
    Printer::getInstance().printStatusUpdate("nnz ratio: " + std::string(str));
    Printer::getInstance().printStatusNewLine();
  }

  prhs[1] = mxCreateDoubleMatrix(n, m, mxREAL);
  double* arrayB = mxGetPr(prhs[1]);

  Printer::getInstance().printStatusUpdate("constructing B");
  for (size_t j = 0; j < m; j++) {
    for (size_t i = 0; i < n; i++) {
      arrayB[i + j * n] = B(i, j);
    }
  }

  Printer::getInstance().printStatusUpdate("calling MATLAB's mldivide");
  const int result = mexCallMATLAB(1, plhs, 2, prhs, "mldivide");

  mxDestroyArray(prhs[0]);
  mxDestroyArray(prhs[1]);

  if (result != 0) {
    Printer::getInstance().printStatusEnd(
        "error: Could not solve via mldivide, error code " + std::to_string(result) +
        ". Are you sure you call SG++ from a MEX file in MATLAB?");
    return false;
  }

  Printer::getInstance().printStatusUpdate("copying result");
  const double* arrayX = mxGetPr(plhs[0]);
  X.resize(n, m);

  for (size_t j = 0; j < m; j++) {
    for (size_t i = 0; i < n; i++) {
      X(i, j) = arrayX[i + j * n];
    }
  }

  mxDestroyArray(plhs[0]);
  Printer::getInstance().printStatusEnd();

  return true;
#else
  std::cerr << "Error in sle_solver::MATLAB::solve: "
            << "SG++ was compiled without MATLAB support!\n";
  return false;
#endif /* USE_MATLAB */
}
}  // namespace sle_solver
}  // namespace optimization
}  // namespace sgpp
