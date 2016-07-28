// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/globaldef.hpp>

#include <sgpp/optimization/sle/solver/MATLABSparse.hpp>
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

MATLABSparse::~MATLABSparse() {}

bool MATLABSparse::solve(SLE& system, base::DataVector& b, base::DataVector& x) const {
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

bool MATLABSparse::solve(SLE& system, base::DataMatrix& B, base::DataMatrix& X) const {
#ifdef USE_MATLAB
  Printer::getInstance().printStatusBegin("Solving linear system (sparse MATLAB)...");

  const size_t n = system.getDimension();
  const size_t m = B.getNcols();

  size_t nnz = 0;
  // non-zero entries (length: nnz)
  std::vector<double> Pr;
  // row indices of non-zero entries (length: nnz)
  std::vector<mwIndex> Ir;
  // Jc[j] contains number of non-zero entries in columns 0, ..., j - 1 (length: n + 1)
  std::vector<mwIndex> Jc(n + 1, 0);

  for (size_t j = 0; j < n; j++) {
    for (size_t i = 0; i < n; i++) {
      if (system.isMatrixEntryNonZero(i, j)) {
        const double entry = system.getMatrixEntry(i, j);

        Pr.push_back(entry);
        Ir.push_back(i);
        nnz++;
      }
    }

    {
      // do cumulative sum to get number of non-zero entries in columns 0, ..., j for Jc[j + 1]
      Jc[j + 1] = nnz;
    }

    // status message
    if (j % 100 == 0) {
      char str[10];
      snprintf(str, sizeof(str), "%.1f%%",
               static_cast<double>(j) / static_cast<double>(n) * 100.0);
      Printer::getInstance().printStatusUpdate("constructing matrix (" + std::string(str) + ")");
    }
  }

  mxArray* plhs[1];
  mxArray* prhs[2];

  prhs[0] = mxCreateSparse(n, n, nnz, mxREAL);
  std::copy(Pr.begin(), Pr.end(), mxGetPr(prhs[0]));
  std::copy(Ir.begin(), Ir.end(), mxGetIr(prhs[0]));
  std::copy(Jc.begin(), Jc.end(), mxGetJc(prhs[0]));

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
