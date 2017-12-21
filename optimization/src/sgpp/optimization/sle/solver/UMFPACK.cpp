// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/globaldef.hpp>

#include <sgpp/optimization/sle/solver/UMFPACK.hpp>
#include <sgpp/optimization/sle/system/CloneableSLE.hpp>
#include <sgpp/optimization/tools/Printer.hpp>

#ifdef USE_UMFPACK
#include <suitesparse/umfpack.h>
#endif /* USE_UMFPACK */

#include <algorithm>
#include <cstddef>
#include <iostream>
#include <string>
#include <vector>

namespace sgpp {
namespace optimization {
namespace sle_solver {

#ifdef USE_UMFPACK
#if SUITESPARSE_MAIN_VERSION >= 4
typedef SuiteSparse_long sslong;
#else
typedef UF_long sslong;
#endif

/**
 * @param       numeric result of umfpack_dl_numeric()
 * @param       Ap      CCS column pointers
 * @param       Ai      CCS row indices
 * @param       Ax      CCS matrix entries
 * @param       b       right-hand side
 * @param[out]  x       solution to the system
 * @return              whether all went well (false if errors occurred)
 */
bool solveInternal(void* numeric, const std::vector<sslong>& Ap, const std::vector<sslong>& Ai,
                   const std::vector<double>& Ax, base::DataVector& b, base::DataVector& x) {
  const size_t n = b.getSize();

  x.resize(n);
  x.setAll(0.0);

  sslong result = umfpack_dl_solve(UMFPACK_A, &Ap[0], &Ai[0], &Ax[0], x.getPointer(),
                                   b.getPointer(), numeric, NULL, NULL);
  return (result == UMFPACK_OK);
}
#endif /* USE_UMFPACK */

UMFPACK::~UMFPACK() {}

bool UMFPACK::solve(SLE& system, base::DataVector& b, base::DataVector& x) const {
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

bool UMFPACK::solve(SLE& system, base::DataMatrix& B, base::DataMatrix& X) const {
#ifdef USE_UMFPACK
  Printer::getInstance().printStatusBegin("Solving linear system (UMFPACK)...");

  const size_t n = system.getDimension();

  size_t nnz = 0;
  std::vector<uint32_t> Ti;
  std::vector<uint32_t> Tj;
  std::vector<double> Tx;

// parallelize only if the system is cloneable
#pragma omp parallel if (system.isCloneable() && (n > 5000)) shared(system, Ti, Tj, Tx, \
                                                                    nnz) default(none)
  {
    SLE* system2 = &system;
#ifdef _OPENMP
    std::unique_ptr<CloneableSLE> clonedSLE;

    if (system.isCloneable() && (omp_get_max_threads() > 1)) {
      dynamic_cast<CloneableSLE&>(system).clone(clonedSLE);
      system2 = clonedSLE.get();
    }

#endif /* _OPENMP */

// get indices and values of nonzero entries
#pragma omp for ordered schedule(dynamic)

    for (uint32_t i = 0; i < n; i++) {
      for (uint32_t j = 0; j < n; j++) {
        double entry = system2->getMatrixEntry(i, j);

        if (entry != 0) {
#pragma omp critical
          {
            Ti.push_back(i);
            Tj.push_back(j);
            Tx.push_back(entry);
            nnz++;
          }
        }
      }

      // status message
      if (i % 100 == 0) {
#pragma omp ordered

        {
          char str[10];
          snprintf(str, sizeof(str), "%.1f%%",
                   static_cast<double>(i) / static_cast<double>(n) * 100.0);
          Printer::getInstance().printStatusUpdate("constructing sparse matrix (" +
                                                   std::string(str) + ")");
        }
      }
    }
  }

  Printer::getInstance().printStatusUpdate("constructing sparse matrix (100.0%)");
  Printer::getInstance().printStatusNewLine();

  // print ratio of nonzero entries
  {
    char str[10];
    double nnz_ratio = static_cast<double>(nnz) / (static_cast<double>(n) * static_cast<double>(n));
    snprintf(str, sizeof(str), "%.1f%%", nnz_ratio * 100.0);
    Printer::getInstance().printStatusUpdate("nnz ratio: " + std::string(str));
    Printer::getInstance().printStatusNewLine();
  }

  std::vector<sslong> Ap(n + 1, 0);
  std::vector<sslong> Ai(nnz, 0);
  std::vector<double> Ax(nnz, 0.0);

  sslong result;

  // convert matrix to CCS
  {
    std::vector<sslong> TiArray(nnz, 0);
    std::vector<sslong> TjArray(nnz, 0);

    for (size_t k = 0; k < nnz; k++) {
      TiArray[k] = static_cast<sslong>(Ti[k]);
      TjArray[k] = static_cast<sslong>(Tj[k]);
    }

    Printer::getInstance().printStatusUpdate("step 1: umfpack_dl_triplet_to_col");

    result = umfpack_dl_triplet_to_col(static_cast<sslong>(n), static_cast<sslong>(n),
                                       static_cast<sslong>(nnz), &TiArray[0], &TjArray[0], &Tx[0],
                                       &Ap[0], &Ai[0], &Ax[0], NULL);

    if (result != UMFPACK_OK) {
      Printer::getInstance().printStatusEnd(
          "error: Could not convert to CCS via "
          "umfpack_dl_triplet_to_col, error code " +
          std::to_string(result));
      return false;
    }
  }

  void *symbolic, *numeric;

  Printer::getInstance().printStatusNewLine();
  Printer::getInstance().printStatusUpdate("step 2: umfpack_dl_symbolic");

  // call umfpack_dl_symbolic
  result = umfpack_dl_symbolic(static_cast<sslong>(n), static_cast<sslong>(n), &Ap[0], &Ai[0],
                               &Ax[0], &symbolic, NULL, NULL);

  if (result != UMFPACK_OK) {
    Printer::getInstance().printStatusEnd(
        "error: Could solve via umfpack_dl_symbolic, "
        "error code " +
        std::to_string(result));
    return false;
  }

  Printer::getInstance().printStatusNewLine();
  Printer::getInstance().printStatusUpdate("step 3: umfpack_dl_numeric");

  // call umfpack_dl_numeric
  result = umfpack_dl_numeric(&Ap[0], &Ai[0], &Ax[0], symbolic, &numeric, NULL, NULL);

  if (result != UMFPACK_OK) {
    Printer::getInstance().printStatusEnd(
        "error: Could solve via umfpack_dl_numeric, "
        "error code " +
        std::to_string(result));
    umfpack_dl_free_symbolic(&symbolic);
    return false;
  }

  umfpack_dl_free_symbolic(&symbolic);

  base::DataVector x(n);
  base::DataVector b(n);
  X.resize(n, B.getNcols());

  // call umfpack_dl_solve for each RHS
  for (size_t i = 0; i < B.getNcols(); i++) {
    B.getColumn(i, b);
    Printer::getInstance().printStatusNewLine();

    if (B.getNcols() == 1) {
      Printer::getInstance().printStatusUpdate("step 4: umfpack_dl_solve");
    } else {
      Printer::getInstance().printStatusUpdate("step 4: umfpack_dl_solve (RHS " +
                                               std::to_string(i + 1) + " of " +
                                               std::to_string(B.getNcols()) + ")");
    }

    if (solveInternal(numeric, Ap, Ai, Ax, b, x)) {
      X.setColumn(i, x);
    } else {
      Printer::getInstance().printStatusEnd(
          "error: Could solve via umfpack_dl_solve, "
          "error code " +
          std::to_string(result));
      umfpack_dl_free_numeric(&numeric);
      return false;
    }
  }

  umfpack_dl_free_numeric(&numeric);
  Printer::getInstance().printStatusEnd();

  return true;
#else
  std::cerr << "Error in sle_solver::UMFPACK::solve: "
            << "SG++ was compiled without UMFPACK support!\n";
  return false;
#endif /* USE_UMFPACK */
}
}  // namespace sle_solver
}  // namespace optimization
}  // namespace sgpp
