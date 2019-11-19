// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/tools/Printer.hpp>
#include <sgpp/base/tools/sle/solver/GaussianElimination.hpp>
#include <sgpp/globaldef.hpp>

#include <cmath>
#include <numeric>

namespace sgpp {
namespace base {
namespace sle_solver {

GaussianElimination::~GaussianElimination() {}

bool GaussianElimination::solve(SLE& system, DataVector& b, DataVector& x) const {
  Printer::getInstance().printStatusBegin("Solving linear system (Gaussian elimination)...");

  // size of the system
  const size_t n = b.getSize();
  // working matrix
  DataMatrix W(n, n + 1);

  // set W := (A, b) at the beginning
  for (size_t i = 0; i < n; i++) {
    for (size_t j = 0; j < n; j++) {
      W(i, j) = system.getMatrixEntry(i, j);
    }

    W(i, n) = b[i];
  }

  // at the beginning of the l-th iteration, W should be of the form
  // +-----------------------------------------+
  // |  (1     * - *  |  *)   \                |
  // |  (  \   | . |  |  *)    | l rows        |
  // |  (    1 * - *  |  *)   /                |
  // |  (0 - 0 * - *  |  *)   \                |
  // |  (|   | | - |  |  *)    | (n-l) rows    |
  // |  (0 - 0 * - *  |  *)   /                |
  // |                                         |
  // |  \----/ \---/    \-/                    |
  // |    l    (n-l)     1    column(s)        |
  // +-----------------------------------------+
  for (size_t l = 0; l < n; l++) {
    Printer::getInstance().printStatusUpdate("k = " + std::to_string(l));

    // search for pivot entry = maximum of the absolute values
    // of the entries w_{l,l}, ..., w_{n,l}
    double maxEntry = 0.0;
    // row index of maximal absolute value
    size_t i = l;

    for (size_t j = l; j < n; j++) {
      double entry = std::abs(W(j, l));

      if (entry > maxEntry) {
        maxEntry = entry;
        i = j;
      }
    }

    // all entries are zero ==> matrices W and A are rank deficient
    if (maxEntry == 0.0) {
      Printer::getInstance().printStatusEnd("error: Could not solve linear system!");
      return false;
    }

    // swap rows l and i
    for (size_t k = l; k <= n; k++) {
      const double entry = W(l, k);
      W(l, k) = W(i, k);
      W(i, k) = entry;
    }

    // divide l-th row by w_{l,l}
    {
      const double wll = W(l, l);

      for (size_t k = l; k <= n; k++) {
        W(l, k) /= wll;
      }
    }

    // subtract w_{j,l} times l-th row from all rows j != l
    for (size_t j = 0; j < n; j++) {
      if (j != l) {
        const double wjl = W(j, l);

        for (size_t k = l; k <= n; k++) {
          W(j, k) -= wjl * W(l, k);
        }
      }
    }
  }

  // x is the last column (right side of the augmented matrix)
  x.resize(n);
  W.getColumn(n, x);

  Printer::getInstance().printStatusUpdate("k = " + std::to_string(n));
  Printer::getInstance().printStatusEnd();
  return true;
}
}  // namespace sle_solver
}  // namespace base
}  // namespace sgpp
