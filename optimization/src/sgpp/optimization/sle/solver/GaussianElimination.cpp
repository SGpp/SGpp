// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/globaldef.hpp>

#include <sgpp/optimization/sle/solver/GaussianElimination.hpp>
#include <sgpp/optimization/tools/Printer.hpp>
#include <sgpp/base/datatypes/DataMatrix.hpp>

#include <cmath>
#include <numeric>

namespace SGPP {
  namespace optimization {
    namespace sle_solver {

      bool GaussianElimination::solve(SLE& system,
                                      base::DataVector& b,
                                      base::DataVector& x) const {
        printer.printStatusBegin(
          "Solving linear system (Gaussian elimination)...");

        // size of the system
        const size_t n = b.getSize();
        // working matrix
        base::DataMatrix W(n, n + 1);

        // set W := (A, b) at the beginning
        for (size_t i = 0; i < n; i++) {
          for (size_t j = 0; j < n; j++) {
            W.set(i, j, system.getMatrixEntry(i, j));
          }

          W.set(i, n, b[i]);
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
          printer.printStatusUpdate("k = " + std::to_string(l));

          // search for pivot entry = maximum of the absolute values
          // of the entries w_{l,l}, ..., w_{n,l}
          float_t maxEntry = 0;
          // row index of maximal absolute value
          size_t i = l;

          for (size_t j = l; j < n; j++) {
            float_t entry = std::abs(W.get(j, l));

            if (entry > maxEntry) {
              maxEntry = entry;
              i = j;
            }
          }

          // all entries are zero ==> matrices W and A are rank deficient
          if (maxEntry == 0) {
            printer.printStatusEnd(
              "error: could not solve linear system!");
            return false;
          }

          // swap rows l and i
          for (size_t k = l; k <= n; k++) {
            const float_t entry = W.get(l, k);
            W.set(l, k, W.get(i, k));
            W.set(i, k, entry);
          }

          // divide l-th row by w_{l,l}
          {
            const float_t wll = W.get(l, l);

            for (size_t k = l; k <= n; k++) {
              W.set(l, k, W.get(l, k) / wll);
            }
          }

          // subtract w_{j,l} times l-th row from all rows j != l
          for (size_t j = 0; j < n; j++) {
            if (j != l) {
              const float_t wjl = W.get(j, l);

              for (size_t k = l; k <= n; k++) {
                W.set(j, k, W.get(j, k) - wjl * W.get(l, k));
              }
            }
          }
        }

        // x is the last column (right side of the augmented matrix)
        x.resize(n);
        W.getColumn(n, x);

        printer.printStatusUpdate("k = " + std::to_string(n));
        printer.printStatusEnd();
        return true;
      }

    }
  }
}
