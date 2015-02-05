// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/globaldef.hpp>

#include <sgpp/optimization/sle/solver/UMFPACK.hpp>
#include <sgpp/optimization/sle/system/Cloneable.hpp>
#include <sgpp/optimization/tools/Printer.hpp>

#ifdef USEUMFPACK
#include <suitesparse/umfpack.h>
#endif

#include <cstddef>
#include <iostream>
#include <algorithm>

namespace SGPP {
  namespace optimization {
    namespace sle {
      namespace solver {

#ifdef USEUMFPACK
        typedef UF_long sslong;

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
                           const std::vector<float_t>& Ax, const std::vector<float_t>& b,
                           std::vector<float_t>& x) {
          const size_t n = b.size();
          x = std::vector<float_t>(n, 0.0);

          sslong result = umfpack_dl_solve(UMFPACK_A, &Ap[0], &Ai[0], &Ax[0], &x[0], &b[0],
                                           numeric, NULL, NULL);
          return (result == UMFPACK_OK);
        }
#endif

        bool UMFPACK::solve(system::System& system, const std::vector<float_t>& b,
                            std::vector<float_t>& x) const {
          // place RHS in its own vector
          std::vector<std::vector<float_t> > B;
          std::vector<std::vector<float_t> > X;
          B.push_back(b);

          // call version for multiple RHSs
          if (solve(system, B, X)) {
            x = X[0];
            return true;
          } else {
            return false;
          }
        }

        bool UMFPACK::solve(system::System& system, const std::vector<std::vector<float_t> >& B,
                            std::vector<std::vector<float_t> >& X) const {
#ifdef USEUMFPACK
          tools::printer.printStatusBegin("Solving linear system (UMFPACK)...");

          const size_t n = system.getDimension();

          size_t nnz = 0;
          std::vector<uint32_t> Ti;
          std::vector<uint32_t> Tj;
          std::vector<float_t> Tx;

          // parallelize only if the system is cloneable
          #pragma omp parallel if (system.isCloneable()) \
          shared(system, Ti, Tj, Tx, nnz, tools::printer) default(none)
          {
            system::System* system2 = &system;
#ifdef _OPENMP
            std::unique_ptr<system::System> clonedSystem;

            if (system.isCloneable() && (omp_get_max_threads() > 1)) {
              clonedSystem = std::unique_ptr<system::System>(
                               dynamic_cast<system::Cloneable&>(system).clone());
              system2 = clonedSystem.get();
            }

#endif

            // get indices and values of nonzero entries
            #pragma omp for ordered schedule(dynamic)

            for (uint32_t i = 0; i < n; i++) {
              for (uint32_t j = 0; j < n; j++) {
                float_t entry = system2->getMatrixEntry(i, j);

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
                  snprintf(str, 10, "%.1f%%",
                  static_cast<float_t>(i) / static_cast<float_t>(n) * 100.0);
                  tools::printer.printStatusUpdate("constructing sparse matrix (" +
                  std::string(str) + ")");
                }
              }
            }
          }

          tools::printer.printStatusUpdate("constructing sparse matrix (100.0%)");
          tools::printer.printStatusNewLine();

          // print ratio of nonzero entries
          {
            char str[10];
            float_t nnz_ratio = static_cast<float_t>(nnz) /
                                (static_cast<float_t>(n) * static_cast<float_t>(n));
            snprintf(str, 10, "%.1f%%", nnz_ratio * 100.0);
            tools::printer.printStatusUpdate("nnz ratio: " + std::string(str));
            tools::printer.printStatusNewLine();
          }

          std::vector<sslong> Ap(n + 1);
          std::vector<sslong> Ai(nnz);
          std::vector<float_t> Ax(nnz);

          sslong result;

          // convert matrix to CCS
          {
            std::vector<sslong> TiArray(nnz);
            std::vector<sslong> TjArray(nnz);

            for (size_t k = 0; k < nnz; k++) {
              TiArray[k] = static_cast<sslong>(Ti[k]);
              TjArray[k] = static_cast<sslong>(Tj[k]);
            }

            tools::printer.printStatusUpdate("step 1: umfpack_dl_triplet_to_col");

            result = umfpack_dl_triplet_to_col(
                       static_cast<sslong>(n), static_cast<sslong>(n), static_cast<sslong>(nnz),
                       &TiArray[0], &TjArray[0], &Tx[0],
                       &Ap[0], &Ai[0], &Ax[0], NULL);

            if (result != UMFPACK_OK) {
              tools::printer.printStatusEnd(
                "error: could not convert to CCS via umfpack_dl_triplet_to_col, error code " +
                toString(result));
              return false;
            }
          }

          void* symbolic, *numeric;

          tools::printer.printStatusNewLine();
          tools::printer.printStatusUpdate("step 2: umfpack_dl_symbolic");

          // call umfpack_dl_symbolic
          result = umfpack_dl_symbolic(static_cast<sslong>(n), static_cast<sslong>(n),
                                       &Ap[0], &Ai[0], &Ax[0], &symbolic, NULL, NULL);

          if (result != UMFPACK_OK) {
            tools::printer.printStatusEnd("error: could solve via umfpack_dl_symbolic, error code " +
                                          toString(result));
            return false;
          }

          tools::printer.printStatusNewLine();
          tools::printer.printStatusUpdate("step 3: umfpack_dl_numeric");

          // call umfpack_dl_numeric
          result = umfpack_dl_numeric(&Ap[0], &Ai[0], &Ax[0], symbolic, &numeric, NULL, NULL);

          if (result != UMFPACK_OK) {
            tools::printer.printStatusEnd("error: could solve via umfpack_dl_numeric, error code " +
                                          toString(result));
            umfpack_dl_free_symbolic(&symbolic);
            return false;
          }

          umfpack_dl_free_symbolic(&symbolic);

          std::vector<float_t> x;
          X.clear();

          // call umfpack_dl_solve for each RHS
          for (size_t i = 0; i < B.size(); i++) {
            const std::vector<float_t>& b = B[i];
            tools::printer.printStatusNewLine();

            if (B.size() == 1) {
              tools::printer.printStatusUpdate("step 4: umfpack_dl_solve");
            } else {
              tools::printer.printStatusUpdate("step 4: umfpack_dl_solve (RHS " + toString(i + 1) +
                                               " of " + toString(B.size()) + ")");
            }

            if (solveInternal(numeric, Ap, Ai, Ax, b, x)) {
              X.push_back(x);
            } else {
              tools::printer.printStatusEnd("error: could solve via umfpack_dl_solve, error code " +
                                            toString(result));
              umfpack_dl_free_numeric(&numeric);
              return false;
            }
          }

          umfpack_dl_free_numeric(&numeric);
          tools::printer.printStatusEnd();

          return true;
#else
          std::cerr << "Error in SGPP::optimization::sle::solver::UMFPACK::solve: "
                    << "SG++ was compiled without UMFPACK support!\n";
          return false;
#endif
        }

      }
    }
  }
}
