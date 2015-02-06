// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/globaldef.hpp>

#include <sgpp/optimization/sle/solver/Gmmpp.hpp>
#include <sgpp/optimization/sle/system/Cloneable.hpp>
#include <sgpp/optimization/tools/Printer.hpp>

#ifdef USEGMMPP
#include <gmm/gmm.h>
#endif

#include <cstddef>
#include <iostream>

namespace SGPP {
  namespace optimization {
    namespace sle {
      namespace solver {

#ifdef USEGMMPP
        /**
         * Gmm++ status printing callback.
         *
         * @param iter  iteration information
         */
        void callback(const gmm::iteration& iter) {
          tools::printer.printStatusUpdate("solving with Gmm++ "
                                           "(k = " + toString(iter.get_iteration()) +
                                           ", residual norm = " + toString(iter.get_res()) + ")");
        }

        /**
         * @param       A   coefficient matrix
         * @param       b   right-hand side
         * @param[out]  x   solution of the linear system
         * @return          whether all went well (false if errors occurred)
         */
        bool solveInternal(const gmm::csr_matrix<float_t>& A,
                           const std::vector<float_t>& b, std::vector<float_t>& x) {
          // allow warnings
          gmm::warning_level::level(1);
          x.assign(b.size(), 0.0);

          // ILU preconditioning
          tools::printer.printStatusUpdate("constructing preconditioner");
          gmm::ilu_precond<gmm::csr_matrix<float_t> > P(A);

          tools::printer.printStatusNewLine();
          tools::printer.printStatusUpdate("solving with Gmm++");

          gmm::iteration iter(1e-6, 0, 100000);
          iter.set_callback(&callback);

          try {
            // call GMRES
            gmm::gmres(A, x, b, P, 50, iter);

            float_t res = iter.get_res();

            if (iter.converged() && (res < 1e3)) {
              // GMRES converged
              tools::printer.printStatusUpdate("solving with Gmm++ "
                                               "(k = " + toString(iter.get_iteration()) +
                                               ", residual norm = " + toString(res) + ")");
              tools::printer.printStatusEnd();
              return true;
            } else {
              // GMRES didn't converge ==> try again without preconditioner
              gmm::identity_matrix P;

              tools::printer.printStatusNewLine();
              tools::printer.printStatusUpdate(
                "solving with preconditioner failed, trying again without one");
              tools::printer.printStatusNewLine();

              // call GMRES again
              gmm::gmres(A, x, b, P, 50, iter);
              res = iter.get_res();

              if (iter.converged() && (res < 1e3)) {
                tools::printer.printStatusUpdate("solving with Gmm++ "
                                                 "(k = " + toString(iter.get_iteration()) +
                                                 ", residual norm = " + toString(res) + ")");
                tools::printer.printStatusEnd();
                return true;
              } else {
                tools::printer.printStatusEnd(
                  "error: could not solve linear system, method didn't converge");
                return false;
              }
            }
          } catch (std::exception& e) {
            tools::printer.printStatusEnd(
              "error: could not solve linear system, what(): " + std::string(e.what()));
            return false;
          }
        }
#endif

        bool Gmmpp::solve(system::System& system, const std::vector<float_t>& b,
                          std::vector<float_t>& x) const {
#ifdef USEGMMPP
          tools::printer.printStatusBegin("Solving linear system (Gmm++)...");

          const size_t n = system.getDimension();
          size_t nnz = 0;
          gmm::csr_matrix<float_t> A2;

          {
            gmm::row_matrix<gmm::rsvector<float_t> > A(n, n);

            // parallelize only if the system is cloneable
            #pragma omp parallel if (system.isCloneable()) \
            shared(system, A, nnz, tools::printer) default(none)
            {
              system::System* system2 = &system;
#ifdef _OPENMP
              std::unique_ptr<system::System> clonedSystem;

              if (system.isCloneable() && (omp_get_max_threads() > 1)) {
                system::Cloneable* system2Cloneable;
                dynamic_cast<system::Cloneable&>(system).clone(system2Cloneable);
                system2 = system2Cloneable;
                clonedSystem = std::unique_ptr<system::System>(system2);
              }

#endif

              // copy system matrix to Gmm++ matrix object
              #pragma omp for ordered schedule(dynamic)

              for (size_t i = 0; i < n; i++) {
                for (size_t j = 0; j < n; j++) {
                  float_t entry = system2->getMatrixEntry(i, j);

                  if (entry != 0) {
                    #pragma omp critical
                    {
                      A(i, j) = entry;
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

            // the Gmm++ manual said to do so... no idea if that's necessary
            gmm::clean(A, 1e-12);
            gmm::copy(A, A2);
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

          bool result = solveInternal(A2, b, x);
          return result;
#else
          std::cerr << "Error in SGPP::optimization::sle::solver::Gmmpp::solve: "
                    << "SG++ was compiled without Gmm++ support!\n";
          return false;
#endif
        }

      }
    }
  }
}
