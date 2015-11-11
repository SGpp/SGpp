// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/globaldef.hpp>

#include <sgpp/optimization/sle/solver/Gmmpp.hpp>
#include <sgpp/optimization/sle/system/CloneableSLE.hpp>
#include <sgpp/optimization/tools/Printer.hpp>

#ifdef USE_GMMPP
#include <gmm/gmm.h>
#endif /* USE_GMMPP */

#include <cstddef>
#include <iostream>

namespace SGPP {
  namespace optimization {
    namespace sle_solver {

#ifdef USE_GMMPP
      /**
       * Gmm++ status printing callback.
       *
       * @param iter  iteration information
       */
      void callback(const gmm::iteration& iter) {
        Printer::getInstance().printStatusUpdate(
          "solving with Gmm++ (k = " +
          std::to_string(iter.get_iteration()) +
          ", residual norm = " + std::to_string(iter.get_res()) + ")");
      }

      /**
       * @param       A   coefficient matrix
       * @param       b   right-hand side
       * @param[out]  x   solution of the linear system
       * @return          whether all went well (false if errors occurred)
       */
      bool solveInternal(const gmm::csr_matrix<float_t>& A,
                         base::DataVector& b,
                         base::DataVector& x) {
        // allow warnings
        gmm::warning_level::level(1);
        const size_t n = b.getSize();
        std::vector<float_t> bVec(b.getPointer(), b.getPointer() + n);
        std::vector<float_t> xVec(n, 0.0);

        // ILU preconditioning
        Printer::getInstance().printStatusUpdate("constructing preconditioner");
        gmm::ilu_precond<gmm::csr_matrix<float_t>> P(A);

        Printer::getInstance().printStatusNewLine();
        Printer::getInstance().printStatusUpdate("solving with Gmm++");

        gmm::iteration iter(1e-6, 0, 10000);
        iter.set_callback(&callback);

        try {
          // call GMRES
          gmm::gmres(A, xVec, bVec, P, 50, iter);

          float_t res = iter.get_res();

          if (iter.converged() && (res < 1e3)) {
            // GMRES converged
            x = base::DataVector(xVec);
            Printer::getInstance().printStatusUpdate(
              "solving with Gmm++ (k = " +
              std::to_string(iter.get_iteration()) +
              ", residual norm = " + std::to_string(res) + ")");
            Printer::getInstance().printStatusEnd();
            return true;
          } else {
            // GMRES didn't converge ==> try again without preconditioner
            gmm::identity_matrix P;

            Printer::getInstance().printStatusNewLine();
            Printer::getInstance().printStatusUpdate(
              "solving with preconditioner failed, trying again without one");
            Printer::getInstance().printStatusNewLine();

            // call GMRES again
            gmm::gmres(A, xVec, bVec, P, 50, iter);
            res = iter.get_res();

            if (iter.converged() && (res < 1e3)) {
              x = base::DataVector(xVec);
              Printer::getInstance().printStatusUpdate(
                "solving with Gmm++ (k = " +
                std::to_string(iter.get_iteration()) +
                ", residual norm = " + std::to_string(res) + ")");
              Printer::getInstance().printStatusEnd();
              return true;
            } else {
              Printer::getInstance().printStatusEnd(
                "error: could not solve linear system, "
                "method didn't converge");
              return false;
            }
          }
        } catch (std::exception& e) {
          Printer::getInstance().printStatusEnd(
            "error: could not solve linear system, what(): " +
            std::string(e.what()));
          return false;
        }
      }
#endif /* USE_GMMPP */

      bool Gmmpp::solve(SLE& system, base::DataVector& b,
                        base::DataVector& x) const {
#ifdef USE_GMMPP
        Printer::getInstance().printStatusBegin("Solving linear system (Gmm++)...");

        const size_t n = system.getDimension();
        size_t nnz = 0;
        gmm::csr_matrix<float_t> A2;

        {
          gmm::row_matrix<gmm::rsvector<float_t>> A(n, n);

          // parallelize only if the system is cloneable
          #pragma omp parallel if (system.isCloneable()) \
          shared(system, A, nnz, printer) default(none)
          {
            SLE* system2 = &system;
#ifdef _OPENMP
            std::unique_ptr<CloneableSLE> clonedSLE;

            if (system.isCloneable() && (omp_get_max_threads() > 1)) {
              dynamic_cast<CloneableSLE&>(system).clone(clonedSLE);
              system2 = clonedSLE.get();
            }

#endif /* _OPENMP */

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
                  Printer::getInstance().printStatusUpdate("constructing sparse matrix (" +
                  std::string(str) + ")");
                }
              }
            }
          }

          // the Gmm++ manual said to do so... no idea if that's necessary
          gmm::clean(A, 1e-12);
          gmm::copy(A, A2);
        }

        Printer::getInstance().printStatusUpdate("constructing sparse matrix (100.0%)");
        Printer::getInstance().printStatusNewLine();

        // print ratio of nonzero entries
        {
          char str[10];
          float_t nnz_ratio = static_cast<float_t>(nnz) /
                              (static_cast<float_t>(n) *
                               static_cast<float_t>(n));
          snprintf(str, 10, "%.1f%%", nnz_ratio * 100.0);
          Printer::getInstance().printStatusUpdate("nnz ratio: " + std::string(str));
          Printer::getInstance().printStatusNewLine();
        }

        x.resize(n);
        bool result = solveInternal(A2, b, x);
        return result;
#else
        std::cerr << "Error in sle_solver::Gmmpp::solve: "
                  << "SG++ was compiled without Gmm++ support!\n";
        return false;
#endif /* USE_GMMPP */
      }

    }
  }
}
