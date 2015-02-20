// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/globaldef.hpp>

#include <sgpp/optimization/sle/solver/Armadillo.hpp>
#include <sgpp/optimization/sle/system/CloneableSLE.hpp>
#include <sgpp/optimization/tools/Printer.hpp>

#ifdef USEARMADILLO
#include <armadillo>
#endif /* USEARMADILLO */

#include <cstddef>
#include <iostream>

namespace SGPP {
  namespace optimization {
    namespace sle_solver {

      bool Armadillo::solve(SLE& system, const std::vector<float_t>& b,
                            std::vector<float_t>& x) const {
        // place RHS in its own vector
        std::vector<std::vector<float_t>> B;
        std::vector<std::vector<float_t>> X;
        B.push_back(b);

        // call version for multiple RHSs
        if (solve(system, B, X)) {
          x = X[0];
          return true;
        } else {
          return false;
        }
      }

      bool Armadillo::solve(SLE& system,
                            const std::vector<std::vector<float_t>>& B,
                            std::vector<std::vector<float_t>>& X) const {
#ifdef USEARMADILLO
        printer.printStatusBegin("Solving linear system (Armadillo)...");

        const arma::uword n = static_cast<arma::uword>(system.getDimension());
        arma::mat A(n, n);
        size_t nnz = 0;

        A.zeros();

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

          // copy system matrix to Armadillo matrix object
          #pragma omp for ordered schedule(dynamic)

          for (arma::uword i = 0; i < n; i++) {
            for (arma::uword j = 0; j < n; j++) {
              A(i, j) = system2->getMatrixEntry(i, j);

              // count nonzero entries
              // (not necessary, you can also remove that if you like)
              if (A(i, j) != 0) {
                #pragma omp atomic
                nnz++;
              }
            }

            // status message
            if (i % 100 == 0) {
              #pragma omp ordered
              {
                char str[10];
                snprintf(str, 10, "%.1f%%",
                static_cast<float_t>(i) / static_cast<float_t>(n) * 100.0);
                printer.printStatusUpdate("constructing matrix (" +
                std::string(str) + ")");
              }
            }
          }
        }

        printer.printStatusUpdate("constructing matrix (100.0%)");
        printer.printStatusNewLine();

        // print ratio of nonzero entries
        {
          char str[10];
          float_t nnzRatio = static_cast<float_t>(nnz) /
                             (static_cast<float_t>(n) *
                                 static_cast<float_t>(n));
          snprintf(str, 10, "%.1f%%", nnzRatio * 100.0);
          printer.printStatusUpdate("nnz ratio: " + std::string(str));
          printer.printStatusNewLine();
        }

        if (B.size() == 1) {
          // only one RHS ==> use vector version of arma::solve
          arma::vec b_armadillo = arma::conv_to<arma::vec>::from(B[0]);
          arma::vec x_armadillo(n);

          printer.printStatusUpdate("solving with Armadillo");

          if (arma::solve(x_armadillo, A, b_armadillo)) {
            X.clear();
            X.push_back(
                arma::conv_to<std::vector<float_t>>::from(x_armadillo));
            printer.printStatusEnd();
            return true;
          } else {
            printer.printStatusEnd("error: could not solve linear system!");
            return false;
          }
        } else {
          // multiple RHSs ==> use matrix version of arma::solve
          const arma::uword B_count = static_cast<arma::uword>(B.size());
          arma::mat BArmadillo(n, B_count);
          arma::mat XArmadillo(n, B_count);

          // copy RHSs to Armadillo matrix
          for (arma::uword i = 0; i < B_count; i++) {
            BArmadillo.col(i) = arma::conv_to<arma::vec>::from(B[i]);
          }

          printer.printStatusUpdate("solving with Armadillo");

          if (arma::solve(XArmadillo, A, BArmadillo)) {
            X.clear();

            // convert solutions to std::vector
            for (arma::uword i = 0; i < B_count; i++) {
              X.push_back(arma::conv_to<std::vector<float_t>>::from(
                  XArmadillo.col(i)));
            }

            printer.printStatusEnd();
            return true;
          } else {
            printer.printStatusEnd("error: could not solve linear system!");
            return false;
          }
        }

#else
        std::cerr << "Error in sle_solver::Armadillo::solve: "
                  << "SG++ was compiled without Armadillo support!\n";
        return false;
#endif /* USEARMADILLO */
      }

    }
  }
}
