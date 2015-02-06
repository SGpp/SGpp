// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/globaldef.hpp>

#include <sgpp/optimization/sle/solver/Armadillo.hpp>
#include <sgpp/optimization/sle/system/Cloneable.hpp>
#include <sgpp/optimization/tools/Printer.hpp>

#ifdef USEARMADILLO
#include <armadillo>
#endif

#include <cstddef>
#include <iostream>

namespace SGPP {
  namespace optimization {
    namespace sle {
      namespace solver {

        bool Armadillo::solve(system::System& system, const std::vector<float_t>& b,
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

        bool Armadillo::solve(system::System& system, const std::vector<std::vector<float_t> >& B,
                              std::vector<std::vector<float_t> >& X) const {
#ifdef USEARMADILLO
          tools::printer.printStatusBegin("Solving linear system (Armadillo)...");

          const arma::uword n = static_cast<arma::uword>(system.getDimension());
          arma::mat A(n, n);
          size_t nnz = 0;

          A.zeros();

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

            // copy system matrix to Armadillo matrix object
            #pragma omp for ordered schedule(dynamic)

            for (arma::uword i = 0; i < n; i++) {
              for (arma::uword j = 0; j < n; j++) {
                A(i, j) = system2->getMatrixEntry(i, j);

                // count nonzero entries (not necessary, you can also remove that if you like)
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
                  tools::printer.printStatusUpdate("constructing matrix (" +
                  std::string(str) + ")");
                }
              }
            }
          }

          tools::printer.printStatusUpdate("constructing matrix (100.0%)");
          tools::printer.printStatusNewLine();

          // print ratio of nonzero entries
          {
            char str[10];
            float_t nnzRatio = static_cast<float_t>(nnz) /
                               (static_cast<float_t>(n) * static_cast<float_t>(n));
            snprintf(str, 10, "%.1f%%", nnzRatio * 100.0);
            tools::printer.printStatusUpdate("nnz ratio: " + std::string(str));
            tools::printer.printStatusNewLine();
          }

          if (B.size() == 1) {
            // only one RHS ==> use vector version of arma::solve
            arma::vec b_armadillo = arma::conv_to<arma::vec>::from(B[0]);
            arma::vec x_armadillo(n);

            tools::printer.printStatusUpdate("solving with Armadillo");

            if (arma::solve(x_armadillo, A, b_armadillo)) {
              X.clear();
              X.push_back(arma::conv_to<std::vector<float_t> >::from(x_armadillo));
              tools::printer.printStatusEnd();
              return true;
            } else {
              tools::printer.printStatusEnd("error: could not solve linear system!");
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

            tools::printer.printStatusUpdate("solving with Armadillo");

            if (arma::solve(XArmadillo, A, BArmadillo)) {
              X.clear();

              // convert solutions to std::vector
              for (arma::uword i = 0; i < B_count; i++) {
                X.push_back(arma::conv_to<std::vector<float_t> >::from(XArmadillo.col(i)));
              }

              tools::printer.printStatusEnd();
              return true;
            } else {
              tools::printer.printStatusEnd("error: could not solve linear system!");
              return false;
            }
          }

#else
          std::cerr << "Error in SGPP::optimization::sle::solver::Armadillo::solve: "
                    << "SG++ was compiled without Armadillo support!\n";
          return false;
#endif
        }

      }
    }
  }
}
