/* ****************************************************************************
* Copyright (C) 2014 Technische Universitaet Muenchen                         *
* This file is part of the SG++ project. For conditions of distribution and   *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp          *
**************************************************************************** */
// @author Julian Valentin (julian.valentin@stud.mathematik.uni-stuttgart.de)

#include "opt/sle/solver/Armadillo.hpp"
#include "opt/sle/system/Cloneable.hpp"
#include "opt/tools/Printer.hpp"

#ifdef USEARMADILLO
#include <armadillo>
#endif

#include <cstddef>
#include <iostream>

namespace sg {
  namespace opt {
    namespace sle {
      namespace solver {

        bool Armadillo::solve(system::System& system, const std::vector<double>& b,
                              std::vector<double>& x) const {
          // place RHS in its own vector
          std::vector<std::vector<double> > B;
          std::vector<std::vector<double> > X;
          B.push_back(b);

          // call version for multiple RHSs
          if (solve(system, B, X)) {
            x = X[0];
            return true;
          } else {
            return false;
          }
        }

        bool Armadillo::solve(system::System& system, const std::vector<std::vector<double> >& B,
                              std::vector<std::vector<double> >& X) const {
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
            system::System* system2;
            tools::SmartPointer<system::System> cloned_system;

            if (system.isCloneable()) {
              cloned_system = dynamic_cast<system::Cloneable&>(system).clone();
              system2 = cloned_system.get();
            } else {
              system2 = &system;
            }

            // copy system matrix to Armadillo matrix object
            #pragma omp for ordered schedule(dynamic)
            for (arma::uword i = 0; i < n; i++) {
              for (arma::uword j = 0; j < n; j++) {
                A(i,j) = system2->getMatrixEntry(i, j);

                // count nonzero entries (not necessary, you can also remove that if you like)
                if (A(i,j) != 0) {
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
                  static_cast<double>(i) / static_cast<double>(n) * 100.0);
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
            double nnz_ratio = static_cast<double>(nnz) /
                               (static_cast<double>(n) * static_cast<double>(n));
            snprintf(str, 10, "%.1f%%", nnz_ratio * 100.0);
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
              X.push_back(arma::conv_to<std::vector<double> >::from(x_armadillo));
              tools::printer.printStatusEnd();
              return true;
            } else {
              tools::printer.printStatusEnd("error: could not solve linear system!");
              return false;
            }
          } else {
            // multiple RHSs ==> use matrix version of arma::solve
            const arma::uword B_count = static_cast<arma::uword>(B.size());
            arma::mat B_armadillo(n, B_count);
            arma::mat X_armadillo(n, B_count);

            // copy RHSs to Armadillo matrix
            for (arma::uword i = 0; i < B_count; i++) {
              B_armadillo.col(i) = arma::conv_to<arma::vec>::from(B[i]);
            }

            tools::printer.printStatusUpdate("solving with Armadillo");

            if (arma::solve(X_armadillo, A, B_armadillo)) {
              X.clear();

              // convert solutions to std::vector
              for (arma::uword i = 0; i < B_count; i++) {
                X.push_back(arma::conv_to<std::vector<double> >::from(X_armadillo.col(i)));
              }

              tools::printer.printStatusEnd();
              return true;
            } else {
              tools::printer.printStatusEnd("error: could not solve linear system!");
              return false;
            }
          }
#else
          std::cerr << "Error in sg::opt::sle::solver::Armadillo::solve: "
                    << "SG++ was compiled without Armadillo support!\n";
          return false;
#endif
        }

      }
    }
  }
}
