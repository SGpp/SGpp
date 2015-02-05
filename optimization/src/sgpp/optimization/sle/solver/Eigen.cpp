// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/globaldef.hpp>

#include <sgpp/optimization/sle/solver/Eigen.hpp>
#include <sgpp/optimization/sle/system/Cloneable.hpp>
#include <sgpp/optimization/tools/Printer.hpp>

#ifdef USEEIGEN
#include <eigen3/Eigen/Dense>
#endif

#include <cstddef>
#include <iostream>

namespace SGPP {
  namespace optimization {
    namespace sle {
      namespace solver {

#ifdef USEEIGEN
        /**
         * @param       A       coefficient matrix
         * @param       A_QR    Householder QR decomposition of coefficient matrix
         * @param       b       right-hand side
         * @param[out]  x       solution of the linear system
         * @return              whether all went well (false if errors occurred)
         */
        bool solveInternal(const ::Eigen::MatrixXd& A,
                           const ::Eigen::HouseholderQR< ::Eigen::MatrixXd>& A_QR,
                           const std::vector<float_t>& b,
                           std::vector<float_t>& x) {
          // solve system
          ::Eigen::VectorXd bEigen = ::Eigen::VectorXd::Map(&b[0], b.size());
          ::Eigen::VectorXd xEigen = A_QR.solve(bEigen);

          // check solution
          if ((A * xEigen).isApprox(bEigen)) {
            x = std::vector<float_t>(xEigen.data(), xEigen.data() + xEigen.size());
            return true;
          } else {
            return false;
          }
        }
#endif

        bool Eigen::solve(system::System& system, const std::vector<float_t>& b,
                          std::vector<float_t>& x) const {
          std::vector<std::vector<float_t> > B;
          std::vector<std::vector<float_t> > X;
          B.push_back(b);

          if (solve(system, B, X)) {
            x = X[0];
            return true;
          } else {
            return false;
          }
        }

        bool Eigen::solve(system::System& system, const std::vector<std::vector<float_t> >& B,
                          std::vector<std::vector<float_t> >& X) const {
#ifdef USEEIGEN
          tools::printer.printStatusBegin("Solving linear system (Eigen)...");

          const size_t n = system.getDimension();
          ::Eigen::MatrixXd A = ::Eigen::MatrixXd::Zero(n, n);
          size_t nnz = 0;

          // parallelize only if the system is cloneable
          #pragma omp parallel shared(system, A, nnz, tools::printer) default(none)
          {
            system::System* system2 = &system;
#ifdef _OPENMP
            std::unique_ptr<system::System> clonedSystem;

            if (omp_get_max_threads() > 1) {
              clonedSystem = std::unique_ptr<system::System>(
                               dynamic_cast<system::Cloneable&>(system).clone());
              system2 = clonedSystem.get();
            }

#endif

            // copy system matrix to Eigen matrix object
            #pragma omp for ordered schedule(dynamic)

            for (size_t i = 0; i < n; i++) {
              for (size_t j = 0; j < n; j++) {
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
            float_t nnz_ratio = static_cast<float_t>(nnz) /
                                (static_cast<float_t>(n) * static_cast<float_t>(n));
            snprintf(str, 10, "%.1f%%", nnz_ratio * 100.0);
            tools::printer.printStatusUpdate("nnz ratio: " + std::string(str));
            tools::printer.printStatusNewLine();
          }

          // calculate QR factorization of system matrix
          tools::printer.printStatusUpdate("step 1: Householder QR factorization");
          ::Eigen::HouseholderQR< ::Eigen::MatrixXd> A_QR = A.householderQr();

          std::vector<float_t> x;
          X.clear();

          // solve system for each RHS
          for (size_t i = 0; i < B.size(); i++) {
            const std::vector<float_t>& b = B[i];
            tools::printer.printStatusNewLine();

            if (B.size() == 1) {
              tools::printer.printStatusUpdate("step 2: solving");
            } else {
              tools::printer.printStatusUpdate("step 2: solving (RHS " + toString(i + 1) +
                                               " of " + toString(B.size()) + ")");
            }

            if (solveInternal(A, A_QR, b, x)) {
              X.push_back(x);
            } else {
              tools::printer.printStatusEnd("error: could not solve linear system!");
              return false;
            }
          }

          tools::printer.printStatusEnd();
          return true;
#else
          std::cerr << "Error in SGPP::optimization::sle::solver::Eigen::solve: "
                    << "SG++ was compiled without Eigen support!\n";
          return false;
#endif
        }

      }
    }
  }
}
