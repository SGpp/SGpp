// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/globaldef.hpp>

#include <sgpp/optimization/sle/solver/Eigen.hpp>
#include <sgpp/optimization/sle/system/CloneableSLE.hpp>
#include <sgpp/optimization/tools/Printer.hpp>

#ifdef USE_EIGEN
#include <eigen3/Eigen/Dense>
#endif /* USE_EIGEN */

#include <cstddef>
#include <iostream>

namespace SGPP {
namespace optimization {
namespace sle_solver {

#ifdef USE_EIGEN

#if USE_DOUBLE_PRECISION
typedef ::Eigen::VectorXd EigenVector;
typedef ::Eigen::MatrixXd EigenMatrix;
#else
typedef ::Eigen::VectorXf EigenVector;
typedef ::Eigen::MatrixXf EigenMatrix;
#endif

/**
 * @param       A     coefficient matrix
 * @param       A_QR  Householder QR decomposition of coefficient matrix
 * @param       b     right-hand side
 * @param[out]  x     solution of the linear system
 * @return            whether all went well (false if errors occurred)
 */
bool solveInternal(
  const EigenMatrix& A,
  const ::Eigen::HouseholderQR<EigenMatrix>& A_QR,
  base::DataVector& b,
  base::DataVector& x) {
  const SGPP::float_t tolerance =
#if USE_DOUBLE_PRECISION
    1e-12;
#else
    1e-4;
#endif

  // solve system
  EigenVector bEigen = EigenVector::Map(b.getPointer(), b.getSize());
  EigenVector xEigen = A_QR.solve(bEigen);

  // check solution
  if ((A * xEigen).isApprox(bEigen, tolerance)) {
    x = base::DataVector(xEigen.data(), xEigen.size());
    return true;
  } else {
    return false;
  }
}
#endif /* USE_EIGEN */

Eigen::~Eigen() {
}

bool Eigen::solve(SLE& system, base::DataVector& b,
                  base::DataVector& x) const {
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

bool Eigen::solve(SLE& system,
                  base::DataMatrix& B,
                  base::DataMatrix& X) const {
#ifdef USE_EIGEN
  Printer::getInstance().printStatusBegin("Solving linear system (Eigen)...");

  const size_t n = system.getDimension();
  EigenMatrix A = EigenMatrix::Zero(n, n);
  size_t nnz = 0;

  // parallelize only if the system is cloneable
  #pragma omp parallel if (system.isCloneable()) \
  shared(system, A, nnz) default(none)
  {
    SLE* system2 = &system;
#ifdef _OPENMP
    std::unique_ptr<CloneableSLE> clonedSLE;

    if (system.isCloneable() && (omp_get_max_threads() > 1)) {
      dynamic_cast<CloneableSLE&>(system).clone(clonedSLE);
      system2 = clonedSLE.get();
    }

#endif /* _OPENMP */

    // copy system matrix to Eigen matrix object
    #pragma omp for ordered schedule(dynamic)

    for (size_t i = 0; i < n; i++) {
      for (size_t j = 0; j < n; j++) {
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
          Printer::getInstance().printStatusUpdate("constructing matrix (" +
          std::string(str) + ")");
        }
      }
    }
  }

  Printer::getInstance().printStatusUpdate("constructing matrix (100.0%)");
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

  // calculate QR factorization of system matrix
  Printer::getInstance().printStatusUpdate("step 1: Householder QR factorization");
  ::Eigen::HouseholderQR<EigenMatrix> A_QR = A.householderQr();

  base::DataVector x(n);
  base::DataVector b(n);
  X.resize(n, B.getNcols());

  // solve system for each RHS
  for (size_t i = 0; i < B.getNcols(); i++) {
    B.getColumn(i, b);
    Printer::getInstance().printStatusNewLine();

    if (B.getNcols() == 1) {
      Printer::getInstance().printStatusUpdate("step 2: solving");
    } else {
      Printer::getInstance().printStatusUpdate("step 2: solving (RHS " +
          std::to_string(i + 1) +
          " of " + std::to_string(B.getNcols()) +
          ")");
    }

    if (solveInternal(A, A_QR, b, x)) {
      X.setColumn(i, x);
    } else {
      Printer::getInstance().printStatusEnd("error: Could not solve linear system!");
      return false;
    }
  }

  Printer::getInstance().printStatusEnd();
  return true;
#else
  std::cerr << "Error in sle_solver::Eigen::solve: "
            << "SG++ was compiled without Eigen support!\n";
  return false;
#endif /* USE_EIGEN */
}

}
}
}
