// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/tools/Printer.hpp>
#include <sgpp/base/tools/sle/solver/Armadillo.hpp>
#include <sgpp/base/tools/sle/system/CloneableSLE.hpp>
#include <sgpp/globaldef.hpp>

#ifdef USE_ARMADILLO
#include <armadillo>

typedef arma::vec ArmadilloVector;
typedef arma::mat ArmadilloMatrix;
#endif /* USE_ARMADILLO */

#include <cstddef>
#include <iostream>
#include <string>

namespace sgpp {
namespace base {
namespace sle_solver {

Armadillo::~Armadillo() {}

bool Armadillo::solve(SLE& system, DataVector& b, DataVector& x) const {
  DataMatrix B(b.getPointer(), b.getSize(), 1);
  DataMatrix X(B.getNrows(), B.getNcols());

  // call version for multiple RHSs
  if (solve(system, B, X)) {
    x.resize(X.getNrows());
    X.getColumn(0, x);
    return true;
  } else {
    return false;
  }
}

bool Armadillo::solve(SLE& system, DataMatrix& B, DataMatrix& X) const {
#ifdef USE_ARMADILLO
  Printer::getInstance().printStatusBegin("Solving linear system (Armadillo)...");

  const arma::uword n = static_cast<arma::uword>(system.getDimension());
  ArmadilloMatrix A(n, n);
  size_t nnz = 0;
  size_t rowsDone = 0;

  A.zeros();

// parallelize only if the system is cloneable
#pragma omp parallel if (system.isCloneable()) shared(system, A, nnz, rowsDone)  // default(none)

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
#pragma omp for ordered schedule(static)

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

#pragma omp atomic
      rowsDone++;

      // status message
      if (rowsDone % 100 == 0) {
        char str[10];
        snprintf(str, sizeof(str), "%.1f%%",
                 static_cast<double>(rowsDone) / static_cast<double>(n) * 100.0);
        Printer::getInstance().printStatusUpdate("constructing matrix (" + std::string(str) + ")");
      }
    }
  }

  Printer::getInstance().printStatusUpdate("constructing matrix (100.0%)");
  Printer::getInstance().printStatusNewLine();

  // print ratio of nonzero entries
  {
    char str[10];
    double nnzRatio = static_cast<double>(nnz) / (static_cast<double>(n) * static_cast<double>(n));
    snprintf(str, sizeof(str), "%.1f%%", nnzRatio * 100.0);
    Printer::getInstance().printStatusUpdate("nnz ratio: " + std::string(str));
    Printer::getInstance().printStatusNewLine();
  }

  if (B.getNcols() == 1) {
    // only one RHS ==> use vector version of arma::solve
    ArmadilloVector bArmadillo(B.getPointer(), n);
    ArmadilloVector xArmadillo(n);

    Printer::getInstance().printStatusUpdate("solving with Armadillo");

    if (arma::solve(xArmadillo, A, bArmadillo)) {
      DataVector x(xArmadillo.memptr(), n);
      X.resize(n, 1);
      X.setColumn(0, x);
      Printer::getInstance().printStatusEnd();
      return true;
    } else {
      Printer::getInstance().printStatusEnd("error: Could not solve linear system!");
      return false;
    }
  } else {
    // multiple RHSs ==> use matrix version of arma::solve
    const arma::uword B_count = static_cast<arma::uword>(B.getNcols());
    ArmadilloMatrix BArmadillo(n, B_count);
    ArmadilloMatrix XArmadillo(n, B_count);
    DataVector b(n);

    // copy RHSs to Armadillo matrix
    for (arma::uword i = 0; i < B_count; i++) {
      B.getColumn(i, b);
      BArmadillo.col(i) = ArmadilloVector(b.getPointer(), n);
    }

    Printer::getInstance().printStatusUpdate("solving with Armadillo");

    if (arma::solve(XArmadillo, A, BArmadillo)) {
      X.resize(n, B_count);

      // convert solutions to DataVector
      for (arma::uword i = 0; i < B_count; i++) {
        DataVector x(XArmadillo.colptr(i), n);
        X.setColumn(i, x);
      }

      Printer::getInstance().printStatusEnd();
      return true;
    } else {
      Printer::getInstance().printStatusEnd("error: Could not solve linear system!");
      return false;
    }
  }

#else
  std::cerr << "Error in sle_solver::Armadillo::solve: "
            << "SG++ was compiled without Armadillo support!\n";
  return false;
#endif /* USE_ARMADILLO */
}
}  // namespace sle_solver
}  // namespace base
}  // namespace sgpp
