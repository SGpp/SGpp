// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/tools/sle/system/SLE.hpp>
#include <sgpp/globaldef.hpp>

#include <vector>

namespace sgpp {
namespace base {
namespace sle_solver {

/**
 * Abstract class for solving systems of linear equations.
 */
class SLESolver {
 public:
  /**
   * Constructor.
   */
  SLESolver() {}

  /**
   * Destructor.
   */
  virtual ~SLESolver() {}

  /**
   * Pure virtual method for a solving linear system.
   *
   * @param       system  system to be solved
   * @param       b       right-hand side
   * @param[out]  x       solution to the system
   * @return              whether all went well
   *                      (false if errors occurred)
   */
  virtual bool solve(SLE& system, DataVector& b, DataVector& x) const = 0;

  /**
   * Virtual method for solving multiple linear systems with
   * different right-hand sides.
   * Defaults to calling the solve() method for a single
   * right-hand side multiple times.
   *
   * @param       system  system to be solved
   * @param       B       matrix of right-hand sides
   * @param[out]  X       matrix of solutions to the systems
   * @return              whether all went well
   *                      (false if errors occurred)
   */
  virtual bool solve(SLE& system, DataMatrix& B, DataMatrix& X) const {
    const size_t n = system.getDimension();
    const size_t m = B.getNcols();
    DataVector b(n);
    DataVector x(n);
    X.resize(n, m);

    for (size_t i = 0; i < m; i++) {
      B.getColumn(i, b);

      if (solve(system, b, x)) {
        X.setColumn(i, x);
      } else {
        return false;
      }
    }

    return true;
  }
};
}  // namespace sle_solver
}  // namespace base
}  // namespace sgpp
