// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/tools/sle/solver/SLESolver.hpp>
#include <sgpp/globaldef.hpp>

#include <cstddef>
#include <vector>

namespace sgpp {
namespace base {
namespace sle_solver {

/**
 * Linear system solver implementing the iterative BiCGStab method.
 */
class BiCGStab : public SLESolver {
 public:
  /// default maximal number of iterations
  static const size_t DEFAULT_MAX_IT_COUNT = 1000;
  /// default tolerance
  static constexpr double DEFAULT_TOLERANCE = 1e-10;

  /**
   * Constructor.
   */
  BiCGStab();

  /**
   * @param maxItCount        maximal number of iterations
   * @param tolerance         tolerance
   * @param startingPoint     starting vector
   */
  BiCGStab(size_t maxItCount, double tolerance, const base::DataVector& startingPoint);

  /**
   * Destructor.
   */
  ~BiCGStab() override;

  /**
   * @param       system  system to be solved
   * @param       b       right-hand side
   * @param[out]  x       solution to the system
   * @return              whether all went well
   *                      (false if errors occurred)
   */
  bool solve(SLE& system, base::DataVector& b, base::DataVector& x) const override;

  /**
   * @return              maximal number of iterations
   */
  size_t getMaxItCount() const;

  /**
   * @param maxItCount   maximal number of iterations
   */
  void setMaxItCount(size_t maxItCount);

  /**
   * @return              tolerance
   */
  double getTolerance() const;

  /**
   * @param tolerance     tolerance
   */
  void setTolerance(double tolerance);

  /**
   * @return                  starting vector
   */
  const base::DataVector& getStartingPoint() const;

  /**
   * @param startingPoint     starting vector
   */
  void setStartingPoint(const base::DataVector& startingPoint);

 protected:
  /// maximal number of iterations
  size_t N;
  /// tolerance
  double tol;
  /// starting vector
  base::DataVector x0;
};
}  // namespace sle_solver
}  // namespace base
}  // namespace sgpp
