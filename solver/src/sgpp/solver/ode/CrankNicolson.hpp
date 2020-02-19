// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef CRANKNICOLSON_HPP
#define CRANKNICOLSON_HPP

#include <sgpp/base/application/ScreenOutput.hpp>
#include <sgpp/solver/ODESolver.hpp>

#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace solver {

/**
 * This class implements the Crank-Nicolson method
 * for solving ordinary partial equations
 *
 * For solving the system of linear equations the
 * already implemented CG-method is used
 *
 */
class CrankNicolson : public ODESolver {
 private:
  /// Pointer to sgpp::base::ScreenOutput object
  sgpp::base::ScreenOutput* myScreen;

 public:
  /**
   * Std-Constructer
   *
   * @param nTimesteps number of maximum executed iterations
   * @param timestepSize the size of one timestep
   * @param screen possible pointer to a sgpp::base::ScreenOutput object
   */
  CrankNicolson(size_t nTimesteps, double timestepSize, sgpp::base::ScreenOutput* screen = nullptr);

  /**
   * Std-Destructor
   */
  virtual ~CrankNicolson();

  virtual void solve(SLESolver& LinearSystemSolver,
                     sgpp::solver::OperationParabolicPDESolverSystem& System,
                     bool bIdentifyLastStep = false, bool verbose = false);
};

}  // namespace solver
}  // namespace sgpp

#endif /* CRANKNICOLSON_HPP */
