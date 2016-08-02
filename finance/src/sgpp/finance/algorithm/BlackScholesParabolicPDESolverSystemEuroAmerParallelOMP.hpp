// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef BLACKSCHOLESPARABOLICPDESOLVERSYSTEMEUROAMERPARALLELOMP_HPP
#define BLACKSCHOLESPARABOLICPDESOLVERSYSTEMEUROAMERPARALLELOMP_HPP

#include <sgpp/finance/algorithm/BlackScholesParabolicPDESolverSystemEuroAmer.hpp>

#include <sgpp/globaldef.hpp>

#include <string>

namespace sgpp {
namespace finance {

/**
 * This class implements the ParabolicPDESolverSystem for the BlackScholes
 * Equation.
 *
 *
 * Here European or American Options with fix Dirichlet boundaries are solved.
 *
 * It's derived from the existing BlackScholesParabolicPDESolverSystemEuropean but uses
 * the OMP task concept to enable further parallelization possibilities
 * in the calculation of the space-discretization operator (L)
 */
class BlackScholesParabolicPDESolverSystemEuroAmerParallelOMP
    : public BlackScholesParabolicPDESolverSystemEuroAmer {
 protected:
  virtual void applyLOperatorInner(sgpp::base::DataVector& alpha, sgpp::base::DataVector& result);

  virtual void applyLOperatorComplete(sgpp::base::DataVector& alpha,
                                      sgpp::base::DataVector& result);

  virtual void applyMassMatrixInner(sgpp::base::DataVector& alpha, sgpp::base::DataVector& result);

  virtual void applyMassMatrixComplete(sgpp::base::DataVector& alpha,
                                       sgpp::base::DataVector& result);

 public:
  /**
   * Std-Constructor
   *
   * @param SparseGrid reference to the sparse grid
   * @param alpha the ansatzfunctions' coefficients
   * @param mu reference to the mus
   * @param sigma reference to the sigmas
   * @param rho reference to the rhos
   * @param r the riskfree interest rate
   * @param TimestepSize the size of one timestep used in the ODE Solver
   * @param OperationMode specifies in which solver this matrix is used, valid values are: ExEul for
   * explicit Euler,
   *                ImEul for implicit Euler, CrNic for Crank Nicolson solver
   * @param dStrike strike
   * @param option_type type of option
   * @param bLogTransform indicates that this system belongs to a log-transformed Black Scholes
   * Equation
   * @param useCoarsen specifies if the grid should be coarsened between timesteps
   * @param coarsenThreshold Threshold to decide, if a grid point should be deleted
   * @param adaptSolveMode adaptive mode during solving: coarsen, refine, coarsenNrefine
   * @param numCoarsenPoints number of point that should be coarsened in one coarsening step
   * !CURRENTLY UNUSED PARAMETER!
   * @param refineThreshold Threshold to decide, if a grid point should be refined
   * @param refineMode refineMode during solving Black Scholes Equation: classic or maxLevel
   * @param refineMaxLevel max. level for refinement during solving
   */
  BlackScholesParabolicPDESolverSystemEuroAmerParallelOMP(
      sgpp::base::Grid& SparseGrid, sgpp::base::DataVector& alpha, sgpp::base::DataVector& mu,
      sgpp::base::DataVector& sigma, sgpp::base::DataMatrix& rho, double r, double TimestepSize,
      std::string OperationMode, double dStrike, std::string option_type,
      bool bLogTransform = false, bool useCoarsen = false, double coarsenThreshold = 0.0,
      std::string adaptSolveMode = "none", int numCoarsenPoints = -1, double refineThreshold = 0.0,
      std::string refineMode = "classic", sgpp::base::GridPoint::level_type refineMaxLevel = 0);

  /**
   * Std-Destructor
   */
  virtual ~BlackScholesParabolicPDESolverSystemEuroAmerParallelOMP();

  /**
   * Multiplicates a vector with the matrix, parallel
   *
   * @param alpha sgpp::base::DataVector that contains the ansatzfunctions' coefficients
   * @param result sgpp::base::DataVector into which the result of the space discretization
   * operation is stored
   */
  virtual void mult(sgpp::base::DataVector& alpha, sgpp::base::DataVector& result);

  /**
   * generates the right hand side of the system, parallel
   *
   * @return returns the rhs
   */
  virtual sgpp::base::DataVector* generateRHS();
};
}  // namespace finance
}  // namespace sgpp

#endif /* BLACKSCHOLESPARABOLICPDESOLVERSYSTEMEUROAMERPARALLELOMP_HPP */
