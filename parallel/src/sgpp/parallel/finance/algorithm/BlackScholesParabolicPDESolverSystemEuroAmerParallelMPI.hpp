// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef BLACKSCHOLESPARABOLICPDESOLVERSYSTEMEUROAMERPARALLELMPI_HPP
#define BLACKSCHOLESPARABOLICPDESOLVERSYSTEMEUROAMERPARALLELMPI_HPP

#include <sgpp/finance/algorithm/BlackScholesParabolicPDESolverSystemEuroAmer.hpp>

#include <sgpp/globaldef.hpp>

#include <string>

namespace SGPP {
namespace parallel {

/**
 * This class implements the ParabolicPDESolverSystem for the BlackScholes
 * Equation.
 *
 * Here a European Option with fix Dirichlet boundaries is solved.
 *
 * It's derived from the existing BlackScholesParabolicPDESolverSystemEuropean but uses
 * the OMP task concept to enable further parallelization possibilities
 * in the calculation of the space-discretization operator (L)
 *
 * Parallelization of the FEM operators is done by using MPI.
 */
class BlackScholesParabolicPDESolverSystemEuroAmerParallelMPI
    : public SGPP::finance::BlackScholesParabolicPDESolverSystemEuroAmer {
 protected:
  virtual void applyLOperatorInner(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result);

  virtual void applyLOperatorComplete(SGPP::base::DataVector& alpha,
                                      SGPP::base::DataVector& result);

  virtual void applyMassMatrixInner(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result);

  virtual void applyMassMatrixComplete(SGPP::base::DataVector& alpha,
                                       SGPP::base::DataVector& result);

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
   * @param dStrike Strike value
   * @param option_type Type of option
   */
  BlackScholesParabolicPDESolverSystemEuroAmerParallelMPI(
      SGPP::base::Grid& SparseGrid, SGPP::base::DataVector& alpha, SGPP::base::DataVector& mu,
      SGPP::base::DataVector& sigma, SGPP::base::DataMatrix& rho, double r, double TimestepSize,
      std::string OperationMode, double dStrike, std::string option_type,
      bool bLogTransform = false, bool useCoarsen = false, double coarsenThreshold = 0.0,
      std::string adaptSolveMode = "none", int numCoarsenPoints = -1, double refineThreshold = 0.0,
      std::string refineMode = "classic", SGPP::base::GridIndex::level_type refineMaxLevel = 0);

  /**
   * Std-Destructor
   */
  virtual ~BlackScholesParabolicPDESolverSystemEuroAmerParallelMPI();

  /**
   * Multiplicates a vector with the matrix, parallel
   *
   * @param alpha DataVector that contains the ansatzfunctions' coefficients
   * @param result DataVector into which the result of the space discretization operation is stored
   */
  virtual void mult(SGPP::base::DataVector& alpha, SGPP::base::DataVector& result);

  /**
   * generates the right hand side of the system, parallel
   *
   * @return returns the rhs
   */
  virtual SGPP::base::DataVector* generateRHS();

  void finishTimestep(bool isLastTimestep = false);
};
}  // namespace parallel
}  // namespace SGPP

#endif /* BLACKSCHOLESPARABOLICPDESOLVERSYSTEMEUROAMERPARALLELMPI_HPP */
