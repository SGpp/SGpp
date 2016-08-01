// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef BLACKSCHOLESPATPARABOLICPDESOLVERSYSTEMEUROAMER_HPP
#define BLACKSCHOLESPATPARABOLICPDESOLVERSYSTEMEUROAMER_HPP

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/pde/operation/hash/OperationParabolicPDESolverSystemDirichlet.hpp>

#include <sgpp/globaldef.hpp>

#include <vector>
#include <string>

namespace sgpp {
namespace finance {

/**
 * This class implements the ParabolicPDESolverSystem for the BlackScholes
 * Equation.
 *
 * Here a European and American Options with fix Dirichlet boundaries are solved.
 *
 * In addition a principal axis transformation is performed in order
 * to reduce the computational effort.
 */
class BlackScholesPATParabolicPDESolverSystemEuroAmer
    : public sgpp::pde::OperationParabolicPDESolverSystemDirichlet {
 protected:
  /// the Laplace Operation, on boundary grid
  sgpp::base::OperationMatrix* OpLaplaceBound;
  /// the LTwoDotProduct Operation (Mass Matrix), on boundary grid
  sgpp::base::OperationMatrix* OpLTwoBound;
  /// the Laplace Operation, on Inner grid
  sgpp::base::OperationMatrix* OpLaplaceInner;
  /// the LTwoDotProduct Operation (Mass Matrix), on Inner grid
  sgpp::base::OperationMatrix* OpLTwoInner;
  /// Pointer to the lambda (coefficients of the needed Laplace operator)
  sgpp::base::DataVector* lambda;
  /// Eigenvectors of the covariance matrix
  sgpp::base::DataMatrix* eigenvecs;
  /// Pointer to the mu_hat (transformed drifts and correlation, needed for constraint of American
  /// options)
  sgpp::base::DataVector* mu_hat;
  /// use coarsening between timesteps in order to reduce gridsize
  bool useCoarsen;
  /// adaptive mode during solving Black Scholes Equation: coarsen, refine, coarsenNrefine
  std::string adaptSolveMode;
  /// number of points the are coarsened in each coarsening-step !CURRENTLY UNUSED PARAMETER!
  int numCoarsenPoints;
  /// Threshold used to decide if a grid point should be deleted
  double coarsenThreshold;
  /// Threshold used to decide if a grid point should be refined
  double refineThreshold;
  /// refine mode during solving Black Scholes Equation: classic or maxLevel
  std::string refineMode;
  /// maxLevel max. Level of refinement
  sgpp::base::GridPoint::level_type refineMaxLevel;
  /// the algorithmic dimensions used in this system
  std::vector<size_t> BSalgoDims;
  /// store number of executed timesteps
  size_t nExecTimesteps;
  /// the strike of the current option
  double dStrike;
  /// the type of the current option
  std::string option_type;
  /// teh riskfree rate
  double r;

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
   * @param lambda reference to the lambdas
   * @param eigenvecs reference to the eigenvectors of the co-variance matrix
   * @param mu_hat reference to transformed drifts and correlation, needed for constraint of
   * American options
   * @param TimestepSize the size of one timestep used in the ODE Solver
   * @param OperationMode specifies in which solver this matrix is used, valid values are: ExEul for
   * explicit Euler,
   *                ImEul for implicit Euler, CrNic for Crank Nicolson solver
   * @param dStrike the strike of the option
   * @param option_type the option to be solved std_amer_put as an special case
   * @param r the riskfree rate
   * @param useCoarsen specifies if the grid should be coarsened between timesteps
   * @param coarsenThreshold Threshold to decide, if a grid point should be deleted
   * @param adaptSolveMode adaptive mode during solving: coarsen, refine, coarsenNrefine
   * @param numCoarsenPoints number of point that should be coarsened in one coarsening step
   * !CURRENTLY UNUSED PARAMETER!
   * @param refineThreshold Threshold to decide, if a grid point should be refined
   * @param refineMode refineMode during solving Black Scholes Equation: classic or maxLevel
   * @param refineMaxLevel max. level of refinement during solving
   */
  BlackScholesPATParabolicPDESolverSystemEuroAmer(
      sgpp::base::Grid& SparseGrid, sgpp::base::DataVector& alpha, sgpp::base::DataVector& lambda,
      sgpp::base::DataMatrix& eigenvecs, sgpp::base::DataVector& mu_hat, double TimestepSize,
      std::string OperationMode, double dStrike, std::string option_type, double r,
      bool useCoarsen = false, double coarsenThreshold = 0.0, std::string adaptSolveMode = "none",
      int numCoarsenPoints = -1, double refineThreshold = 0.0, std::string refineMode = "classic",
      sgpp::base::GridPoint::level_type refineMaxLevel = 0);

  /**
   * Std-Destructor
   */
  virtual ~BlackScholesPATParabolicPDESolverSystemEuroAmer();

  virtual void finishTimestep();

  /**
   * @param isLastTimestep specify if last time step
   */
  virtual void coarsenAndRefine(bool isLastTimestep = false);

  void startTimestep();
};
}  // namespace finance
}  // namespace sgpp

#endif /* BLACKSCHOLESPATPARABOLICPDESOLVERSYSTEMEUROAMER_HPP */
