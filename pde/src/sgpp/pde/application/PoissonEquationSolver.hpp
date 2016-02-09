// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef POISSONEQUATIONSOLVER_HPP
#define POISSONEQUATIONSOLVER_HPP


#include <sgpp/pde/application/EllipticPDESolver.hpp>

#include <sgpp/base/grid/type/LinearGrid.hpp>
#include <sgpp/base/grid/common/BoundingBox.hpp>

#include <sgpp/base/tools/StdNormalDistribution.hpp>

#include <sgpp/base/application/ScreenOutput.hpp>

#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <cmath>

#include <sgpp/globaldef.hpp>
#include "../../../../../base/src/sgpp/base/grid/type/LinearBoundaryGrid.hpp"


namespace SGPP {
namespace pde {

/**
 * This class provides a simple-to-use solver of the multi dimensional
 * Poisson Equation on Sparse Grids.
 *
 * The class's aim is, to hide all complex details of solving the
 * Poisson Equation on Sparse Grids!
 *
 */
class PoissonEquationSolver : public EllipticPDESolver {
 private:
  /// screen object used in this solver
  SGPP::base::ScreenOutput* myScreen;

 public:
  /**
   * Std-Constructor of the solver
   */
  PoissonEquationSolver();

  /**
   * Std-Destructor of the solver
   */
  virtual ~PoissonEquationSolver();

  void constructGrid(SGPP::base::BoundingBox& myBoundingBox, int level);

  void solvePDE(SGPP::base::DataVector& alpha, SGPP::base::DataVector& rhs,
                size_t maxCGIterations, float_t epsilonCG, bool verbose = false);

  /**
   * Inits the grid with a smooth heat distribution (based on
   * a std-normal distribution) on its boundaries
   *
   * Coefficients of inner grid points are set to zero
   * since an elliptical PDE is solved
   *
   * @param alpha reference to the coefficients vector
   * @param mu the exspected value of the normal distribution
   * @param sigma the sigma of the normal distribution
   * @param factor a factor that is used to stretch the function values
   */
  void initGridWithSmoothHeat(SGPP::base::DataVector& alpha, float_t mu,
                              float_t sigma, float_t factor);

  /**
   * Inits the grid with a smooth heat distribution (based on
   * a std-normal distribution) on its boundaries
   *
   * Coefficients of inner grid points aren't set to zero
   * since they are used to hint an adaptive refinement
   * of the grid BEFORE solving the PDE.
   *
   * @param alpha reference to the coefficients vector
   * @param mu the exspected value of the normal distribution
   * @param sigma the sigma of the normal distribution
   * @param factor a factor that is used to stretch the function values
   */
  void initGridWithSmoothHeatFullDomain(SGPP::base::DataVector& alpha, float_t mu,
                                        float_t sigma, float_t factor);

  /**
   * Inits the grid with a heat distribution based on
   * the e-function
   *
   * The e-function is shifted in that way the right boundary
   * values becomes 1 (in case of factor = 1)
   *
   * @param alpha reference to the coefficient's vector
   * @param factor a constant factor used to enlarge the exp-functions input parameter
   */
  void initGridWithExpHeat(SGPP::base::DataVector& alpha, float_t factor = 1.0);

  /**
   * Inits the grid with a heat distribution based on
   * the e-function
   *
   * The e-function is shifted in that way the right boundary
   * values becomes 1 (in case of factor = 1)
   *
   * @param alpha reference to the coefficient's vector
   * @param factor a constant factor used to enlarge the exp-functions input parameter
   */
  void initGridWithExpHeatFullDomain(SGPP::base::DataVector& alpha,
                                     float_t factor = 1.0);

  /**
   * Routine to export the RHS of the inner system which has to be
   * solved in order to solve the Poisson equation
   *
   * @param alpha the start solution
   * @param tFilename file into which the rhs is written
   */
  void storeInnerRHS(SGPP::base::DataVector& alpha, std::string tFilename);

  /**
   * Routine to export the solution of the inner system which
   * has been calculated by Up/Down scheme
   *
   * @param alpha the start solution
   * @param maxCGIterations the maximum of interation in the CG solver
   * @param epsilonCG the epsilon used in the C
   * @param tFilename file into which the rhs is written
   */
  void storeInnerSolution(SGPP::base::DataVector& alpha, size_t maxCGIterations,
                          float_t epsilonCG, std::string tFilename);

  /**
   * Inits the screen object
   */
  void initScreen();
};

}
}

#endif /* POISSONEQUATIONSOLVER_HPP */
