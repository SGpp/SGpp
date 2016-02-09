// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef HEATEQUATIONSOLVER_HPP
#define HEATEQUATIONSOLVER_HPP

#include <sgpp/pde/application/ParabolicPDESolver.hpp>

#include <sgpp/base/grid/type/LinearGrid.hpp>
#include <sgpp/base/grid/common/BoundingBox.hpp>

#include <sgpp/base/tools/StdNormalDistribution.hpp>

#include <sgpp/base/application/ScreenOutput.hpp>
#include <sgpp/base/tools/SGppStopwatch.hpp>

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
 * Heat Equation on Sparse Grids.
 *
 * The class's aim is, to hide all complex details of solving the
 * Heat Equation on Sparse Grids!
 *
 */
class HeatEquationSolver : public ParabolicPDESolver {
 protected:
  /// the heat coefficient
  float_t a;
  /// screen object used in this solver
  SGPP::base::ScreenOutput* myScreen;

 public:
  /**
   * Std-Constructor of the solver
   */
  HeatEquationSolver();

  /**
   * Std-Destructor of the solver
   */
  virtual ~HeatEquationSolver();

  void constructGrid(SGPP::base::BoundingBox& myBoundingBox, int level);

  virtual void solveExplicitEuler(size_t numTimesteps, float_t timestepsize,
                                  size_t maxCGIterations, float_t epsilonCG, SGPP::base::DataVector& alpha,
                                  bool verbose = false, bool generateAnimation = false,
                                  size_t numEvalsAnimation = 20);

  virtual void solveImplicitEuler(size_t numTimesteps, float_t timestepsize,
                                  size_t maxCGIterations, float_t epsilonCG, SGPP::base::DataVector& alpha,
                                  bool verbose = false, bool generateAnimation = false,
                                  size_t numEvalsAnimation = 20);

  virtual void solveCrankNicolson(size_t numTimesteps, float_t timestepsize,
                                  size_t maxCGIterations, float_t epsilonCG, SGPP::base::DataVector& alpha,
                                  size_t NumImEul = 0);

  /**
   * This method sets the heat coefficient of the regarded material
   *
   * @param a the heat coefficient
   */
  void setHeatCoefficient(float_t a);

  /**
   * Inits the grid with a smooth heat distribution based on the
   * normal distribution formula
   *
   * @param alpha reference to the coefficient's vector
   * @param mu the exspected value of the normal distribution
   * @param sigma the sigma of the normal distribution
   * @param factor a factor that is used to stretch the function values
   */
  void initGridWithSmoothHeat(SGPP::base::DataVector& alpha, float_t mu,
                              float_t sigma, float_t factor);

  /**
   * Inits the screen object
   */
  virtual void initScreen();

  /**
   * Routine to export the RHS of the inner system which has to be
   * solved in order to solve the Poisson equation
   *
   * @param alpha the start solution
   * @param tFilename file into which the rhs is written
   * @param timestepsize the size of the timesteps
   */
  void storeInnerRHS(SGPP::base::DataVector& alpha, std::string tFilename,
                     float_t timestepsize);

  /**
   * Routine to export the solution of the inner system which
   * has been calculated by Up/Down scheme
   *
   * @param alpha the start solution
   * @param numTimesteps number timesteps
   * @param timestepsize size of timesteps
   * @param maxCGIterations the maximum of interation in the CG solver
   * @param epsilonCG the epsilon used in the C
   * @param tFilename file into which the rhs is written
   */
  void storeInnerSolution(SGPP::base::DataVector& alpha, size_t numTimesteps,
                          float_t timestepsize, size_t maxCGIterations, float_t epsilonCG,
                          std::string tFilename);
};

}
}

#endif /* HEATEQUATIONSOLVER_HPP */
