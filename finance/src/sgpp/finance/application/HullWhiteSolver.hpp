// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef HULLWHITESOLVER_HPP
#define HULLWHITESOLVER_HPP

#include <sgpp/pde/application/ParabolicPDESolver.hpp>

#include <sgpp/base/grid/type/LinearGrid.hpp>
#include <sgpp/base/grid/common/BoundingBox.hpp>
#include <sgpp/base/application/ScreenOutput.hpp>

#include <sgpp/globaldef.hpp>
#include <sgpp/base/grid/type/LinearBoundaryGrid.hpp>

#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <cmath>
#include <algorithm>

namespace sgpp {
namespace finance {

/**
 * This class provides a simple-to-use solver of the "multi" dimensional Hull
 * White Equation that uses Sparse Grids.
 *
 * The class's aim is, to hide all complex details of solving the Hull White
 * Equation with Sparse Grids!
 *
 */

class HullWhiteSolver : public sgpp::pde::ParabolicPDESolver {
 private:
  ///  the theta value
  double theta;
  /// the sigma value
  double sigma;
  /// the a value
  double a;
  /// the current time
  // double t;
  /// stores if the stochastic asset data was passed to the solver
  bool bStochasticDataAlloc;
  /// screen object used in this solver
  sgpp::base::ScreenOutput* myScreen;
  /// use coarsening between timesteps in order to reduce gridsize
  bool useCoarsen;
  /// Threshold used to decide if a grid point should be deleted
  double coarsenThreshold;
  /// Threshold used to decide if a grid point should be refined
  double refineThreshold;
  /// adaptive mode during solving Black Scholes Equation: none, coarsen, refine, coarsenNrefine
  std::string adaptSolveMode;
  /// refine mode during solving Black Scholes Equation: classic or maxLevel
  std::string refineMode;
  /// number of points the are coarsened in each coarsening-step
  int numCoarsenPoints;
  /// max. level for refinement during solving
  sgpp::base::GridIndex::level_type refineMaxLevel;
  /// variable to store needed solving iterations

 public:
  /**
   * Std-Constructor of the solver
   */
  HullWhiteSolver();

  /**
   * Std-Destructor of the solver
   */
  virtual ~HullWhiteSolver();

  void constructGrid(sgpp::base::BoundingBox& myBoundingBox, size_t level);

  void setStochasticData(double theta, double sigma, double a);

  void solveImplicitEuler(size_t numTimesteps, double timestepsize, size_t maxCGIterations,
                          double epsilonCG, sgpp::base::DataVector& alpha, bool verbose = false,
                          bool generateAnimation = false);

  void solveExplicitEuler(size_t numTimesteps, double timestepsize, size_t maxCGIterations,
                          double epsilonCG, sgpp::base::DataVector& alpha, bool verbose = false,
                          bool generateAnimation = false);

  void solveCrankNicolson(size_t numTimesteps, double timestepsize, size_t maxCGIterations,
                          double epsilonCG, sgpp::base::DataVector& alpha, size_t NumImEul = 0);

  /**
   * Inits the alpha vector with a payoff function of an European call option
   *
   * @param alpha the coefficient vector of the grid's ansatzfunctions
   * @param strike the option's strike
   * @param payoffType specifies the type of the combined payoff function; std_euro_call or
   * std_euro_put are available
   * @param sigma the sigma value in HullWhite model
   * @param a the value of a in HullWhite model
   * @param t the current time
   * @param T the maturity time
   */
  void initGridWithPayoff(sgpp::base::DataVector& alpha, double strike, std::string payoffType,
                          double sigma, double a, double t, double T);

  /**
   * Inits the screen object
   */
  void initScreen();

  /**
   * returns the algorithmic dimensions (the dimensions in which the Up Down
   * operations (need for space discretization) should be applied)
   *
   * @return the algorithmic dimensions
   */
  std::vector<size_t> getAlgorithmicDimensions();

  /**
   * sets the algorithmic dimensions (the dimensions in which the Up Down
   * operations (need for space discretization) should be applied)
   *
   * @param newAlgoDims std::vector containing the algorithmic dimensions
   */
  void setAlgorithmicDimensions(std::vector<size_t> newAlgoDims);
};
}  // namespace finance
}  // namespace sgpp

#endif /* BLACKSCHOLESSOLVER_HPP */
