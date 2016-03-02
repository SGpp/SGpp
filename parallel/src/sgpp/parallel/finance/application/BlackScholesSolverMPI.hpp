// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef BLACKSCHOLESSOLVERMPI_HPP
#define BLACKSCHOLESSOLVERMPI_HPP

#include <sgpp/pde/application/ParabolicPDESolver.hpp>

#include <sgpp/base/grid/type/LinearGrid.hpp>
#include <sgpp/base/grid/common/BoundingBox.hpp>
#include <sgpp/solver/ODESolver.hpp>

#include <sgpp/base/grid/generation/functors/SurplusRefinementFunctor.hpp>
#include <sgpp/base/grid/generation/functors/SurplusCoarseningFunctor.hpp>

#include <sgpp/base/tools/StdNormalDistribution.hpp>

#include <sgpp/base/application/ScreenOutput.hpp>
#include <sgpp/base/tools/SGppStopwatch.hpp>

#include <sgpp/globaldef.hpp>
#include <sgpp/base/grid/type/LinearBoundaryGrid.hpp>

#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <cmath>
#include <algorithm>

namespace SGPP {

namespace parallel {

/**
 * This class provides a simple-to-use solver of the multi dimensional Black
 * Scholes Equation that uses Sparse Grids.
 *
 * The class's aim is, to hide all complex details of solving the Black Scholes
 * Equation with Sparse Grids!
 *
 * This version of Black-Scholes Solver provides MPI parallelization!
 *
 * Only European options are support so far!
 *
 */
class BlackScholesSolverMPI : public SGPP::pde::ParabolicPDESolver {
 protected:
  /// vector that contains the assets' weight
  SGPP::base::DataVector* mus;
  /// vector that contains the standard deviations
  SGPP::base::DataVector* sigmas;
  /// Matrix that contains the correlations
  SGPP::base::DataMatrix* rhos;
  /// the riskfree rate
  double r;
  /// stores if the stochastic asset data was passed to the solver
  bool bStochasticDataAlloc;
  /// screen object used in this solver
  SGPP::base::ScreenOutput* myScreen;
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
  /// identifies if the Black Scholes Equation should be solved on a log-transformed grid
  bool useLogTransform;
  /// identifies if the Black Scholes Equation should be solved by using a principal axis
  /// transformation
  bool usePAT;
  /// max. level for refinement during solving
  int refineMaxLevel;
  /// variable to store needed solving iterations
  size_t nNeededIterations;
  /// variable to store the solving time
  double dNeededTime;
  /// variable to store start grid size (Inner Grid)
  size_t staInnerGridSize;
  /// variable to store final grid size (Inner Grid)
  size_t finInnerGridSize;
  /// variable to store average grid size (Inner Grid)
  size_t avgInnerGridSize;
  /// Type of the Option to solve
  std::string tBoundaryType;
  /// Eigenvectors of the co-variance matrix
  SGPP::base::DataMatrix* eigvec_covar;
  /// Eigenvalues of the co-variance matrix
  SGPP::base::DataVector* eigval_covar;
  /// mu hat, tanslation coefficient needed if PAT is used
  SGPP::base::DataVector* mu_hat;
  /// stores the current time until which the option has been solved
  double current_time;
  /// stores the strike of the current option
  double dStrike;
  /// stores the option type of the current option
  std::string payoffType;

  /**
   * returns the option value (payoff value) for an European call option
   *
   * @param assetValue the current asset's value
   * @param strike the strike price of the option
   *
   * @return the call premium
   */
  virtual double get1DEuroCallPayoffValue(double assetValue, double strike);

  /**
   * Inits the alpha vector with a payoff function of an European call option or put option.
   * The grid is initialized based on Cartesian coordinates!
   *
   * @param alpha the coefficient vector of the grid's ansatzfunctions
   * @param strike the option's strike
   * @param payoffType specifies the type of the combined payoff function; std_euro_call or
   * std_euro_put are available
   */
  virtual void initCartesianGridWithPayoff(SGPP::base::DataVector& alpha, double strike,
                                           std::string payoffType);

  /**
   * Inits the alpha vector with a payoff function of an European call option or put option
   * The grid is initialized based on log-transformed coordinates!
   *
   * @param alpha the coefficient vector of the grid's ansatzfunctions
   * @param strike the option's strike
   * @param payoffType specifies the type of the combined payoff function; std_euro_call or
   * std_euro_put are available
   */
  virtual void initLogTransformedGridWithPayoff(SGPP::base::DataVector& alpha, double strike,
                                                std::string payoffType);

  /**
   * Inits the alpha vector with a payoff function of an European call option or put option
   * The grid is initialized based on log-transformed and a principal axis transformation
   * coordinates!
   *
   * @param alpha the coefficient vector of the grid's ansatzfunctions
   * @param strike the option's strike
   * @param payoffType specifies the type of the combined payoff function; std_euro_call or
   * std_euro_put are available
   */
  virtual void initPATTransformedGridWithPayoff(SGPP::base::DataVector& alpha, double strike,
                                                std::string payoffType);

  /**
   * This function calculates for every grid point the value
   * of a normal distribution given by norm_mu and norm_sigma.
   * The result is stored dehierarchized in alpha.
   *
   * This method is overwritten in order to support grids with logarithmic coordinates.
   *
   * @param alpha contains dehierarchized sparse grid coefficients containing the values of the
   * multi dimensional normal distribution after call
   * @param norm_mu the expected values of the normal distribution for every grid dimension
   * @param norm_sigma the standard deviation of the normal distribution for every grid dimension
   */
  virtual void getGridNormalDistribution(SGPP::base::DataVector& alpha,
                                         std::vector<double>& norm_mu,
                                         std::vector<double>& norm_sigma);

 public:
  /**
   * Std-Constructor of the solver
   *
   * @param useLogTransform speciefies if a log transformed formulation should be used for solving
   * BlackScholes Equation
   * @param usePAT speciefies if a principal axis transformation (also enabling a
   * log-transformation) should be used for solving BlackScholes Equation
   */
  explicit BlackScholesSolverMPI(bool useLogTransform = false, bool usePAT = false);

  /**
   * Std-Destructor of the solver
   */
  virtual ~BlackScholesSolverMPI();

  virtual void constructGrid(SGPP::base::BoundingBox& myBoundingBox, int level);

  /**
   * In order to solve the multi dimensional Black Scholes Equation you have to provided
   * some statistical data about the underlying (assets' weight, standard deviation
   * and the correlation between them). This function allows you to set this data.
   *
   * @param mus a DataVector that contains the underlyings' weight
   * @param sigmas a DataVector that contains the underlyings' standard deviations
   * @param rhos a DataMatrix that contains the correlations between the underlyings
   * @param r the riskfree rate used in the market model
   */
  virtual void setStochasticData(SGPP::base::DataVector& mus, SGPP::base::DataVector& sigmas,
                                 SGPP::base::DataMatrix& rhos, double r);

  void solveImplicitEuler(size_t numTimesteps, double timestepsize, size_t maxCGIterations,
                          double epsilonCG, SGPP::base::DataVector& alpha, bool verbose = false,
                          bool generateAnimation = false, size_t numEvalsAnimation = 20);

  void solveExplicitEuler(size_t numTimesteps, double timestepsize, size_t maxCGIterations,
                          double epsilonCG, SGPP::base::DataVector& alpha, bool verbose = false,
                          bool generateAnimation = false, size_t numEvalsAnimation = 20);

  void solveCrankNicolson(size_t numTimesteps, double timestepsize, size_t maxCGIterations,
                          double epsilonCG, SGPP::base::DataVector& alpha, size_t NumImEul = 0);

  void solveX(size_t numTimesteps, double timestepsize, size_t maxCGIterations, double epsilonCG,
              SGPP::base::DataVector& alpha, bool verbose = false, void* myODESolverV = NULL,
              std::string Solver = "ImEul");

  void solveAdamsBashforth(size_t numTimesteps, double timestepsize, size_t maxCGIterations,
                           double epsilonCG, SGPP::base::DataVector& alpha, bool verbose = false);

  void solveSCAC(size_t numTimesteps, double timestepsize, double epsilon, size_t maxCGIterations,
                 double epsilonCG, SGPP::base::DataVector& alpha, bool verbose = false);

  void solveSCH(size_t numTimesteps, double timestepsize, double epsilon, size_t maxCGIterations,
                double epsilonCG, SGPP::base::DataVector& alpha, bool verbose = false);

  void solveSCBDF(size_t numTimesteps, double timestepsize, double epsilon, size_t maxCGIterations,
                  double epsilonCG, SGPP::base::DataVector& alpha, bool verbose = false);

  void solveSCEJ(size_t numTimesteps, double timestepsize, double epsilon, double myAlpha,
                 size_t maxCGIterations, double epsilonCG, SGPP::base::DataVector& alpha,
                 bool verbose = false);

  /**
   * Inits the alpha vector with a payoff function of an European call option or put option
   *
   * @param alpha the coefficient vector of the grid's ansatzfunctions
   * @param strike the option's strike
   * @param payoffType specifies the type of the combined payoff function; std_euro_call or
   * std_euro_put are available
   */
  virtual void initGridWithPayoff(SGPP::base::DataVector& alpha, double strike,
                                  std::string payoffType);

  /**
   * In order to enable some options the payoff type has to be known by all ranks. Call
   * this method on all ranks to allow optimizations. Otherwise your application will deadlock.
   *
   * @param payoffType the payoff type of the current option.
   */
  virtual void setPayoffType(std::string payoffType);

  /**
   * Inits the screen object
   */
  virtual void initScreen();

  /**
   * returns the algorithmic dimensions (the dimensions in which the Up Down
   * operations (need for space discretization) should be applied)
   *
   * @return the algorithmic dimensions
   */
  virtual std::vector<size_t> getAlgorithmicDimensions();

  /**
   * sets the algorithmic dimensions (the dimensions in which the Up Down
   * operations (need for space discretization) should be applied)
   *
   * @param newAlgoDims std::vector containing the algorithmic dimensions
   */
  virtual void setAlgorithmicDimensions(std::vector<size_t> newAlgoDims);

  /**
   *  enables coarsening of grid during solving the Black Scholes
   *  Equation. The coarsening settings have to be specified in order to
   *  enable coarsening.
   *
   *  @param coarsenThreshold Threshold needed to determine if a grid point should be removed
   *  @param refineMode the Mode used for refining the grid: classic or maxLevel
   *  @param refineMaxLevel max. level for refinement during solving
   *  @param adaptSolveMode adaptive mode during solving equation: coarsen, refine, coarsenNrefine
   *  @param numCoarsenPoints number of points coarsened, -1 all coarsenable points are coarsened
   *  @param refineThreshold Threshold needed to determine if a grid point should be refined
   */
  virtual void setEnableCoarseningData(std::string adaptSolveMode, std::string refineMode,
                                       int refineMaxLevel, int numCoarsenPoints,
                                       double coarsenThreshold, double refineThreshold);

  /**
   * Evaluates the current option value
   * at a point given in Cartesian coordinates
   *
   * @param eval_point the point at with the option price should be determined
   * @param alpha the grid's coefficients
   *
   * @return the option price at the given point
   */
  virtual double evalOption(std::vector<double>& eval_point, SGPP::base::DataVector& alpha);

  /**
   * This method transforms a point given
   * in Cartesian coordinates into the coordinates used by the
   * current instance of BlackScholesSolver
   *
   * @param point point given in Cartision coordinates that should be transformed
   */
  virtual void transformPoint(SGPP::base::DataVector& point);

  /**
   * Resets the current solving time.
   *
   * use this in order to get the discounting right when using one
   * instance of multiple option pricings
   */
  virtual void resetSolveTime();

  /**
   * Prints the SGPP::base::Grid Points of the Sparse SGPP::base::Grid either with their node basis
   * value
   * or their hierarchical surplus
   *
   * This function is available for all dimensions
   *
   * @param alpha the coefficients of the grid's ansatzfunctions
   * @param tfilename absoulte path to the file the grid is written into
   * @param bSurplus specifies whether the surplus (true) or the node basis value (false) is written
   */
  void printSparseGridPAT(SGPP::base::DataVector& alpha, std::string tfilename,
                          bool bSurplus) const;

  /**
   * gets the number needed iterations to solve Black Scholes Equation
   *
   * @return number of iterations needed to solve Black Scholes Equation, if called before solving 0
   * is returned
   */
  virtual size_t getNeededIterationsToSolve();

  /**
   * gets needed time in seconds to solve Black Scholes Equation
   *
   * @return needed time in seconds to solve Black Scholes Equation, if called before solving 0 is
   * returned
   */
  virtual double getNeededTimeToSolve();

  /**
   * gets the number of points in start grid
   *
   * @returns the number of points in start grid, if called before constructing grid, 0 is returned
   */
  virtual size_t getStartInnerGridSize();

  /**
   * gets the number of points in final grid
   *
   * @returns the number of points in final grid, if called before solving, 0 is returned
   */
  virtual size_t getFinalInnerGridSize();

  /**
   * gets the number of average gridpoints
   *
   * @returns the number of average gridpoints, if called before solving, 0 is returned
   */
  virtual size_t getAverageInnerGridSize();
};
}  // namespace parallel
}  // namespace SGPP

#endif /* BLACKSCHOLESSOLVER_HPP */
