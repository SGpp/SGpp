// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp_base.hpp>
#include <sgpp_pde.hpp>
#include <sgpp_finance.hpp>
#include <sgpp_parallel.hpp>
#include <sgpp_solver.hpp>
#include <sgpp_datadriven.hpp>
#include <stdlib.h>

#include <iostream>
#include <string>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <vector>

std::string tFileEvalCuboid = "evalCuboid.data";
std::string tFileEvalCuboidValues = "evalCuboidValues.data";

/// default number of Implicit Euler steps before starting with Crank Nicolson approach
#define CRNIC_IMEUL_STEPS 3
/// default value for epsilon in gridpoints at money
#define DFLT_EPS_AT_MONEY 0.0
/// default value for sigma of refinement normal distribution
#define DFLT_SIGMA_REFINE_NORMDIST 0.15

/**
 * reads the values of mu, sigma and rho of all assets from
 * a file and stores them into three separated DataVectors
 *
 * @param tFile the file that contains the stochastic data
 * @param numAssests the of Assets stored in the file
 * @param mu DataVector for the exspected values
 * @param sigma DataVector for standard deviation
 * @param rho DataMatrix for the correlations
 *
 * @return returns 0 if the file was successfully read, otherwise -1
 */
int readStochasticData(std::string tFile, size_t numAssests, sgpp::base::DataVector& mu,
                       sgpp::base::DataVector& sigma, sgpp::base::DataMatrix& rho) {
  std::fstream file;
  double cur_mu;
  double cur_sigma;
  double cur_rho;

  file.open(tFile.c_str());

  if (!file.is_open()) {
    std::cout << "Error cannot read file: " << tFile << std::endl;
    return -1;
  }

  for (size_t i = 0; i < numAssests; i++) {
    file >> cur_mu;
    file >> cur_sigma;
    mu.set(i, cur_mu);
    sigma.set(i, cur_sigma);

    for (size_t j = 0; j < numAssests; j++) {
      file >> cur_rho;
      rho.set(i, j, cur_rho);
    }
  }

  file.close();

  return 0;
}

/**
 * reads the values of the Bounding Box
 *
 * @param tFile the file that contains the stochastic data
 * @param numAssests the of Assets stored in the file
 * @param BoundaryArray Pointer to the Bounding Box array
 *
 * @return returns 0 if the file was successfully read, otherwise -1
 */
int readBoudingBoxData(std::string tFile, size_t numAssests,
                       sgpp::base::DimensionBoundary* BoundaryArray) {
  std::fstream file;
  double cur_right;
  double cur_left;

  file.open(tFile.c_str());

  if (!file.is_open()) {
    std::cout << "Error cannot read file: " << tFile << std::endl;
    return -1;
  }

  for (size_t i = 0; i < numAssests; i++) {
    file >> cur_left;
    file >> cur_right;

    BoundaryArray[i].leftBoundary = cur_left;
    BoundaryArray[i].rightBoundary = cur_right;
    BoundaryArray[i].bDirichletLeft = true;
    BoundaryArray[i].bDirichletRight = true;
  }

  file.close();

  return 0;
}

/**
 * reads a cuboid defined by several points from a file. These points are stored in the
 * cuboid DataMatrix
 *
 * @param cuboid DataMatrix into which the evaluations points are stored
 * @param tFile file that contains the cuboid
 * @param dim the dimensions of cuboid
 */
int readEvalutionCuboid(sgpp::base::DataMatrix& cuboid, std::string tFile, size_t dim) {
  std::fstream file;
  double cur_coord;

  file.open(tFile.c_str());

  if (cuboid.getNcols() != dim) {
    std::cout << "Cuboid-definition file doesn't match: " << tFile << std::endl;
    return -1;
  }

  if (!file.is_open()) {
    std::cout << "Error cannot read file: " << tFile << std::endl;
    return -1;
  }

  // Get number of lines and resize DataMatrix
  size_t i = 0;

  while (!file.eof()) {
    for (size_t d = 0; d < dim; d++) {
      file >> cur_coord;
    }

    i++;
  }

  file.close();
  cuboid.resize(i);

  // Read data from file
  file.open(tFile.c_str());
  i = 0;

  while (!file.eof()) {
    sgpp::base::DataVector line(dim);
    line.setAll(0.0);

    for (size_t d = 0; d < dim; d++) {
      file >> cur_coord;
      line.set(d, cur_coord);
    }

    cuboid.setRow(i, line);
    i++;
  }

  file.close();

  return 0;
}

/**
 * reads function values (here option prices) from a file
 *
 * @param values DataVector into which the values will be stored
 * @param tFile file from which the values are read
 */
int readOptionsValues(sgpp::base::DataVector& values, std::string tFile) {
  std::fstream file;
  double cur_value;

  file.open(tFile.c_str());

  if (!file.is_open()) {
    std::cout << "Error cannot read file: " << tFile << std::endl;
    return -1;
  }

  // Count number of lines
  size_t i = 0;

  while (!file.eof()) {
    file >> cur_value;
    i++;
  }

  values.resize(i);
  file.close();

  // Read data from File
  file.open(tFile.c_str());
  i = 0;

  while (!file.eof()) {
    file >> cur_value;
    values.set(i, cur_value);
    i++;
  }

  file.close();

  return 0;
}

/**
 * calculate the theta value for the Hull White model
 *
 * @param a is the mean reversion rate
 * @param sigma is the volatility
 * @param T is the maturity
 * @param t is the current time
 *
 * @return returns 0 if the file was successfully read, otherwise -1
 */

double calculatetheta(double a, double sigma, double T, double t) {
  double theta = 0;
  return theta = 0.04 * a + pow(sigma, 2.0) * (1 - exp(-2 * a * (T - t))) / (2 * a);
}
/**
 * Combine Hull White solver and Black Scholes solver with European call option
 *
 * @param d dimensions
 * @param l the number of levels used in the Sparse Grid
 * @param sigma
 * @param a
 * @param fileStoch the stochastic data (theta, sigmahw, sigmabs, a)
 * @param fileBound the grid's bounding box - domain boundary(min,max)
 * @param payoffType method that is used to determine the multidimensional payoff function
 * @param timeSt the number of timesteps that are executed during the solving process
 *               (not the number of timesteps in total, but the number of timesteps after which
 *               BS and HW solver are alternatingly combined)
 * @param dt the size of delta t in the ODE solver
 * @param CGIt the maximum number of Iterations that are executed by the CG/BiCGStab
 * @param CGeps the epsilon used in the CG/BiCGStab
 * @param Solver specifies the solver that should be used, ExEul, ImEul and CrNic are the
 * possibilities
 * @param T
 * @param dStrike
 * @param isLogSolve
 */
void testBSHW(size_t d, int l, double sigma, double a, std::string fileStoch, std::string fileBound,
              std::string payoffType, size_t timeSt, double dt, size_t CGIt, double CGeps,
              std::string Solver, double T, double dStrike, bool isLogSolve) {
  size_t dim = d;
  int level = l;
  size_t timesteps_innerCall = timeSt;
  double stepsize_general = dt;
  size_t CGiterations = CGIt;
  double CGepsilon = CGeps;
  sgpp::base::DataVector mu(dim);
  sgpp::base::DataVector sigmabs(dim);
  sgpp::base::DataMatrix rho(dim, dim);

  if (readStochasticData(fileStoch, dim, mu, sigmabs, rho) != 0) {
    return;
  }

  sgpp::base::DimensionBoundary* myBoundaries = new sgpp::base::DimensionBoundary[dim];

  if (readBoudingBoxData(fileBound, dim, myBoundaries) != 0) {
    return;
  }

  // std::cout<<myBoundaries[0].bDirichletLeft << std::endl;
  // std::cout<<myBoundaries[1].bDirichletLeft << std::endl;
  sgpp::finance::BlackScholesHullWhiteSolver* myBSHWSolver;

  if (isLogSolve == true) {
    myBSHWSolver = new sgpp::finance::BlackScholesHullWhiteSolver(true);
  } else {
    myBSHWSolver = new sgpp::finance::BlackScholesHullWhiteSolver(false);
  }

  sgpp::base::BoundingBox* myBoundingBox = new sgpp::base::BoundingBox(dim, myBoundaries);
  // delete[] myBoundaries; // we need them for calculating the evaluation point later!

  // init Screen Object
  myBSHWSolver->initScreen();

  //  // set BS and HW dimension
  //  int dim_BS = 0;  // default: 0
  //  int dim_HW = 1;  // default: 1
  //  // if this method is not called, the default values are set, i.e. BS=0, HW=1
  //  myBSHWSolver->setProcessDimensions(dim_BS, dim_HW);

  // Construct a grid
  myBSHWSolver->constructGrid(*myBoundingBox, level);

  // init the basis functions' coefficient vector
  sgpp::base::DataVector* alpha = new sgpp::base::DataVector(myBSHWSolver->getNumberGridPoints());

  std::cout << "Grid has " << level << " Levels" << std::endl;
  std::cout << "Initial Grid size: " << myBSHWSolver->getNumberGridPoints() << std::endl;
  std::cout << "Initial Grid size (inner): " << myBSHWSolver->getNumberInnerGridPoints()
            << std::endl
            << std::endl
            << std::endl;

  // Init the grid with on payoff function
  myBSHWSolver->initGridWithPayoffBSHW(*alpha, dStrike, payoffType, a, sigma);

  // Gridpoints @Money
  std::cout << "Gridpoints @Money: "
            << myBSHWSolver->getGridPointsAtMoney(payoffType, dStrike, DFLT_EPS_AT_MONEY)
            << std::endl
            << std::endl
            << std::endl;

  std::stringstream level_string;
  level_string << l;

  // Print the payoff function into a gnuplot file
  if (dim < 3) {
    if (dim == 2) {
      sgpp::base::DimensionBoundary* myAreaBoundaries = new sgpp::base::DimensionBoundary[dim];

      for (size_t i = 0; i < 2; i++) {
        myAreaBoundaries[i].leftBoundary = 0.9;
        myAreaBoundaries[i].rightBoundary = 1.1;
      }

      sgpp::base::BoundingBox* myGridArea = new sgpp::base::BoundingBox(dim, myAreaBoundaries);

      myBSHWSolver->printGridDomain(*alpha, 50, *myGridArea,
                                    "payoff_area.level_" + level_string.str() + ".gnuplot");

      delete[] myAreaBoundaries;
      delete myGridArea;
    }

    myBSHWSolver->printGrid(*alpha, 21, ("payoffBSHW.level_" + level_string.str() + ".gnuplot"));
    myBSHWSolver->printSparseGrid(
        *alpha, "payoffBSHW_surplus.grid.level_" + level_string.str() + ".gnuplot", true);
    myBSHWSolver->printSparseGrid(
        *alpha, "payoffBSHW_nodal.grid.level_" + level_string.str() + ".gnuplot", false);

    // myBSHWSolver->printPayoffInterpolationError2D(*alpha,
    // "payoff_interpolation_error.grid.gnuplot", 10000, dStrike);
    if (isLogSolve == true) {
      myBSHWSolver->printSparseGridExpTransform(
          *alpha, "payoffBSHW_surplus_cart.grid.level_" + level_string.str() + ".gnuplot", true);
      myBSHWSolver->printSparseGridExpTransform(
          *alpha, "payoffBSHW_nodal_cart.grid.level_" + level_string.str() + ".gnuplot", false);
    }
  }

  // Set stochastic data
  // myBSHWSolver->setStochasticData(mu, sigmabs, rho, 0.0,theta, sigma, a);

  // Start combining the Black Scholes and Hull White Equation
  if (Solver == "ImEul") {
    double theta = 0;
    int count = 0;
    double dt_outerCall = stepsize_general * static_cast<double>(timesteps_innerCall);
    double t_local = 0.0;

    sgpp::base::SGppStopwatch* myStopwatch = new sgpp::base::SGppStopwatch();
    myStopwatch->start();

    for (int i = 0; i < T / dt_outerCall; i++) {
      theta = calculatetheta(a, sigma, T, t_local);
      myBSHWSolver->setStochasticData(mu, sigmabs, rho, 0.0, theta, sigma, a);
      myBSHWSolver->solveImplicitEuler(timesteps_innerCall, stepsize_general, CGiterations,
                                       CGepsilon, *alpha, false, false, 20);
      count = count + 1;
      t_local += dt_outerCall;

      std::cout << "solved at t = " << t_local << " of T = " << T << std::endl;
    }

    double neededTime = myStopwatch->stop();
    std::cout << "Time to solve in total: " << neededTime << " seconds" << std::endl;
  } else {
    std::cout << "!!!! You have chosen an unsupported solver type !!!!" << std::endl;
  }

  if (dim < 3) {
    // Print the solved Black Scholes Equation into a gnuplot file
    myBSHWSolver->printGrid(*alpha, 21, "solvedBSHW.level_" + level_string.str() + ".gnuplot");
    myBSHWSolver->printSparseGrid(
        *alpha, "solvedBSHW_surplus.grid.level_" + level_string.str() + ".gnuplot", true);
    myBSHWSolver->printSparseGrid(
        *alpha, "solvedBSHW_nodal.grid.level_" + level_string.str() + ".gnuplot", false);
    /*if (isLogSolve == true)
    {
      myBSSolver->printSparseGridExpTransform(*alpha, "solvedBSHW_surplus_cart.grid.gnuplot", true);
      myBSSolver->printSparseGridExpTransform(*alpha, "solvedBSHW_nodal_cart.grid.gnuplot", false);
    }*/
  }

  // Test call @ the money
  std::vector<double> point;

  for (size_t i = 0; i < d; i++) {
    //    if (isLogSolve == true)
    //    {
    //      point.push_back(log(1.0));
    //    }
    //    else
    //    {
    //      point.push_back(1.0);
    //    }

    // take the mean of the given domain
    double point_i = (myBoundaries[i].rightBoundary + myBoundaries[i].leftBoundary) / 2.0;
    point.push_back(point_i);
  }

  delete[] myBoundaries;
  std::cout << "Optionprice at [" << point[0] << ", " << point[1]
            << "] (center): " << myBSHWSolver->evaluatePoint(point, *alpha) << std::endl
            << std::endl;

  // for evaluation at strike
  double PB = 0;
  double PT = 0;

  int timeT = 12;
  int endtime = 32;
  double c = 0.06;
  double rr = 0.05;  // change the risk-free rate here!!

  for (int k = (timeT + 1); k <= endtime; k++) {
    PT = exp(0.04 * (timeT - k) + 0.04 * (1 - exp(-a * (k - timeT))) / a -
             pow(sigma, 2.0) * pow((exp(-a * k) - exp(-a * timeT)), 2.0) *
                 (exp(2 * a * timeT) - 1) / (4 * pow(a, 3.0)) -
             (1 - exp(-a * (k - timeT))) * rr / a);
    PB = PB + c * PT;
  }

  double dS = 0.01;

  for (int j = -10; j <= 10; j++) {
    double S = PB + j * dS;

    std::vector<double> point_strike;
    point_strike.push_back(S);
    point_strike.push_back(rr);

    if (j == 0) {
      std::cout << "Optionprice at [" << point_strike[0] << ", " << point_strike[1]
                << "] (at-the-money!!): " << myBSHWSolver->evaluatePoint(point_strike, *alpha)
                << std::endl;
    } else {
      std::cout << "Optionprice at [" << point_strike[0] << ", " << point_strike[1]
                << "]: " << myBSHWSolver->evaluatePoint(point_strike, *alpha) << std::endl;
    }
  }

  delete alpha;
  delete myBSHWSolver;
  delete myBoundingBox;
}

/**
 * Combine Hull White solver and Black Scholes solver with European call option
 *
 * @param d dimensions
 * @param l the number of levels used in the Sparse Grid
 * @param sigma is the volatility
 * @param a is the mean reversion rate
 * @param fileStoch the stochastic data (theta, sigmahw, sigmabs, a)
 * @param fileBound the grid's bounding box - domain boundary(min,max)
 * @param payoffType method that is used to determine the multidimensional payoff function
 * @param timeSt the number of timesteps that are executed during the solving process
 *               (not the number of timesteps in total, but the number of timesteps after which
 *               BS and HW solver are alternatingly combined)
 * @param dt the size of delta t in the ODE solver
 * @param CGIt the maximum number of Iterations that are executed by the CG/BiCGStab
 * @param CGeps the epsilon used in the CG/BiCGStab
 * @param Solver specifies the sovler that should be used, ExEul, ImEul and CrNic are the
 * possibilities
 * @param T
 * @param dStrike
 * @param isLogSolve set this to true if the log-transformed Black Scholes Equation should be solved
 * @param refinementMode the mode selected for surplus refinement: available: classic, maxLevel
 * @param maxRefineLevel ignored for refinement mode classic, in maxLevel: max. level to which the
 * grid is refined
 * @param numRefinePoints number of points that should be refined in each refine iteration before
 * Black Scholes Equation is solved: -1 try to refine all points steered by threshold
 * @param nIterAdaptSteps number of the iterative Grid Refinement that should be executed
 * @param dRefineThreshold Threshold for a point's surplus for refining this point
 * @param useCoarsen specifies if the grid should be coarsened between timesteps
 * @param adaptSolvingMode specifies which adaptive methods are applied during solving the BS
 * Equation
 * @param coarsenThreshold Threshold to decide, if a grid point should be deleted
 * @param useNormalDist enable local initial refinement based on a normal distribution
 */
void testBSHW_adaptive(size_t d, int l, double sigma, double a, std::string fileStoch,
                       std::string fileBound, std::string payoffType, size_t timeSt, double dt,
                       size_t CGIt, double CGeps, std::string Solver, double T, double dStrike,
                       bool isLogSolve, std::string refinementMode, int numRefinePoints,
                       sgpp::base::GridIndex::level_type maxRefineLevel, size_t nIterAdaptSteps,
                       double dRefineThreshold, bool useCoarsen, std::string adaptSolvingMode,
                       double coarsenThreshold, bool useNormalDist) {
  size_t dim = d;
  int level = l;
  size_t timesteps_innerCall = timeSt;
  double stepsize_general = dt;
  size_t CGiterations = CGIt;
  double CGepsilon = CGeps;
  sgpp::base::DataVector mu(dim);
  sgpp::base::DataVector sigmabs(dim);
  sgpp::base::DataMatrix rho(dim, dim);

  if (readStochasticData(fileStoch, dim, mu, sigmabs, rho) != 0) {
    return;
  }

  sgpp::base::DimensionBoundary* myBoundaries = new sgpp::base::DimensionBoundary[dim];

  if (readBoudingBoxData(fileBound, dim, myBoundaries) != 0) {
    return;
  }

  // std::cout<<myBoundaries[0].bDirichletLeft << std::endl;
  // std::cout<<myBoundaries[1].bDirichletLeft << std::endl;
  sgpp::finance::BlackScholesHullWhiteSolver* myBSHWSolver;

  if (isLogSolve == true) {
    myBSHWSolver = new sgpp::finance::BlackScholesHullWhiteSolver(true);
  } else {
    myBSHWSolver = new sgpp::finance::BlackScholesHullWhiteSolver(false);
  }

  sgpp::base::BoundingBox* myBoundingBox = new sgpp::base::BoundingBox(dim, myBoundaries);
  // delete[] myBoundaries; // we need them for calculating the evaluation point later!

  // init Screen Object
  myBSHWSolver->initScreen();

  // Construct a grid
  myBSHWSolver->constructGrid(*myBoundingBox, level);

  // Enable Coarsening
  if (useCoarsen == true) {
    myBSHWSolver->setEnableCoarseningData(adaptSolvingMode, refinementMode, maxRefineLevel, -1,
                                          coarsenThreshold, dRefineThreshold);
  }

  // init the basis functions' coefficient vector
  sgpp::base::DataVector* alpha = new sgpp::base::DataVector(myBSHWSolver->getNumberGridPoints());

  std::cout << "Grid has " << level << " Levels" << std::endl;
  std::cout << "Initial Grid size: " << myBSHWSolver->getNumberGridPoints() << std::endl;
  std::cout << "Initial Grid size (inner): " << myBSHWSolver->getNumberInnerGridPoints()
            << std::endl
            << std::endl
            << std::endl;

  // Init the grid with on payoff function
  myBSHWSolver->initGridWithPayoffBSHW(*alpha, dStrike, payoffType, a, sigma);

  std::vector<double> norm_mu;
  std::vector<double> norm_sigma;
  double refineSigma = DFLT_SIGMA_REFINE_NORMDIST;

  // estimate refine sigma from evaluation cuboid
  // read Evaluation cuboid
  sgpp::base::DataMatrix EvalCuboid(1, dim);
  int retCuboid = readEvalutionCuboid(EvalCuboid, tFileEvalCuboid, dim);

  // read reference values for evaluation cuboid
  sgpp::base::DataVector EvalCuboidValues(1);
  int retCuboidValues = readOptionsValues(EvalCuboidValues, tFileEvalCuboidValues);

  if (EvalCuboid.getNrows() != EvalCuboidValues.getSize()) {
    retCuboid = 1;
    retCuboidValues = 1;
  }

  if (retCuboid == 0 && retCuboidValues == 0) {
    refineSigma = EvalCuboid.max(0) - EvalCuboid.min(0);
  }

  if (useNormalDist == true) {
    for (size_t i = 0; i < d; i++) {
      norm_mu.push_back(dStrike);
      norm_sigma.push_back(refineSigma);
    }
  }

  // Gridpoints @Money
  std::cout << "Gridpoints @Money: "
            << myBSHWSolver->getGridPointsAtMoney(payoffType, dStrike, DFLT_EPS_AT_MONEY)
            << std::endl
            << std::endl
            << std::endl;

  std::cout << "Initial Grid size: " << myBSHWSolver->getNumberGridPoints() << std::endl;
  std::cout << "Initial Grid size (inner): " << myBSHWSolver->getNumberInnerGridPoints()
            << std::endl
            << std::endl
            << std::endl;

  // refine the grid to approximate the singularity in the start solution better
  if (refinementMode == "classic") {
    for (size_t i = 0; i < nIterAdaptSteps; i++) {
      std::cout << "Refining Grid..." << std::endl;

      if (useNormalDist == true) {
        myBSHWSolver->refineInitialGridSurplusSubDomain(*alpha, numRefinePoints, dRefineThreshold,
                                                        norm_mu, norm_sigma);
      } else {
        myBSHWSolver->refineInitialGridSurplus(*alpha, numRefinePoints, dRefineThreshold);
      }

      myBSHWSolver->initGridWithPayoffBSHW(*alpha, dStrike, payoffType, a, sigma);
      std::cout << "Refined Grid size: " << myBSHWSolver->getNumberGridPoints() << std::endl;
      std::cout << "Refined Grid size (inner): " << myBSHWSolver->getNumberInnerGridPoints()
                << std::endl;
    }

  } else if (refinementMode == "maxLevel") {
    size_t oldGridSize = 0;
    size_t newGridSize = myBSHWSolver->getNumberGridPoints();
    size_t addedGridPoint = 0;
    size_t stepCounter = 0;

    if (nIterAdaptSteps > 0) {
      do {
        oldGridSize = newGridSize;
        std::cout << "Refining Grid..." << std::endl;

        if (useNormalDist == true) {
          myBSHWSolver->refineInitialGridSurplusToMaxLevelSubDomain(
              *alpha, dRefineThreshold, maxRefineLevel, norm_mu, norm_sigma);
        } else {
          myBSHWSolver->refineInitialGridSurplusToMaxLevel(*alpha, dRefineThreshold,
                                                           maxRefineLevel);
        }

        myBSHWSolver->initGridWithPayoffBSHW(*alpha, dStrike, payoffType, a, sigma);
        std::cout << "Refined Grid size: " << myBSHWSolver->getNumberGridPoints() << std::endl;
        std::cout << "Refined Grid size (inner): " << myBSHWSolver->getNumberInnerGridPoints()
                  << std::endl;
        newGridSize = myBSHWSolver->getNumberGridPoints();
        addedGridPoint = newGridSize - oldGridSize;
        stepCounter++;
      } while ((addedGridPoint > 0) && (stepCounter < nIterAdaptSteps));
    }
  } else {
    std::cout << "An unsupported refinement mode has be chosen!" << std::endl;
    std::cout << "Skipping initial grid refinement!" << std::endl;
  }

  std::cout << std::endl << std::endl << std::endl;

  // Print the payoff function into a gnuplot file
  if (dim < 3) {
    if (dim == 2) {
      sgpp::base::DimensionBoundary* myAreaBoundaries = new sgpp::base::DimensionBoundary[dim];

      for (size_t i = 0; i < 2; i++) {
        myAreaBoundaries[i].leftBoundary = 0.9;
        myAreaBoundaries[i].rightBoundary = 1.1;
      }

      sgpp::base::BoundingBox* myGridArea = new sgpp::base::BoundingBox(dim, myAreaBoundaries);

      myBSHWSolver->printGridDomain(*alpha, 50, *myGridArea, "payoff_area.adaptive.gnuplot");

      delete[] myAreaBoundaries;
      delete myGridArea;
    }

    myBSHWSolver->printGrid(*alpha, 21, "payoffBSHW.adaptive.gnuplot");
    myBSHWSolver->printSparseGrid(*alpha, "payoffBSHW_surplus.grid.adaptive.gnuplot", true);
    myBSHWSolver->printSparseGrid(*alpha, "payoffBSHW_nodal.grid.adaptive.gnuplot", false);

    // myBSHWSolver->printPayoffInterpolationError2D(*alpha,
    // "payoff_interpolation_error.grid.gnuplot", 10000, dStrike);
    if (isLogSolve == true) {
      myBSHWSolver->printSparseGridExpTransform(
          *alpha, "payoffBSHW_surplus_cart.grid.adaptive.gnuplot", true);
      myBSHWSolver->printSparseGridExpTransform(
          *alpha, "payoffBSHW_nodal_cart.grid.adaptive.gnuplot", false);
    }
  }

  // Set stochastic data
  // myBSHWSolver->setStochasticData(mu, sigmabs, rho, 0.0,theta, sigma, a);

  // Start combining the Black Scholes and Hull White Equation
  if (Solver == "ImEul") {
    double theta = 0;
    int count = 0;
    double dt_outerCall = stepsize_general * static_cast<double>(timesteps_innerCall);
    double t_local = 0.0;

    sgpp::base::SGppStopwatch* myStopwatch = new sgpp::base::SGppStopwatch();
    myStopwatch->start();

    for (int i = 0; i < T / dt_outerCall; i++) {
      theta = calculatetheta(a, sigma, T, t_local);
      myBSHWSolver->setStochasticData(mu, sigmabs, rho, 0.0, theta, sigma, a);
      myBSHWSolver->solveImplicitEuler(timesteps_innerCall, stepsize_general, CGiterations,
                                       CGepsilon, *alpha, false, false, 20);
      count = count + 1;
      t_local += dt_outerCall;
      std::cout << "solved at t = " << t_local << " of T = " << T << std::endl;
    }

    double neededTime = myStopwatch->stop();
    std::cout << "Time to solve in total: " << neededTime << " seconds" << std::endl;
  } else {
    std::cout << "!!!! You have chosen an unsupported solver type !!!!" << std::endl;
  }

  if (dim < 3) {
    // Print the solved Black Scholes Equation into a gnuplot file
    myBSHWSolver->printGrid(*alpha, 21, "solvedBSHW.gnuplot");
    myBSHWSolver->printSparseGrid(*alpha, "solvedBSHW_surplus.grid.adaptive.gnuplot", true);
    myBSHWSolver->printSparseGrid(*alpha, "solvedBSHW_nodal.grid.adaptive.gnuplot", false);
    /*if (isLogSolve == true)
    {
      myBSSolver->printSparseGridExpTransform(*alpha, "solvedBSHW_surplus_cart.grid.gnuplot", true);
      myBSSolver->printSparseGridExpTransform(*alpha, "solvedBSHW_nodal_cart.grid.gnuplot", false);
    }*/
  }

  // Test call @ the money
  std::vector<double> point;

  for (size_t i = 0; i < d; i++) {
    //    if (isLogSolve == true)
    //    {
    //      point.push_back(log(1.0));
    //    }
    //    else
    //    {
    //      point.push_back(1.0);
    //    }

    // take the mean of the given domain
    double point_i = (myBoundaries[i].rightBoundary + myBoundaries[i].leftBoundary) / 2.0;
    point.push_back(point_i);
  }

  delete[] myBoundaries;
  std::cout << "Optionprice at [" << point[0] << ", " << point[1]
            << "] (center of domain): " << myBSHWSolver->evaluatePoint(point, *alpha) << std::endl
            << std::endl;

  // for evaluation at strike
  double PB = 0;
  double PT = 0;

  int timeT = 12;
  int endtime = 32;
  double c = 0.06;
  double rr = 0.05;  // change the risk-free rate here!!

  for (int k = (timeT + 1); k <= endtime; k++) {
    PT = exp(0.04 * (timeT - k) + 0.04 * (1 - exp(-a * (k - timeT))) / a -
             pow(sigma, 2.0) * pow((exp(-a * k) - exp(-a * timeT)), 2.0) *
                 (exp(2 * a * timeT) - 1) / (4 * pow(a, 3.0)) -
             (1 - exp(-a * (k - timeT))) * rr / a);
    PB = PB + c * PT;
  }

  std::vector<double> point_strike;
  point_strike.push_back(PB);
  point_strike.push_back(rr);

  std::cout << "Optionprice at [" << point_strike[0] << ", " << point_strike[1]
            << "] (at-the-money!!): " << ::std::setprecision(12)
            << myBSHWSolver->evaluatePoint(point_strike, *alpha) << std::endl
            << std::endl;

  delete alpha;
  delete myBSHWSolver;
  delete myBoundingBox;
}

/**
 * Calls the writeHelp method in the BlackScholesHullWhiteSolver Object
 * after creating a screen.
 */
void writeHelp() {
  std::stringstream mySStream;

  mySStream << "Some instructions for the use of combing Hull White and Black Scholes Solver:"
            << std::endl;
  mySStream << "------------------------------------------------------" << std::endl << std::endl;
  mySStream << "Available execution modes are:" << std::endl;
  mySStream << "  CombineBSHW                   Combines Hull-White and Black-Scholes" << std::endl;
  mySStream << "  CombineBSHW_adapt             Combines Hull-White and Black-Scholes,"
            << std::endl;
  mySStream << "                                adaptive version                      "
            << std::endl;
  mySStream << "  CombineBSHW_adapt_subdomain   Combines Hull-White and Black-Scholes,"
            << std::endl;
  mySStream << "                                adaptive version with local restricted adaptivity"
            << std::endl;

  mySStream << "Execution modes descriptions:" << std::endl;
  mySStream << "-----------------------------------------------------" << std::endl;
  mySStream << "CombineBSHW" << std::endl << "------" << std::endl;
  mySStream << "the following options must be specified:" << std::endl;
  mySStream << "  dim: the number of dimensions of Sparse Grid" << std::endl;
  mySStream << "  level: number of levels within the Sparse Grid" << std::endl;
  mySStream << "  value of sigma: sigma value-determine overall level of volatility for hull white"
            << std::endl;
  mySStream << "  value of a: a" << std::endl;
  mySStream << "  file_Stochdata: file with the asset's mu, sigma, rho" << std::endl;
  mySStream << "  file_Boundaries: file that contains the bounding box" << std::endl;
  mySStream << "  payoff_func: function for n-d payoff: std_euro_call, GMIB" << std::endl;
  mySStream << "  simeSt: number of time steps of doing HW and BS when calling "
               "alternatively, default: 1"
            << std::endl;
  mySStream << "  dt: timestep size" << std::endl;
  mySStream << "  CGIterations: Maxmimum number of iterations used in CG method" << std::endl;
  mySStream << "  CGEpsilon: Epsilon used in CG" << std::endl;
  mySStream << "  Solver: the solver to use: ExEul, ImEul or CrNic" << std::endl;
  mySStream << "  T: time to maturity" << std::endl;
  mySStream << "  Strike: the strike" << std::endl;
  mySStream << "  Coordinates: cart: cartisian coordinates; log: log coords, this is only "
               "related to Black-Scholes"
            << std::endl;

  mySStream << std::endl;
  mySStream << "Example:" << std::endl;
  mySStream << "2 5 0.01 0.1 stoch.data  bound.data "
            << " std_euro_call "
            << "1.0 "
            << "0.01 "
            << "400 "
            << "0.000001 "
            << "ImEul "
            << "1.0 "
            << "0.3 "
            << "cart " << std::endl;
  mySStream << std::endl;

  mySStream << "CombineBSHW_adapt, CombineBSHW_adapt_subdomain" << std::endl
            << "------" << std::endl;
  mySStream << "the following options must be specified:" << std::endl;
  mySStream << "  dim: the number of dimensions of Sparse Grid" << std::endl;
  mySStream << "  level: number of levels within the Sparse Grid" << std::endl;
  mySStream << "  value of sigma: sigma value-determine overall level of volatility for hull white"
            << std::endl;
  mySStream << "  value of a: a" << std::endl;
  mySStream << "  file_Stochdata: file with the asset's mu, sigma, rho" << std::endl;
  mySStream << "  file_Boundaries: file that contains the bounding box" << std::endl;
  mySStream << "  payoff_func: function for n-d payoff: std_euro_call, GMIB" << std::endl;
  mySStream << "  simeSt: number of time steps of doing HW and BS when calling "
               "allternatively, default: 1"
            << std::endl;
  mySStream << "  dt: timestep size" << std::endl;
  mySStream << "  CGIterations: Maxmimum number of iterations used in CG method" << std::endl;
  mySStream << "  CGEpsilon: Epsilon used in CG" << std::endl;
  mySStream << "  Solver: the solver to use: ExEul, ImEul or CrNic" << std::endl;
  mySStream << "  T: time to maturity" << std::endl;
  mySStream << "  Strike: the strike" << std::endl;
  mySStream << "  Coordinates: cart: cartisian coordinates; log: log coords, this is only "
               "related to Black-Scholes"
            << std::endl;
  mySStream << "  RefinementMode: classic or maxLevel" << std::endl;
  mySStream << "  MaxRefinement Level: Max. Level for refinement" << std::endl;
  mySStream << "  numAdaptRefinement: Number of adaptive refinements at the beginning" << std::endl;
  mySStream << "  refinementThreshold: Threshold of point's surplus to refine point" << std::endl;
  mySStream << "  adapt-mode during solving: none, coarsen, refine, coarsenNrefine" << std::endl;
  mySStream << "  Coarsening Threshold: Threshold of point's surplus to remove point" << std::endl;

  mySStream << std::endl;
  mySStream << "Example:" << std::endl;
  mySStream << "2 5 0.01 0.1 stoch.data  bound.data "
            << " std_euro_call "
            << "1.0 "
            << "0.01 "
            << "400 "
            << "0.000001 "
            << "ImEul "
            << "1.0 "
            << "0.3 "
            << "cart "
            << "maxLevel 10 5 1e-10 coarsen 1e-6" << std::endl;
  mySStream << std::endl;

  mySStream << std::endl;

  mySStream << std::endl << std::endl;
  std::cout << mySStream.str() << std::endl;
}

/**
 * main routine of the application, do some first cli
 * correction test and branches to right solver configuration
 *
 * @param argc contains the number of cli arguments
 * @param argv contains the cli arguments as C-Strings
 */
int main(int argc, char* argv[]) {
  std::string option;

  if (argc == 1) {
    writeHelp();

    return 0;
  }

  option.assign(argv[1]);

  if (option == "CombineBSHW") {
    if (argc != 17) {
      writeHelp();
    } else {
      std::string solver;
      std::string payoff;
      double sigma;
      double a;
      double dStrike;
      std::string fileStoch;
      std::string fileBound;
      sigma = atof(argv[4]);
      a = atof(argv[5]);
      fileStoch.assign(argv[6]);
      fileBound.assign(argv[7]);
      payoff.assign(argv[8]);
      solver.assign(argv[13]);
      dStrike = atof(argv[15]);
      std::string coordsType;
      bool coords = false;
      coordsType.assign(argv[16]);

      if (coordsType == "cart") {
        coords = false;
      } else if (coordsType == "log") {
        coords = true;
      } else {
        std::cout << "Unsupported coordinate option! cart or log are supported!" << std::endl;
        std::cout << std::endl << std::endl;
        writeHelp();
      }

      testBSHW(atoi(argv[2]), atoi(argv[3]), sigma, a, fileStoch, fileBound, payoff, atoi(argv[9]),
               atof(argv[10]), atoi(argv[11]), atof(argv[12]), solver, atof(argv[14]), dStrike,
               coords);
    }

  } else if (option == "CombineBSHW_adapt" || option == "CombineBSHW_adapt_subdomain") {
    if (argc != 23) {
      writeHelp();
    } else {
      bool isNormalDist = false;

      if (option == "CombineBSHW_adapt_subdomain") {
        isNormalDist = true;
      }

      std::string solver;
      std::string payoff;
      double sigma;
      double a;
      double dStrike;
      std::string fileStoch;
      std::string fileBound;
      sigma = atof(argv[4]);
      a = atof(argv[5]);
      fileStoch.assign(argv[6]);
      fileBound.assign(argv[7]);
      payoff.assign(argv[8]);
      solver.assign(argv[13]);
      dStrike = atof(argv[15]);
      std::string coordsType;
      bool coords = false;
      coordsType.assign(argv[16]);
      std::string refinementMode;
      std::string adaptSolveMode;

      if (coordsType == "cart") {
        coords = false;
      } else if (coordsType == "log") {
        coords = true;
      } else {
        std::cout << "Unsupported coordinate option! cart or log are supported!" << std::endl;
        std::cout << std::endl << std::endl;
        writeHelp();
      }

      refinementMode.assign(argv[17]);
      adaptSolveMode.assign(argv[21]);

      if (refinementMode != "maxLevel" && refinementMode != "classic") {
        std::cout << "Unsupported refinement type! classic or maxLevel are supported!" << std::endl;
        std::cout << std::endl << std::endl;
        writeHelp();
        return 0;
      }

      bool useAdaptSolve = false;

      if (adaptSolveMode == "coarsen" || adaptSolveMode == "refine" ||
          adaptSolveMode == "coarsenNrefine") {
        useAdaptSolve = true;
      } else if (adaptSolveMode == "none") {
        useAdaptSolve = false;
      } else {
        std::cout << "Unsupported adapt solve mode! none, coarsen, refine or coarsenNrefine are "
                     "supported!"
                  << std::endl;
        std::cout << std::endl << std::endl;
        writeHelp();
        return 0;
      }

      testBSHW_adaptive(atoi(argv[2]), atoi(argv[3]), sigma, a, fileStoch, fileBound, payoff,
                        atoi(argv[9]), atof(argv[10]), atoi(argv[11]), atof(argv[12]), solver,
                        atof(argv[14]), dStrike, coords, refinementMode, -1, atoi(argv[18]),
                        atoi(argv[19]), atof(argv[20]), useAdaptSolve, adaptSolveMode,
                        atof(argv[22]), isNormalDist);
    }

  } else {
    writeHelp();
  }

  return 0;
}
