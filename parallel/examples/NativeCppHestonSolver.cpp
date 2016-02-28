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
#include <complex>
#include <limits>
#include <vector>

// default number of Implicit Euler steps before starting with Crank Nicolson approach
#define CRNIC_IMEUL_STEPS 3
// default value for epsilon in gridpoints @money
#define DFLT_EPS_AT_MONEY 0.0
// Resolution (num points per dimension) of the created gnuplots
#define PLOT_RESOLUTION 40

double alphaDone;
double vProbe;
double sProbe;
double s2Probe;
double v2Probe;
double refinementThresh;
size_t numGridPoints;

/**
 * Calls the writeHelp method in the HestonSolver Object
 * after creating a screen.
 */
void writeHelp() {
  sgpp::finance::HestonSolver* myHestonSolver = new sgpp::finance::HestonSolver();

  myHestonSolver->initScreen();

  delete myHestonSolver;

  std::stringstream mySStream;

  mySStream << "Some instructions for the use of the Heston Solver:" << std::endl;
  mySStream << "------------------------------------------------------" << std::endl << std::endl;
  mySStream << "Available execution modes are:" << std::endl;
  mySStream << "  solveND             Solves a European Call option" << std::endl;
  mySStream << "                      for N assets on a regular sparse grid" << std::endl
            << std::endl;

  mySStream << "Two files are needed to specify input parameters:" << std::endl;
  mySStream << "-----------------------------------------------------" << std::endl;
  mySStream << "file_Boundaries:  this file contains the grid's bounding box" << std::endl;
  mySStream << "                  for every dimension this file contains a" << std::endl;
  mySStream << "                  tuple with the boundaries. The first dimension is the "
            << std::endl;
  mySStream << "                  stock price for the first asset. The second dimension"
            << std::endl;
  mySStream << "                  is the variance for the first asset. The third dimension"
            << std::endl;
  mySStream << "                  is the stock price for the second asset. The fourth dimension"
            << std::endl;
  mySStream << "                  is the variance for the second asset, and so on." << std::endl;
  mySStream << "Example (two assets (four dimensions)):" << std::endl;
  mySStream << "                  0.0 2.5" << std::endl;
  mySStream << "                  0.01 0.61" << std::endl;
  mySStream << "                  0.0 3.0" << std::endl;
  mySStream << "                  0.02 0.7" << std::endl << std::endl << std::endl;

  mySStream << "file_Stochdata:   this file contains the stochastic Heston" << std::endl;
  mySStream << "                  parameters for the assets." << std::endl;
  mySStream << "                  The i-th line contains the details of the i-th asset."
            << std::endl;
  mySStream << "                  This is the following data (on one line):" << std::endl;
  mySStream << "                  xi_i theta_i kappa_i rho(Si,S1) rho(Si,v1) ... rho(Si,Si) "
               "rho(Si,vi) ... rho(Si,SM) rho(Si,vM)"
            << std::endl;
  mySStream << "                             rho(vi,S1) rho(vi,v1) ... "
               "rho(vi,Si) rho(vi,vi) ... rho(vi,SM) rho(vi,vM)"
            << std::endl;
  mySStream << "Example (2 assets (four dimensions)):" << std::endl;
  mySStream << "                  0.3 0.2 2.0 1.0 -0.5 0.2 0.0" << std::endl;
  mySStream << "                  0.4 0.3 1.5 0.2 0.0 1.0 -0.5" << std::endl;

  mySStream << "Execution modes descriptions:" << std::endl;
  mySStream << "-----------------------------------------------------" << std::endl;
  mySStream << "solveND" << std::endl << "------" << std::endl;
  mySStream << "the following options must be specified:" << std::endl;
  mySStream << "  Coordinates: cart: cartisian coordinates; log: log coords" << std::endl;
  mySStream << "  dim: the number of assets (half the number of dimensions of Sparse Grid)"
            << std::endl;
  mySStream << "  level: number of levels within the Sparse Grid" << std::endl;
  mySStream << "  file_Boundaries: file that contains the bounding box" << std::endl;
  mySStream << "  file_Stochdata: file with the asset's mu, sigma, rho" << std::endl;
  mySStream << "  Strike: the strike" << std::endl;
  mySStream << "  payoff_func: function for n-d payoff: std_euro_{call|put}" << std::endl;
  mySStream << "  r: the riskfree rate" << std::endl;
  mySStream << "  T: time to maturity" << std::endl;
  mySStream << "  dT: timestep size" << std::endl;
  mySStream << "  Solver: the solver to use: CrNic" << std::endl;
  mySStream << "          (for explanations of the options, see end of help!)" << std::endl;
  mySStream << "  CGIterations: Maxmimum number of iterations used in CG mehtod" << std::endl;
  mySStream << "  CGEpsilon: Epsilon used in CG" << std::endl;
  mySStream << std::endl;
  mySStream << "Example:" << std::endl;
  mySStream << "cart 2 5 "
            << "bound.data stoch.data 1.0 std_euro_call "
            << "0.05 "
            << "1.0 "
            << "0.01 ImEul "
            << "400 "
            << "0.000001" << std::endl;
  mySStream << std::endl;
  mySStream << "Remark: This test generates following files (dim<=2):" << std::endl;
  mySStream << "  payoff.gnuplot: the start condition" << std::endl;
  mySStream << "  solvedHeston.gnuplot: the numerical solution" << std::endl;
  mySStream << std::endl << std::endl;

  mySStream << "options for time-stepping:" << std::endl << "------" << std::endl;
  mySStream << "   * CrNic                  Crank-Nicoloson" << std::endl;
  mySStream << std::endl << std::endl;
  mySStream << std::endl << std::endl;
  std::cout << mySStream.str() << std::endl;
}

/**
 * reads the values of mu, sigma and rho of all assets from
 * a file and stores them into three separated DataVectors
 *
 * @param tFile the file that contains the stochastic data
 * @param numAssets the of Assets stored in the file
 * @param xi volatility of the volatility
 * @param theta long-run variance
 * @param kappa mean-reversion rate
 * @param hMatrix correlation matrix
 *
 * @return returns 0 if the file was successfully read, otherwise -1
 */
int readStochasticData(std::string tFile, size_t numAssets, sgpp::base::DataVector& xi,
                       sgpp::base::DataVector& theta, sgpp::base::DataVector& kappa,
                       sgpp::base::DataMatrix& hMatrix) {
  // For the Heston model we need the following stochastic process data
  // xi: each xi value represents the volatility of the volatility (also called volatility of the
  // variance) for a particular asset. For d assets we have a vector of size d here.
  // theta: each theta value represents the long-run variance for a particular asset. For M assets
  // we have a vector of size M here.
  // kappa: each kappa value represents the mean-reversion rate (i.e. the speed which the variance
  // goes back to its long-run value whenever it's not equal to it). For M assets we have a vector
  // of size M here.
  // H: the correlation matrix, of size 2Mx2M, where M is the number of assets

  // Explanation: Correlation in the Heston model:
  // For one asset, we have two processes, dS and dv, the Wiener processes of which are correlated
  // by a factor p.
  // So that means that for one asset we essentially have a 2x2 H matrix. The h11 entry of this
  // matrix is the correlation between the S Wiener process and the
  // S Wiener process (this is usually one). The h12 entry of the matrix is the correlation between
  // the S Wiener process and the v wiener process, and so on for h21 and h22.
  // The correlation matrix is symmetric.
  // For two assets h would be 4x4 and so on.
  // Example of the file for one and two assets:
  //
  // (for one asset)
  // xi theta kappa h11 h12 h21 h22
  //
  // (for two assets)
  // xi1 theta1 kappa1 h11 h12 h13 h14 h21 h22 h23 h24
  // xi2 theta2 kappa2 h31 h32 h33 h34 h41 h42 h43 h44
  //
  // ...and so on for more assets

  std::fstream file;
  double cur_xi;
  double cur_theta;
  double cur_kappa;
  double cur_h;
  size_t hMatrixDim = 2 * numAssets;

  // Count the number of elements in the file
  file.open(tFile.c_str());

  if (!file.is_open()) {
    std::cout << "Error cannot read file: " << tFile << std::endl;
    return -1;
  }

  // There must be a xi, theta and kappa for each asset
  // In addition, there must be a total of (2xnumAssets)^2 hmatrix elements
  // Thus, the element count must be (3xnumAssets)+(2xnumAssets)^2
  size_t t = 0;
  double test;

  do {
    file >> test;
    t++;
  } while (!file.eof());

  file.close();

  if (t < ((hMatrixDim * hMatrixDim) + (3 * numAssets))) {
    std::cout << "Invalid stoch file: " << tFile << " Last Value:" << test << std::endl;
    return -1;
  }

  // Read the file (outermost loop here iterates over each asset line of text)
  file.open(tFile.c_str());

  for (size_t i = 0; i < numAssets; i++) {
    // Read in the xi, theta and kappa values for asset i from the file
    file >> cur_xi;
    file >> cur_theta;
    file >> cur_kappa;

    // Set the xi, theta and kappa values for asset i
    xi.set(i, cur_xi);
    theta.set(i, cur_theta);
    kappa.set(i, cur_kappa);

    // Now we deal with the h matrix values.
    // On each text line (i.e. for each asset), the number of matrix entries is equal
    // to 2*(2*numAssets), i.e. two matrix rows of (2*numAssets) each.
    // We will process each of these rows in a separate loop. The first loop handles the first row
    // (i.e. row index 2*i)
    // and the second loop handles the second row (i.e. row index 2*i + 1).
    for (size_t j = 0; j < hMatrixDim; j++) {
      file >> cur_h;
      hMatrix.set(2 * i, j, cur_h);
    }

    for (size_t j = 0; j < hMatrixDim; j++) {
      file >> cur_h;
      hMatrix.set(2 * i + 1, j, cur_h);
    }
  }

  file.close();

  return 0;
}

/**
 * reads the values of the Bounding Box
 *
 * @param tFile the file that contains the stochastic data
 * @param numDims the number of dimensions
 * @param BoundaryArray Pointer to the Bounding Box array
 *
 * @return returns 0 if the file was successfully read, otherwise -1
 */
int readBoudingBoxData(std::string tFile, size_t numDims,
                       sgpp::base::DimensionBoundary* BoundaryArray) {
  std::fstream file;
  double cur_right;
  double cur_left;

  file.open(tFile.c_str());

  if (!file.is_open()) {
    std::cout << "Error cannot read file: " << tFile << std::endl;
    return -1;
  }

  // Get number of elements in bound file, must be 2*numDims (i.e. an upper and lower bound for each
  // dimension)
  size_t j = 0;
  double test;

  do {
    file >> test;
    j++;
  } while (!file.eof());

  file.close();

  if (j < (numDims * 2)) {
    std::cout << "Invalid boundary file (j=" << j << "): " << tFile << " Last Value:" << test
              << std::endl;
    return -1;
  }

  // Read the boundary data...two elements (left and right) for each dimension.
  file.open(tFile.c_str());

  for (size_t i = 0; i < numDims; i++) {
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
 * Do a Heston solver test with n assets (2N Dimensional Sparse Grid) European call / put option
 *
 * @param numAssets the number of assets
 * @param l the number of levels used in the Sparse Grid
 * @param fileStoch filename of the file that contains the stochastic data (mu, sigma, rho)
 * @param fileBound filename of the file that contains the grid's bounding box
 * @param dStrike the strike of the option
 * @param payoffType method that is used to determine the multidimensional payoff function
 * @param riskfree the riskfree rate of the marketmodel
 * @param timeSt the number of timesteps that are executed during the solving process
 * @param dt the size of delta t in the ODE solver
 * @param CGIt the maximum number of Iterations that are executed by the CG/BiCGStab
 * @param CGeps the epsilon used in the CG/BiCGStab
 * @param Solver specifies the sovler that should be used, ExEul, ImEul and CrNic are the
 * possibilities
 * @param coordsType set the type of coordinates that should be used: cart, log, PAT
 */
void testNUnderlyings(size_t numAssets, int l, std::string fileStoch, std::string fileBound,
                      double dStrike, std::string payoffType, double riskfree, size_t timeSt,
                      double dt, size_t CGIt, double CGeps, std::string Solver,
                      std::string coordsType) {
  size_t numberOfAssets = numAssets;
  size_t pdeDim = numberOfAssets * 2;
  int level = l;
  size_t timesteps = timeSt;
  double stepsize = dt;
  size_t CGiterations = CGIt;
  double CGepsilon = CGeps;
  //  double maxStock = 0.0;

  sgpp::base::DataVector theta(numberOfAssets);
  sgpp::base::DataVector xi(numberOfAssets);
  sgpp::base::DataVector kappa(numberOfAssets);
  sgpp::base::DataMatrix hMatrix(pdeDim, pdeDim);

  double r = riskfree;

  if (readStochasticData(fileStoch, numberOfAssets, xi, theta, kappa, hMatrix) != 0) {
    return;
  }

  // We have boundary data for each dimension in the PDE
  sgpp::base::DimensionBoundary* myBoundaries = new sgpp::base::DimensionBoundary[pdeDim];

  if (readBoudingBoxData(fileBound, pdeDim, myBoundaries) != 0) {
    return;
  }

  sgpp::finance::HestonSolver* myHestonSolver;

  if (coordsType == "log") {
    myHestonSolver = new sgpp::finance::HestonSolver(true);
  } else if (coordsType == "cart") {
    myHestonSolver = new sgpp::finance::HestonSolver(false);
  } else {
    // Write Error
    std::cout << "Unsupported grid transformation!" << std::endl;
    std::cout << std::endl << std::endl;
    writeHelp();
    return;
  }

  sgpp::base::BoundingBox* myBoundingBox = new sgpp::base::BoundingBox(pdeDim, myBoundaries);
  delete[] myBoundaries;

  // init Screen Object
  myHestonSolver->initScreen();

  // Construct a grid
  myHestonSolver->constructGrid(*myBoundingBox, level);

  std::string adaptSolvingMode = "refine";
  std::string refinementMode = "classic";
  //  size_t maxRefineLevel = 10;
  //  double coarsenThreshold = 0.0;
  //  double dRefineThreshold = 0.00001;// See Alex's second thesis
  //  double dRefineThreshold = refinementThresh;

  // Set coarsening dat
  //    myHestonSolver->setEnableCoarseningData(adaptSolvingMode, refinementMode, maxRefineLevel,
  //    -1, coarsenThreshold, dRefineThreshold);

  // init the basis functions' coefficient vector
  sgpp::base::DataVector* alpha = new sgpp::base::DataVector(myHestonSolver->getNumberGridPoints());

  std::cout << "Grid has " << level << " Levels" << std::endl;
  std::cout << "Initial Grid size: " << myHestonSolver->getNumberGridPoints() << std::endl;
  std::cout << "Initial Grid size (inner): " << myHestonSolver->getNumberInnerGridPoints()
            << std::endl
            << std::endl
            << std::endl;

  //  size_t nIterAdaptSteps = 5;
  //  bool useNormalDist = true;
  //  int numRefinePoints = -1;
  std::vector<double> norm_mu;
  std::vector<double> norm_sigma;
  norm_mu.push_back(0.5);
  norm_mu.push_back(5);
  norm_sigma.push_back(0.5);
  norm_sigma.push_back(5);

  // refine the grid to approximate the singularity in the start solution better
  // uncomment if required
  //    if (refinementMode == "classic")
  //    {
  //      for (size_t i = 0 ; i < nIterAdaptSteps; i++)
  //      {
  //        std::cout << "Refining Grid..." << std::endl;
  //        if (useNormalDist == true)
  //        {
  //          myHestonSolver->refineInitialGridSurplusSubDomain(*alpha, numRefinePoints,
  //          dRefineThreshold, norm_mu, norm_sigma);
  //        }
  //        else
  //        {
  //          myHestonSolver->refineInitialGridSurplus(*alpha, numRefinePoints, dRefineThreshold);
  //        }
  //        myHestonSolver->initGridWithPayoff(*alpha, dStrike, payoffType);
  //        std::cout << "Refined Grid size: " << myHestonSolver->getNumberGridPoints() <<
  //        std::endl;
  //        std::cout << "Refined Grid size (inner): " << myHestonSolver->getNumberInnerGridPoints()
  //        << std::endl;
  //      }
  //    }
  //    else
  //    {
  //      std::cout << "An unsupported refinement mode has be chosen!" << std::endl;
  //      std::cout << "Skipping initial grid refinement!" << std::endl;
  //    }

  numGridPoints = myHestonSolver->getNumberGridPoints();

  // Set stochastic data
  myHestonSolver->setStochasticData(theta, kappa, xi, hMatrix, r);

  // Init the grid with on payoff function
  myHestonSolver->initGridWithPayoff(*alpha, dStrike, payoffType);

  // Gridpoints @Money
  std::cout << "Gridpoints @Money: "
            << myHestonSolver->getGridPointsAtMoney(payoffType, dStrike, DFLT_EPS_AT_MONEY)
            << std::endl
            << std::endl
            << std::endl;

  if (numberOfAssets < 2) {
    myHestonSolver->printGrid(*alpha, 100, "payoff.gnuplot");
  }

  if (numberOfAssets < 2) {
    myHestonSolver->printSparseGrid(*alpha, "payoff_surplus.grid.gnuplot", true);
    myHestonSolver->printSparseGrid(*alpha, "payoff_nodal.grid.gnuplot", false);

    if (coordsType == "log") {
      myHestonSolver->printSparseGridExpTransform(*alpha, "payoff_surplus_cart.grid.gnuplot", true);
      myHestonSolver->printSparseGridExpTransform(*alpha, "payoff_nodal_cart.grid.gnuplot", false);
    }
  }

  // Start solving the Heston Equation
  if (Solver == "CrNic") {
    myHestonSolver->solveCrankNicolson(timesteps, stepsize, CGiterations, CGepsilon, *alpha,
                                       CRNIC_IMEUL_STEPS);
  } else {
    std::cout << "!!!! You have chosen an unsupported solver type !!!!" << std::endl;
  }

  if (numberOfAssets < 2) {
    // Print the solved Heston Equation into a gnuplot file
    myHestonSolver->printGrid(*alpha, PLOT_RESOLUTION, "solvedHeston.gnuplot");
  }

  if (numberOfAssets < 2) {
    myHestonSolver->printSparseGrid(*alpha, "solvedHeston_surplus.grid.gnuplot", true);
    myHestonSolver->printSparseGrid(*alpha, "solvedHeston_nodal.grid.gnuplot", false);

    if (coordsType == "log") {
      myHestonSolver->printSparseGridExpTransform(*alpha, "solvedHeston_surplus_cart.grid.gnuplot",
                                                  true);
      myHestonSolver->printSparseGridExpTransform(*alpha, "solvedHeston_nodal_cart.grid.gnuplot",
                                                  false);
    }
  }

  // For the at-the-money price, use the variance value in the middle of the domain
  // Test option @ the money
  std::vector<double> point;

  for (size_t i = 0; i < numAssets; i++) {
    point.push_back(1.0);  // strike
    double middleVol = (myBoundingBox->getBoundary(2 * i + 1).leftBoundary +
                        myBoundingBox->getBoundary(2 * i + 1).rightBoundary) /
                       2.0;
    point.push_back(middleVol);  // middle volatility in the range
  }

  alphaDone = myHestonSolver->evalOption(point, *alpha);
  std::cout << "Optionprice at testpoint (at-the-money, and midpoint variance value):" << alphaDone
            << std::endl
            << std::endl;

  if (numAssets == 1) {
    if (payoffType == "std_euro_call") {
      std::cout << "Analytical solution: (" << point[0] << ", " << point[1] << ") "
                << myHestonSolver->EvaluateHestonPriceExact(
                       point[0], point[1], static_cast<double>(timesteps) * stepsize)
                << std::endl
                << std::endl;
    } else {
      std::cout << "Analytical solution: (" << point[0] << ", " << point[1] << ") "
                << myHestonSolver->EvaluateHestonPriceExactPut(
                       point[0], point[1], static_cast<double>(timesteps) * stepsize)
                << std::endl
                << std::endl;
    }
  }

  delete alpha;
  delete myHestonSolver;
  delete myBoundingBox;
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

  if (option == "solveND") {
    if (argc != 15) {
      writeHelp();
    } else {
      std::string fileStoch;
      std::string fileBound;
      double dStrike;
      std::string ani;
      std::string solver;
      std::string payoff;

      fileStoch.assign(argv[6]);
      fileBound.assign(argv[5]);
      dStrike = atof(argv[7]);
      payoff.assign(argv[8]);
      solver.assign(argv[12]);

      std::string coordsType;
      coordsType.assign(argv[2]);

      testNUnderlyings(atoi(argv[3]), atoi(argv[4]), fileStoch, fileBound, dStrike, payoff,
                       atof(argv[9]), (size_t)(atof(argv[10]) / atof(argv[11])), atof(argv[11]),
                       atoi(argv[13]), atof(argv[14]), solver, coordsType);
    }
  } else {
    writeHelp();
  }

  return 0;
}
