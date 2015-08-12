// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <mpi.h>
#include "sgpp_base.hpp"
#include "sgpp_pde.hpp"
#include "sgpp_finance.hpp"
#include "sgpp_parallel.hpp"
#include "sgpp_solver.hpp"
#include "sgpp_datadriven.hpp"
#include "sgpp_mpi.hpp"

#include <iostream>
#include <string>
#include <stdlib.h>
#include <fstream>
#include <iomanip>
#include <cmath>

std::string tFileEvalCuboid = "evalCuboid.MPI.data";
std::string tFileEvalCuboidValues = "evalCuboidValues.MPI.data";

/// default number of Implicit Euler steps before starting with Crank Nicolson approach
#define CRNIC_IMEUL_STEPS 3
/// default value for epsilon in gridpoints at money
#define DFLT_EPS_AT_MONEY 0.0
/// default value for sigma of refinement normal distribution
#define DFLT_SIGMA_REFINE_NORMDIST 0.15

/**
 * Calls the writeHelp method in the BlackScholesSolver Object
 * after creating a screen.
 */
void writeHelp() {
  SGPP::parallel::BlackScholesSolverMPI* myBSSolver = new SGPP::parallel::BlackScholesSolverMPI();

  myBSSolver->initScreen();

  delete myBSSolver;

  std::stringstream mySStream;

  mySStream << "Some instructions for the use of Black Scholes Solver:" << std::endl;
  mySStream << "------------------------------------------------------" << std::endl << std::endl;
  mySStream << "Available execution modes are:" << std::endl;
  mySStream << "  solveNDanalyze      same as solveND, but the option is" << std::endl;
  mySStream << "                      solved for several regular grids with" << std::endl;
  mySStream << "                      different numbers of levels" << std::endl << std::endl;
  mySStream << "  solveNDadaptSurplus Solves an European Call/Up option" << std::endl;
  mySStream << "                      on a refined grid based on" << std::endl;
  mySStream << "                      the hierarchical surplus" << std::endl << std::endl;
  mySStream << "  solveNDadaptSurplusSubDomain   Same as above but" << std::endl;
  mySStream << "						a normal distribution is used" << std::endl;
  mySStream << "						to do refinement just near the strike!" << std::endl << std::endl;

  mySStream << "Several files are needed to specify inputs:" << std::endl;
  mySStream << "-----------------------------------------------------" << std::endl;
  mySStream << "file_Boundaries:  this file contains the grid's bounding box" << std::endl;
  mySStream << "                  for every dimension this file contains a" << std::endl;
  mySStream << "                  tuple with the boundaries." << std::endl;
  mySStream << "Example (3 dimensions):" << std::endl;
  mySStream << "                  0.0 2.5" << std::endl;
  mySStream << "                  0.0 2.5" << std::endl;
  mySStream << "                  0.0 2.5" << std::endl << std::endl << std::endl;

  mySStream << "file_Stochdata:   this file contains the option's asset's" << std::endl;
  mySStream << "                  expected values, standard deviations and" << std::endl;
  mySStream << "                  correlations. The i-th line contains" << std::endl;
  mySStream << "                  followning data:" << std::endl;
  mySStream << "                  mu_i sigma_i rho_i_0 ... rhi_i_d" << std::endl;
  mySStream << "Example (3 dimensions):" << std::endl;
  mySStream << "                  0.05 0.4 1.0 0.1 0.2" << std::endl;
  mySStream << "                  0.05 0.5 0.1 1.0 0.3" << std::endl;
  mySStream << "                  0.05 0.6 0.2 0.3 1.0" << std::endl << std::endl << std::endl;

  mySStream << "file_analyze:     this file contains the options for" << std::endl;
  mySStream << "                  the analyzing runs. This file contains" << std::endl;
  mySStream << "                  two parts: The first lines is the " << std::endl;
  mySStream << "                  evaluation cuboid as bounding box. " << std::endl;
  mySStream << "                  The second one is the number of points" << std::endl;
  mySStream << "                  in every dimension in the evaluation" << std::endl;
  mySStream << "                  cuboid." << std::endl;
  mySStream << "Example (3 dimensions):" << std::endl;
  mySStream << "                  0.0 1.0" << std::endl;
  mySStream << "                  0.0 1.0" << std::endl;
  mySStream << "                  0.0 1.0" << std::endl;
  mySStream << "                  20" << std::endl << std::endl << std::endl;

  mySStream << "Execution modes descriptions:" << std::endl;
  mySStream << "-----------------------------------------------------" << std::endl;
  mySStream << "solveNDanalyze" << std::endl << "------" << std::endl;
  mySStream << "the following options must be specified:" << std::endl;
  mySStream << "	Coordinates: cart: cartisian coordinates; log: log coords; PAT principal axis transform" << std::endl;
  mySStream << "	dim: the number of dimensions of Sparse Grid" << std::endl;
  mySStream << "	level_start: number of levels within the Sparse Grid (start)" << std::endl;
  mySStream << "	level_end: number of levels within the Sparse Grid (end)" << std::endl;
  mySStream << "	file_Boundaries: file that contains the bounding box" << std::endl;
  mySStream << "	file_Stochdata: file with the asset's mu, sigma, rho" << std::endl;
  mySStream << "	Strikes: the strike" << std::endl;
  mySStream << "	payoff_func: function for n-d payoff: std_euro_{call|put}" << std::endl;
  mySStream << "	r: the riskfree rate" << std::endl;
  mySStream << "	T: time to maturity" << std::endl;
  mySStream << "	dT: timestep size" << std::endl;
  mySStream << "	Solver: the solver to use: ExEul, ImEul or CrNic" << std::endl;
  mySStream << "	CGIterations: Maxmimum number of iterations used in CG mehtod" << std::endl;
  mySStream << "	CGEpsilon: Epsilon used in CG" << std::endl;
  mySStream << "	file_analyze: file containing the analyzing options" << std::endl;
  mySStream << std::endl;
  mySStream << "Example:" << std::endl;
  mySStream << "cart 3 2 5 " << "bound.data stoch.data 1.0 std_euro_call " << "0.05 " << "1.0 " << "0.01 ImEul " << "400 " << "0.000001 anal.data" << std::endl;
  mySStream << std::endl;
  mySStream << "Remark: This test generates following files (dim<=2):" << std::endl;
  mySStream << "	payoff.gnuplot: the start condition" << std::endl;
  mySStream << "	solvedBS.gnuplot: the numerical solution" << std::endl;
  mySStream << "For all cases following files are generated:" << std::endl;
  mySStream << "	EvalCuboidPoints.data: containing the evaluation" << std::endl;
  mySStream << "		cuboid" << std::endl;
  mySStream << "	HighLevelOptionValue.data: containing the option's" << std::endl;
  mySStream << "		for the highest leveled grid." << std::endl;
  mySStream << std::endl << std::endl;

  mySStream << "solveNDadaptSurplus/solveNDadaptSurplusSubDomain" << std::endl << "------" << std::endl;
  mySStream << "the following options must be specified:" << std::endl;
  mySStream << "	Coordinates: cart: cartisian coordinates; log: log coords" << std::endl;
  mySStream << "	dim: the number of dimensions of Sparse Grid" << std::endl;
  mySStream << "	level: number of levels within the Sparse Grid" << std::endl;
  mySStream << "	file_Boundaries: file that contains the bounding box" << std::endl;
  mySStream << "	file_Stochdata: file with the asset's mu, sigma, rho" << std::endl;
  mySStream << "	Strike: the strike" << std::endl;
  mySStream << "	payoff_func: function for n-d payoff: std_euro_{call|put}" << std::endl;
  mySStream << "	r: the riskfree rate" << std::endl;
  mySStream << "	T: time to maturity" << std::endl;
  mySStream << "	dT: timestep size" << std::endl;
  mySStream << "	Solver: the solver to use: ExEul, ImEul or CrNic" << std::endl;
  mySStream << "	CGIterations: Maxmimum number of iterations used in CG mehtod" << std::endl;
  mySStream << "	CGEpsilon: Epsilon used in CG" << std::endl;
  mySStream << "	RefinementMode: classic or maxLevel" << std::endl;
  mySStream << "	MaxRefinement Level: Max. Level for refinement" << std::endl;
  mySStream << "	numAdaptRefinement: Number of adaptive refinements at the beginning" << std::endl;
  mySStream << "	refinementThreshold: Threshold of point's surplus to refine point" << std::endl;
  mySStream << "	adapt-mode during solving: none, coarsen, refine, coarsenNrefine" << std::endl;
  mySStream << "	Coarsening Threshold: Threshold of point's surplus to remove point" << std::endl;
  mySStream << std::endl;
  mySStream << "Example:" << std::endl;
  mySStream << "cart 3 5 " << "bound.data stoch.data 1.0 std_euro_call " << "0.05 " << "1.0 " << "0.01 ImEul " << "400 " << "0.000001 classic 0 5 1e-10 coarsen 1e-6" << std::endl;
  mySStream << std::endl;
  mySStream << "Remark: This test generates following files (dim<=2):" << std::endl;
  mySStream << "	payoff.gnuplot: the start condition" << std::endl;
  mySStream << "	solvedBS.gnuplot: the numerical solution" << std::endl;
  mySStream << std::endl << std::endl;

  mySStream << std::endl << std::endl;
  std::cout << mySStream.str() << std::endl;
}

/**
 * reads the values of mu, sigma and rho of all assets from
 * a file and stores them into three separated SGPP::base::DataVectors
 *
 * @param tFile the file that contains the stochastic data
 * @param numAssets the of Assets stored in the file
 * @param mu SGPP::base::DataVector for the exspected values
 * @param sigma SGPP::base::DataVector for standard deviation
 * @param rho SGPP::base::DataMatrix for the correlations
 *
 * @return returns 0 if the file was successfully read, otherwise -1
 */
int readStochasticData(std::string tFile, size_t numAssets, SGPP::base::DataVector& mu, SGPP::base::DataVector& sigma, SGPP::base::DataMatrix& rho) {
  std::fstream file;
  double cur_mu;
  double cur_sigma;
  double cur_rho;

  file.open(tFile.c_str());

  if (!file.is_open()) {
    std::cout << "Error cannot read file: " << tFile << std::endl;
    return -1;
  }

  // Get number of elements in stoch file, must be numAssests*(numAssests+2)
  size_t t = 0;
  double test;

  do {
    file >> test;
    t++;
  } while (!file.eof());

  file.close();

  if (t < ((numAssets * numAssets) + (2 * numAssets))) {
    std::cout << "Invalid stoch file: " << tFile << " Last Value:" << test << std::endl;
    return -1;
  }

  file.open(tFile.c_str());

  for (size_t i = 0; i < numAssets; i++) {
    file >> cur_mu;
    file >> cur_sigma;
    mu.set(i, cur_mu);
    sigma.set(i, cur_sigma);

    for (size_t j = 0; j < numAssets; j++) {
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
 * @param numAssets the of Assets stored in the file
 * @param BoundaryArray Pointer to the Bounding Box array
 *
 * @return returns 0 if the file was successfully read, otherwise -1
 */
int readBoudingBoxData(std::string tFile, size_t numAssets, SGPP::base::DimensionBoundary* BoundaryArray) {
  std::fstream file;
  double cur_right;
  double cur_left;

  file.open(tFile.c_str());

  if (!file.is_open()) {
    std::cout << "Error cannot read file: " << tFile << std::endl;
    return -1;
  }

  // Get number of elements in bound file, must be 2*numAssests
  size_t j = 0;
  double test;

  do {
    file >> test;
    j++;
  } while (!file.eof());

  file.close();

  if (j < (numAssets * 2)) {
    std::cout << "Invalid boundary file: " << tFile << " Last Value:" << test << std::endl;
    return -1;
  }

  file.open(tFile.c_str());

  for (size_t i = 0; i < numAssets; i++) {
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
 * reads the analyze configuration from a file
 *
 * @param tFile the file that contains the analyze data
 * @param numAssets the number of assets
 * @param BoundaryArray Pointer to the Bounding Box array
 * @param points variable to store the number of points in every dimension
 *
 * @return returns 0 if the file was successfully read, otherwise -1
 */
int readAnalyzeData(std::string tFile, size_t numAssets, SGPP::base::DimensionBoundary* BoundaryArray, size_t& points) {
  std::fstream file;
  double cur_right;
  double cur_left;

  file.open(tFile.c_str());

  if (!file.is_open()) {
    std::cout << "Error cannot read file: " << tFile << std::endl;
    return -1;
  }

  // Get number of elements in analyze file, must be 2*numAssests+1
  size_t j = 0;
  double test;

  do {
    file >> test;
    j++;
  } while (!file.eof());

  file.close();

  if (j < ((numAssets * 2) + 1)) {
    std::cout << "Invalid analyze file: " << tFile << " Last Value:" << test << std::endl;
    return -1;
  }

  file.open(tFile.c_str());


  for (size_t i = 0; i < numAssets; i++) {
    file >> cur_left;
    file >> cur_right;

    BoundaryArray[i].leftBoundary = cur_left;
    BoundaryArray[i].rightBoundary = cur_right;
    BoundaryArray[i].bDirichletLeft = true;
    BoundaryArray[i].bDirichletRight = true;
  }

  file >> points;

  file.close();

  return 0;
}

/**
 * reads a cuboid defined by several points from a file. These points are stored in the
 * cuboid SGPP::base::DataMatrix
 *
 * @param cuboid SGPP::base::DataMatrix into which the evaluations points are stored
 * @param tFile file that contains the cuboid
 * @param dim the dimensions of cuboid
 */
int readEvalutionCuboid(SGPP::base::DataMatrix& cuboid, std::string tFile, size_t dim) {
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

  // Get number of lines and resize SGPP::base::DataMatrix
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
    SGPP::base::DataVector line(dim);
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
 * @param values SGPP::base::DataVector into which the values will be stored
 * @param tFile file from which the values are read
 */
int readOptionsValues(SGPP::base::DataVector& values, std::string tFile) {
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
 * Writes a SGPP::base::DataMatrix into a file
 *
 * @param data the SGPP::base::DataMatrix that should be written into a file
 * @param tFile the file into which the data is written
 *
 * @return error code
 */
int writeDataMatrix(SGPP::base::DataMatrix& data, std::string tFile) {
  std::ofstream file;
  file.open(tFile.c_str());

  if (!file.is_open()) {
    std::cout << "Error cannot write file: " << tFile << std::endl;
    return -1;
  }

  for (size_t i = 0; i < data.getNrows(); i++) {
    for (size_t j = 0; j < data.getNcols(); j++) {
      file << std::scientific << std::setprecision( 16 ) << data.get(i, j) << " ";
    }

    file << std::endl;
  }

  file.close();

  return 0;
}


/**
 * Writes a SGPP::base::DataVector into a file
 *
 * @param data the SGPP::base::DataVector that should be written into a file
 * @param tFile the file into which the data is written
 *
 * @return error code
 */
int writeDataVector(SGPP::base::DataVector& data, std::string tFile) {
  std::ofstream file;
  file.open(tFile.c_str());

  if (!file.is_open()) {
    std::cout << "Error cannot write file: " << tFile << std::endl;
    return -1;
  }

  for (size_t i = 0; i < data.getSize(); i++) {

    file << std::scientific << std::setprecision( 16 ) << data.get(i) << " " << std::endl;
  }

  file.close();

  return 0;
}

/**
 * Do a Black Scholes solver test with n assets (ND Sparse Grid) European call option
 *
 * @param d the number of dimensions used in the Sparse Grid
 * @param start_l the number of levels used in the Sparse Grid (first test)
 * @param end_l the number of level used in the Sparse Grid (last test)
 * @param fileStoch filename of the file that contains the stochastic data (mu, sigma, rho)
 * @param fileBound filename of the file that contains the grid's bounding box
 * @param dStrike the strike of the option
 * @param payoffType method that is used to determine the multidimensional payoff function
 * @param riskfree the riskfree rate of the marketmodel
 * @param timeSt the number of timesteps that are executed during the solving process
 * @param dt the size of delta t in the ODE solver
 * @param CGIt the maximum number of Iterations that are executed by the CG/BiCGStab
 * @param CGeps the epsilon used in the CG/BiCGStab
 * @param Solver specifies the sovler that should be used, ExEul, ImEul and CrNic are the possibilities
 * @param fileAnalyze filename of the file that contains the analyze data
 * @param coordsType set the type of coordinates that should be used: cart, log, PAT
 */
void testNUnderlyingsAnalyze(size_t d, size_t start_l, size_t end_l, std::string fileStoch, std::string fileBound, double dStrike, std::string payoffType,
                             double riskfree, size_t timeSt, double dt, size_t CGIt, double CGeps, std::string Solver, std::string fileAnalyze, std::string coordsType) {
  size_t dim = d;
  size_t timesteps = timeSt;
  double stepsize = dt;
  size_t CGiterations = CGIt;
  double CGepsilon = CGeps;

  SGPP::base::DataVector mu(dim);
  SGPP::base::DataVector sigma(dim);
  SGPP::base::DataMatrix rho(dim, dim);

  SGPP::base::DataMatrix EvalPoints(1, d);

  double r = riskfree;

  SGPP::base::DataVector* alpha;

  std::vector<SGPP::base::DataVector> results;

  // read process configuration
  if (readStochasticData(fileStoch, dim, mu, sigma, rho) != 0) {
    return;
  }

  // create Black Scholes Solver Object
  SGPP::parallel::BlackScholesSolverMPI* myBSSolver = NULL;

  if (coordsType == "log") {
    myBSSolver = new SGPP::parallel::BlackScholesSolverMPI(true);
  } else if (coordsType == "cart") {
    myBSSolver = new SGPP::parallel::BlackScholesSolverMPI(false);
  } else if (coordsType == "PAT") {
    myBSSolver = new SGPP::parallel::BlackScholesSolverMPI(true, true);
  } else {
    // Write Error
    std::cout << "Unsupported grid transformation!" << std::endl;
    std::cout << std::endl << std::endl;

    if (SGPP::parallel::myGlobalMPIComm->getMyRank() == 0) {
      writeHelp();
    }

    SGPP::parallel::myGlobalMPIComm->Abort();
  }

  // Screen initialization only on rank 0
  if (SGPP::parallel::myGlobalMPIComm->getMyRank() == 0) {
    // init Screen Object
    myBSSolver->initScreen();
  }

  for (size_t i = start_l; i <= end_l; i++) {
    int level = (int)i;

    // Reset Solve Time
    myBSSolver->resetSolveTime();

    if (SGPP::parallel::myGlobalMPIComm->getMyRank() == 0) {
      SGPP::base::DimensionBoundary* myBoundaries = new SGPP::base::DimensionBoundary[dim];

      if (readBoudingBoxData(fileBound, dim, myBoundaries) != 0) {
        return;
      }

      SGPP::base::BoundingBox* myBoundingBox = new SGPP::base::BoundingBox(dim, myBoundaries);
      delete[] myBoundaries;

      // Construct a grid
      myBSSolver->constructGrid(*myBoundingBox, level);
      delete myBoundingBox;

      // in first iteration -> calculate the evaluation points
      if (i == start_l) {
        size_t points = 0;
        SGPP::base::DimensionBoundary* myEvalBoundaries = new SGPP::base::DimensionBoundary[dim];

        if (readAnalyzeData(fileAnalyze, dim, myEvalBoundaries, points) != 0) {
          return;
        }

        SGPP::base::BoundingBox* myEvalBoundingBox = new SGPP::base::BoundingBox(dim, myEvalBoundaries);
        delete[] myEvalBoundaries;

        SGPP::base::EvalCuboidGenerator* myEvalCuboidGen = new SGPP::base::EvalCuboidGenerator();

        myEvalCuboidGen->getEvaluationCuboid(EvalPoints, *myEvalBoundingBox, points);
        writeDataMatrix(EvalPoints, tFileEvalCuboid);

        std::cout << "=====================================================================" << std::endl;
        std::cout << "=====================================================================" << std::endl << std::endl;
        std::cout << "Calculating norms of relative errors to a grid" << std::endl;
        std::cout << "with the bounding box:" << std::endl;

        for (size_t j = 0; j < d; j++) {
          std::cout << myEvalBoundingBox->getBoundary(j).leftBoundary << " " << myEvalBoundingBox->getBoundary(j).rightBoundary << std::endl;
        }

        std::cout << "=====================================================================" << std::endl;
        std::cout << "=====================================================================" << std::endl << std::endl << std::endl;

        delete myEvalBoundingBox;
        delete myEvalCuboidGen;
      }

      // init the basis functions' coefficient vector
      alpha = new SGPP::base::DataVector(myBSSolver->getNumberGridPoints());

      std::cout << "Grid has " << level << " Levels" << std::endl;
      std::cout << "Initial Grid size: " << myBSSolver->getNumberGridPoints() << std::endl;
      std::cout << "Initial Grid size (inner): " << myBSSolver->getNumberInnerGridPoints() << std::endl << std::endl << std::endl;

      // Set stochastic data
      myBSSolver->setStochasticData(mu, sigma, rho, r);

      // Init the grid with on payoff function
      myBSSolver->initGridWithPayoff(*alpha, dStrike, payoffType);

      // Print the payoff function into a gnuplot file
      if (dim < 3) {
        myBSSolver->printGrid(*alpha, 20, "payoff.MPI.gnuplot");
      }

      if (dim < 4) {
        myBSSolver->printSparseGrid(*alpha, "payoff_surplus.grid.MPI.gnuplot", true);
        myBSSolver->printSparseGrid(*alpha, "payoff_nodal.grid.MPI.gnuplot", false);

        if (coordsType == "log") {
          myBSSolver->printSparseGridExpTransform(*alpha, "payoff_surplus_cart.grid.MPI.gnuplot", true);
          myBSSolver->printSparseGridExpTransform(*alpha, "payoff_nodal_cart.grid.MPI.gnuplot", false);
        }

        if (coordsType == "PAT") {
          myBSSolver->printSparseGridPAT(*alpha, "payoff_surplus_cart.PAT.grid.MPI.gnuplot", true);
          myBSSolver->printSparseGridPAT(*alpha, "payoff_nodal_cart.PAT.grid.MPI.gnuplot", false);
        }
      }
    }

    // Communicate grid
    if (SGPP::parallel::myGlobalMPIComm->getMyRank() == 0) {
      std::string serialized_grid = myBSSolver->getGrid();

      SGPP::parallel::myGlobalMPIComm->broadcastGrid(serialized_grid);
    } else {
      // Now receive the grid
      std::string serialized_grid = "";

      SGPP::parallel::myGlobalMPIComm->receiveGrid(serialized_grid);
      myBSSolver->setGrid(serialized_grid);

      alpha = new SGPP::base::DataVector(myBSSolver->getNumberGridPoints());

      // Set stochastic data
      myBSSolver->setStochasticData(mu, sigma, rho, r);

      // Set Payoff type
      myBSSolver->setPayoffType(payoffType);
    }

    // Communicate coefficients
    SGPP::parallel::myGlobalMPIComm->broadcastGridCoefficientsFromRank0(*alpha);

    // Start solving the Black Scholes Equation
    if (Solver == "ExEul") {
      myBSSolver->solveExplicitEuler(timesteps, stepsize, CGiterations, CGepsilon, *alpha, false, false, 20);
    } else if (Solver == "ImEul") {
      myBSSolver->solveImplicitEuler(timesteps, stepsize, CGiterations, CGepsilon, *alpha, false, false, 20);
    } else if (Solver == "CrNic") {
      myBSSolver->solveCrankNicolson(timesteps, stepsize, CGiterations, CGepsilon, *alpha, CRNIC_IMEUL_STEPS);
    } else {
      std::cout << "!!!! You have chosen an unsupported solver type !!!!" << std::endl;
    }

    // print solved grid only on rank 0
    if (SGPP::parallel::myGlobalMPIComm->getMyRank() == 0) {
      if (dim < 3) {
        // Print the solved Black Scholes Equation into a gnuplot file
        myBSSolver->printGrid(*alpha, 20, "solvedBS.MPI.gnuplot");
      }

      if (dim < 4) {
        myBSSolver->printSparseGrid(*alpha, "solvedBS_surplus.grid.MPI.gnuplot", true);
        myBSSolver->printSparseGrid(*alpha, "solvedBS_nodal.grid.MPI.gnuplot", false);

        if (coordsType == "log") {
          myBSSolver->printSparseGridExpTransform(*alpha, "solvedBS_surplus_cart.grid.MPI.gnuplot", true);
          myBSSolver->printSparseGridExpTransform(*alpha, "solvedBS_nodal_cart.grid.MPI.gnuplot", false);
        }

        if (coordsType == "PAT") {
          myBSSolver->printSparseGridPAT(*alpha, "solvedBS_surplus_cart.PAT.grid.MPI.gnuplot", true);
          myBSSolver->printSparseGridPAT(*alpha, "solvedBS_nodal_cart.PAT.grid.MPI.gnuplot", false);
        }
      }

      // Test option @ the money
      std::vector<double> point;

      for (size_t x = 0; x < dim; x++) {
        point.push_back(dStrike);
      }

      std::cout << "Optionprice at testpoint (Strike): " << myBSSolver->evalOption(point, *alpha) << std::endl << std::endl;

      if (i == start_l) {
        // if transformed Black Scholes Equation is used -> transform evaluation domain
        for (size_t v = 0; v < EvalPoints.getNrows(); v++) {
          SGPP::base::DataVector r(dim);
          EvalPoints.getRow(v, r);
          myBSSolver->transformPoint(r);
          EvalPoints.setRow(v, r);
        }
      }

      // Evaluate Cuboid
      SGPP::base::DataVector Prices(EvalPoints.getNrows());
      myBSSolver->evaluateCuboid(*alpha, Prices, EvalPoints);
      results.push_back(Prices);

      // write solution in a additional file
      std::stringstream level_string;
      level_string << i;
      writeDataVector(Prices, tFileEvalCuboidValues + ".level_" + level_string.str());
      writeDataVector(Prices, tFileEvalCuboidValues);

      if (i > start_l) {
        std::cout << "=====================================================================" << std::endl;
        std::cout << "=====================================================================" << std::endl << std::endl;
        std::cout << "Calculating norms of relative errors to a grid" << std::endl;
        std::cout << "with " << i << " levels and testing-coboid" << std::endl << std::endl;

        double oldMaxNorm = 0.0;
        double oldTwoNorm = 0.0;

        // Calculate relative errors and some norms
        std::cout << "Level, max-norm(rel-error), two-norm(rel-error), rate max-norm,  rate two-norm" << std::endl;
        std::cout << "------------------------------------------------------------------------------" << std::endl;

        for (size_t j = 0; j < i - start_l; j++) {
          SGPP::base::DataVector maxLevel(results[i - start_l]);
          SGPP::base::DataVector relError(results[j]);
          double maxNorm = 0.0;
          double l2Norm = 0.0;

          // calculate relative error
          relError.sub(maxLevel);
          relError.componentwise_div(maxLevel);

          // calculate max. norm of relative error
          maxNorm = relError.maxNorm();

          // calculate two norm of relative error
          l2Norm = relError.RMSNorm();

          // Printing norms
          //std::cout << "Level " << j + start_l << ": max-norm(rel-error)=" << maxNorm << "; two-norm(rel-error)=" << l2Norm << "; rate max-norm: " << log(oldMaxNorm/maxNorm) << "; rate two-norm: " << log(oldTwoNorm/l2Norm) << std::endl;
          std::cout << j + start_l << ", " << maxNorm << ", " << l2Norm << ", " << log(oldMaxNorm / maxNorm) << ", " << log(oldTwoNorm / l2Norm) << std::endl;

          oldMaxNorm = maxNorm;
          oldTwoNorm = l2Norm;
        }
      }

      std::cout << std::endl << std::endl << std::endl;
    }

    myBSSolver->deleteGrid();
    delete alpha;
    alpha = NULL;
  }

  delete myBSSolver;
}

/**
 * Do a Black Scholes solver test with n assets (ND Sparse Grid) European call option, with Initial
 * Grid Refinement
 *
 * @param d the number of dimensions used in the Sparse Grid
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
 * @param Solver specifies the sovler that should be used, ExEul, ImEul and CrNic are the possibilities
 * @param refinementMode the mode selected for surplus refinement: available: classic, maxLevel
 * @param maxRefineLevel ignored for refinement mode classic, in maxLevel: max. level to which the grid is refined
 * @param numRefinePoints number of points that should be refined in each refine iteration before Black Scholes Equation is solved: -1 try to refine all points steered by threshold
 * @param nIterAdaptSteps number of the iterative Grid Refinement that should be executed
 * @param dRefineThreshold Threshold for a point's surplus for refining this point
 * @param useCoarsen specifies if the grid should be coarsened between timesteps
 * @param adaptSolvingMode specifies which adaptive methods are applied during solving the BS Equation
 * @param coarsenThreshold Threshold to decide, if a grid point should be deleted
 * @param coordsType set the type of coordinates that should be used: cart, log, PAT
 * @param useNormalDist enable local initial refinement based on a normal distribution
 */
void testNUnderlyingsAdaptSurplus(size_t d, size_t l, std::string fileStoch, std::string fileBound, double dStrike,
                                  std::string payoffType, double riskfree, size_t timeSt, double dt, size_t CGIt, double CGeps,
                                  std::string Solver, std::string refinementMode, int numRefinePoints, size_t maxRefineLevel, size_t nIterAdaptSteps, double dRefineThreshold,
                                  bool useCoarsen, std::string adaptSolvingMode, double coarsenThreshold, std::string coordsType, bool useNormalDist) {
  size_t dim = d;
  int level = (int)l;
  size_t timesteps = timeSt;
  double stepsize = dt;
  size_t CGiterations = CGIt;
  double CGepsilon = CGeps;

  double refineSigma = DFLT_SIGMA_REFINE_NORMDIST;

  SGPP::base::DataVector mu(dim);
  SGPP::base::DataVector sigma(dim);
  SGPP::base::DataMatrix rho(dim, dim);

  SGPP::base::DataVector* alpha = NULL;

  SGPP::base::DataVector EvalCuboidValues(1);
  SGPP::base::DataMatrix EvalCuboid(1, dim);

  int retCuboid = -1;
  int retCuboidValues = -1;

  double r = riskfree;

  if (readStochasticData(fileStoch, dim, mu, sigma, rho) != 0) {
    return;
  }

  SGPP::parallel::BlackScholesSolverMPI* myBSSolver = NULL;

  if (coordsType == "log") {
    myBSSolver = new SGPP::parallel::BlackScholesSolverMPI(true);
  } else if (coordsType == "cart") {
    myBSSolver = new SGPP::parallel::BlackScholesSolverMPI(false);
  } else if (coordsType == "PAT") {
    myBSSolver = new SGPP::parallel::BlackScholesSolverMPI(true, true);
  } else {
    // Write Error
    std::cout << "Unsupported grid transformation!" << std::endl;
    std::cout << std::endl << std::endl;

    if (SGPP::parallel::myGlobalMPIComm->getMyRank() == 0) {
      writeHelp();
    }

    SGPP::parallel::myGlobalMPIComm->Abort();
  }

  // Enable Coarsening
  if (useCoarsen == true) {
    myBSSolver->setEnableCoarseningData(adaptSolvingMode, refinementMode, (int)maxRefineLevel, -1, coarsenThreshold, dRefineThreshold);
  }

  if (SGPP::parallel::myGlobalMPIComm->getMyRank() == 0) {
    SGPP::base::DimensionBoundary* myBoundaries = new SGPP::base::DimensionBoundary[dim];

    if (readBoudingBoxData(fileBound, dim, myBoundaries) != 0) {
      return;
    }

    SGPP::base::BoundingBox* myBoundingBox = new SGPP::base::BoundingBox(dim, myBoundaries);
    delete[] myBoundaries;

    // init Screen Object
    myBSSolver->initScreen();

    // Construct a grid
    myBSSolver->constructGrid(*myBoundingBox, level);
    delete myBoundingBox;

    // init the basis functions' coefficient vector
    alpha = new SGPP::base::DataVector(myBSSolver->getNumberGridPoints());

    // Set stochastic data
    myBSSolver->setStochasticData(mu, sigma, rho, r);

    // Init the grid with on payoff function
    myBSSolver->initGridWithPayoff(*alpha, dStrike, payoffType);

    std::vector<double> norm_mu;
    std::vector<double> norm_sigma;

    // estimate refine sigma from evaluation cuboid
    // read Evaluation cuboid
    retCuboid = readEvalutionCuboid(EvalCuboid, tFileEvalCuboid, dim);

    // read reference values for evaluation cuboid
    retCuboidValues = readOptionsValues(EvalCuboidValues, tFileEvalCuboidValues);

    if (EvalCuboid.getNrows() != EvalCuboidValues.getSize()) {
      retCuboid = -1;
      retCuboidValues = -1;
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

    std::cout << "Initial Grid size: " << myBSSolver->getNumberGridPoints() << std::endl;
    std::cout << "Initial Grid size (inner): " << myBSSolver->getNumberInnerGridPoints() << std::endl << std::endl << std::endl;

    // refine the grid to approximate the singularity in the start solution better
    if (refinementMode == "classic") {
      for (size_t i = 0 ; i < nIterAdaptSteps; i++) {
        std::cout << "Refining Grid..." << std::endl;

        if (useNormalDist == true) {
          myBSSolver->refineInitialGridSurplusSubDomain(*alpha, numRefinePoints, dRefineThreshold, norm_mu, norm_sigma);
        } else {
          myBSSolver->refineInitialGridSurplus(*alpha, numRefinePoints, dRefineThreshold);
        }

        myBSSolver->initGridWithPayoff(*alpha, dStrike, payoffType);
        std::cout << "Refined Grid size: " << myBSSolver->getNumberGridPoints() << std::endl;
        std::cout << "Refined Grid size (inner): " << myBSSolver->getNumberInnerGridPoints() << std::endl;
      }

    } else if (refinementMode == "maxLevel") {
      size_t oldGridSize = 0;
      size_t newGridSize = myBSSolver->getNumberGridPoints();
      size_t addedGridPoint = 0;
      size_t stepCounter = 0;

      if (nIterAdaptSteps > 0) {
        do {
          oldGridSize = newGridSize;
          std::cout << "Refining Grid..." << std::endl;

          if (useNormalDist == true) {
            myBSSolver->refineInitialGridSurplusToMaxLevelSubDomain(*alpha, dRefineThreshold, (int)maxRefineLevel, norm_mu, norm_sigma);
          } else {
            myBSSolver->refineInitialGridSurplusToMaxLevel(*alpha, dRefineThreshold, (int)maxRefineLevel);
          }

          myBSSolver->initGridWithPayoff(*alpha, dStrike, payoffType);
          std::cout << "Refined Grid size: " << myBSSolver->getNumberGridPoints() << std::endl;
          std::cout << "Refined Grid size (inner): " << myBSSolver->getNumberInnerGridPoints() << std::endl;
          newGridSize = myBSSolver->getNumberGridPoints();
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
      myBSSolver->printGrid(*alpha, 20, "payoff.MPI.gnuplot");
    }

    if (dim < 4) {
      myBSSolver->printSparseGrid(*alpha, "payoff_surplus.grid.MPI.gnuplot", true);
      myBSSolver->printSparseGrid(*alpha, "payoff_nodal.grid.MPI.gnuplot", false);

      if (coordsType == "log") {
        myBSSolver->printSparseGridExpTransform(*alpha, "payoff_surplus_cart.grid.MPI.gnuplot", true);
        myBSSolver->printSparseGridExpTransform(*alpha, "payoff_nodal_cart.grid.MPI.gnuplot", false);
      }

      if (coordsType == "PAT") {
        myBSSolver->printSparseGridPAT(*alpha, "payoff_surplus_cart.PAT.grid.MPI.gnuplot", true);
        myBSSolver->printSparseGridPAT(*alpha, "payoff_nodal_cart.PAT.grid.MPI.gnuplot", false);
      }
    }
  }

  // Communicate grid
  if (SGPP::parallel::myGlobalMPIComm->getMyRank() == 0) {
    std::string serialized_grid = myBSSolver->getGrid();

    SGPP::parallel::myGlobalMPIComm->broadcastGrid(serialized_grid);
  } else {
    // Now receive the grid
    std::string serialized_grid = "";

    SGPP::parallel::myGlobalMPIComm->receiveGrid(serialized_grid);
    myBSSolver->setGrid(serialized_grid);

    alpha = new SGPP::base::DataVector(myBSSolver->getNumberGridPoints());

    // Set stochastic data
    myBSSolver->setStochasticData(mu, sigma, rho, r);

    // Set Payoff type
    myBSSolver->setPayoffType(payoffType);
  }

  // Communicate coefficients
  SGPP::parallel::myGlobalMPIComm->broadcastGridCoefficientsFromRank0(*alpha);

  // Start solving the Black Scholes Equation
  if (Solver == "ExEul") {
    myBSSolver->solveExplicitEuler(timesteps, stepsize, CGiterations, CGepsilon, *alpha, false, false, 20);
  } else if (Solver == "ImEul") {
    myBSSolver->solveImplicitEuler(timesteps, stepsize, CGiterations, CGepsilon, *alpha, false, false, 20);
  } else if (Solver == "CrNic") {
    myBSSolver->solveCrankNicolson(timesteps, stepsize, CGiterations, CGepsilon, *alpha, CRNIC_IMEUL_STEPS);
  } else {
    std::cout << "!!!! You have chosen an unsupported solver type !!!!" << std::endl;
  }

  if (SGPP::parallel::myGlobalMPIComm->getMyRank() == 0) {
    if (dim < 3) {
      // Print the solved Black Scholes Equation into a gnuplot file
      myBSSolver->printGrid(*alpha, 20, "solvedBS.MPI.gnuplot");
    }

    if (dim < 4) {
      myBSSolver->printSparseGrid(*alpha, "solvedBS_surplus.grid.MPI.gnuplot", true);
      myBSSolver->printSparseGrid(*alpha, "solvedBS_nodal.grid.MPI.gnuplot", false);

      if (coordsType == "log") {
        myBSSolver->printSparseGridExpTransform(*alpha, "solvedBS_surplus_cart.grid.MPI.gnuplot", true);
        myBSSolver->printSparseGridExpTransform(*alpha, "solvedBS_nodal_cart.grid.MPI.gnuplot", false);
      }

      if (coordsType == "PAT") {
        myBSSolver->printSparseGridPAT(*alpha, "solvedBS_surplus_cart.PAT.grid.MPI.gnuplot", true);
        myBSSolver->printSparseGridPAT(*alpha, "solvedBS_nodal_cart.PAT.grid.MPI.gnuplot", false);
      }
    }

    // Test option @ the money
    std::vector<double> point;

    for (size_t i = 0; i < d; i++) {
      point.push_back(dStrike);
    }

    std::cout << "Optionprice at testpoint (Strike): " << myBSSolver->evalOption(point, *alpha) << std::endl << std::endl;

    // calculate relative errors
    ////////////////////////////
    double maxNorm = 0.0;
    double l2Norm = 0.0;

    if (retCuboid == 0 && retCuboidValues == 0) {
      // If the transformed Black Scholes Equation is used -> transform Eval-cuboid
      for (size_t v = 0; v < EvalCuboid.getNrows(); v++) {
        SGPP::base::DataVector r(dim);
        EvalCuboid.getRow(v, r);
        myBSSolver->transformPoint(r);
        EvalCuboid.setRow(v, r);
      }

      std::cout << "Calculating relative errors..." << std::endl;
      // Evaluate Cuboid
      SGPP::base::DataVector Prices(EvalCuboid.getNrows());
      myBSSolver->evaluateCuboid(*alpha, Prices, EvalCuboid);

      SGPP::base::DataVector relError(Prices);

      // calculate relative error
      relError.sub(EvalCuboidValues);
      relError.componentwise_div(EvalCuboidValues);

      // calculate max. norm of relative error
      maxNorm = relError.maxNorm();

      // calculate two norm of relative error
      l2Norm = relError.RMSNorm();

      // Printing norms
      std::cout << "Results: max-norm(rel-error)=" << maxNorm << "; two-norm(rel-error)=" << l2Norm << std::endl;

      // reprint data with prefix -> can be easily grep-ed
      std::cout << std::endl << std::endl;
    } else {
      std::cout << "Couldn't open evaluation cuboid data -> skipping tests!" << std::endl << std::endl;
    }

    std::cout << "$ Startlevel: " << level << "; RefineMode: " << refinementMode << "; MaxRefLevel: " << maxRefineLevel << std::endl;
    std::string normDistrefine;

    if (useNormalDist == true) {
      std::stringstream normDistRefineStream;
      normDistRefineStream << "solveNDadaptSurplusSubDomain;" << dStrike << ";" << refineSigma;
      normDistrefine = normDistRefineStream.str();
      std::cout << "$ AdaptSurplus-Mode: solveNDadaptSurplusSubDomain" << std::endl;
      std::cout << "$ Refine mu = " << dStrike << "; Refine sigma = " << refineSigma << std::endl;
    } else {
      normDistrefine = "solveNDadaptSurplus;-1.0;1.0";
      std::cout << "$ AdaptSurplus-Mode: solveNDadaptSurplus" << std::endl;
    }

    std::cout << "$ NumRefinements: " << nIterAdaptSteps << "; RefineThreshd: " << dRefineThreshold << std::endl;
    std::cout << "$ AdpatSolveMode: " << adaptSolvingMode << "; CoarsenThreshd: " << coarsenThreshold << std::endl;
    std::cout << "$ Start #gridpoints (inner): " << myBSSolver->getStartInnerGridSize() << std::endl;
    std::cout << "$ Final #gridpoints (inner): " << myBSSolver->getFinalInnerGridSize() << std::endl;
    std::cout << "$ Average #gridpoints (inner): " << myBSSolver->getAverageInnerGridSize() << std::endl;
    std::cout << "$ Needed iterations: " << myBSSolver->getNeededIterationsToSolve() << "; Needed time: " << myBSSolver->getNeededTimeToSolve() << std::endl;
    std::cout << "$ Results: max-norm(rel-error)=" << maxNorm << "; two-norm(rel-error)=" << l2Norm << std::endl;
    std::cout << "$ Optionprice at testpoint (Strike): " << myBSSolver->evalOption(point, *alpha) << std::endl;
    std::cout << "$ CSV-DATA: " << level << ";" << refinementMode << ";" << maxRefineLevel << ";" << nIterAdaptSteps
              << ";" << dRefineThreshold << ";" << normDistrefine << ";" << adaptSolvingMode << ";" << coarsenThreshold
              << ";" << myBSSolver->getStartInnerGridSize() << ";" << myBSSolver->getFinalInnerGridSize()
              << ";" << myBSSolver->getAverageInnerGridSize() << ";" << myBSSolver->getNeededIterationsToSolve()
              << ";" << myBSSolver->getNeededTimeToSolve() << ";" << maxNorm << ";" << l2Norm << std::endl;
    std::cout << std::endl << std::endl;
  }

  delete myBSSolver;
  delete alpha;
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
  int mpi_myid;
  int mpi_ranks;

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_ranks);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_myid);
  SGPP::parallel::myGlobalMPIComm = new SGPP::parallel::MPICommunicator(mpi_myid, mpi_ranks);

  std::streambuf* stdoutBuf = std::cout.rdbuf();
  std::ofstream dummy_out("/dev/null");

  if (mpi_myid != 0) { // disable output for all processes but proc 0
    std::cout.rdbuf(dummy_out.rdbuf());
  }

  if (argc == 1) {
    if (mpi_myid == 0) {
      writeHelp();
    }

    SGPP::parallel::myGlobalMPIComm->Abort();
    return 0;
  }

  option.assign(argv[1]);

  if (option == "solveNDanalyze") {
    if (argc != 17) {
      if (mpi_myid == 0) {
        writeHelp();
      }

      SGPP::parallel::myGlobalMPIComm->Abort();
      return 0;
    } else {
      std::string fileStoch;
      std::string fileBound;
      double dStrike;
      std::string fileAnalyze;
      std::string ani;
      std::string solver;
      std::string payoff;

      fileStoch.assign(argv[7]);
      fileBound.assign(argv[6]);
      dStrike = atof(argv[8]);
      fileAnalyze.assign(argv[16]);
      payoff.assign(argv[9]);
      solver.assign(argv[13]);

      std::string coordsType;
      coordsType.assign(argv[2]);

      testNUnderlyingsAnalyze(atoi(argv[3]), atoi(argv[4]), atoi(argv[5]), fileStoch, fileBound, dStrike, payoff, atof(argv[10]), (size_t)(atof(argv[11]) / atof(argv[12])), atof(argv[12]), atoi(argv[14]), atof(argv[15]), solver, fileAnalyze, coordsType);
    }
  } else if (option == "solveNDadaptSurplus" || option == "solveNDadaptSurplusSubDomain") {
    if (argc != 21) {
      if (mpi_myid == 0) {
        writeHelp();
      }

      SGPP::parallel::myGlobalMPIComm->Abort();
      return 0;
    } else {
      bool isNormalDist = false;

      if (option == "solveNDadaptSurplusSubDomain") {
        isNormalDist = true;
      }

      std::string fileStoch;
      std::string fileBound;
      double dStrike;
      std::string ani;
      std::string solver;
      std::string payoff;
      std::string refinementMode;
      std::string adaptSolveMode;

      fileStoch.assign(argv[6]);
      fileBound.assign(argv[5]);
      dStrike = atof(argv[7]);
      payoff.assign(argv[8]);
      solver.assign(argv[12]);
      refinementMode.assign(argv[15]);
      adaptSolveMode.assign(argv[19]);

      std::string coordsType;
      coordsType.assign(argv[2]);

      if (refinementMode != "maxLevel" && refinementMode != "classic") {
        std::cout << "Unsupported refinement type! classic or maxLevel are supported!" << std::endl;
        std::cout << std::endl << std::endl;

        if (mpi_myid == 0) {
          writeHelp();
        }

        SGPP::parallel::myGlobalMPIComm->Abort();
        return 0;
      }

      bool useAdaptSolve = false;

      if (adaptSolveMode == "coarsen" || adaptSolveMode == "refine" || adaptSolveMode == "coarsenNrefine") {
        useAdaptSolve = true;
      } else if (adaptSolveMode == "none") {
        useAdaptSolve = false;
      } else {
        std::cout << "Unsupported adapt solve mode! none, coarsen, refine or coarsenNrefine are supported!" << std::endl;
        std::cout << std::endl << std::endl;

        if (mpi_myid == 0) {
          writeHelp();
        }

        SGPP::parallel::myGlobalMPIComm->Abort();
        return 0;
      }

      testNUnderlyingsAdaptSurplus(atoi(argv[3]), atoi(argv[4]), fileStoch, fileBound, dStrike, payoff, atof(argv[9]), (size_t)(atof(argv[10]) / atof(argv[11])), atof(argv[11]), atoi(argv[13]), atof(argv[14]), solver, refinementMode, -1, atoi(argv[16]), atoi(argv[17]), atof(argv[18]), useAdaptSolve, adaptSolveMode, atof(argv[20]), coordsType, isNormalDist);
    }
  } else {
    if (mpi_myid == 0) {
      writeHelp();
    }

    SGPP::parallel::myGlobalMPIComm->Abort();
    return 0;
  }

  delete SGPP::parallel::myGlobalMPIComm;

  if (mpi_myid != 0) { // restore stdout buffer
    std::cout.rdbuf(stdoutBuf);
  }

  MPI_Finalize();

  return 0;
}
