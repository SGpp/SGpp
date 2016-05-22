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

#define CRNIC_IMEUL_STEPS 3

int readBoudingBoxData(std::string tFile, size_t numAssests,
                       sgpp::base::BoundingBox1D* BoundaryArray) {
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

double calculatetheta(double a, double sigma, double T, int count, double stepsize) {
  double theta = 0;
  return theta = 0.04 * a + pow(sigma, 2.0) * (1 - exp(-2 * a * (T - count * stepsize))) / (2 * a);
}
/**
 * Do a Hull White solver test with n assets (ND Sparse Grid) European call option
 *
 * @param l the number of levels used in the Sparse Grid
 * @param sigma sigma
 * @param a a
 * @param fileBound the grid's bounding box - domain boundary(min,max)
 * @param payoffType method that is used to determine the multidimensional payoff function
 * @param timeSt the number of timesteps that are executed during the solving process
 * @param dt the size of delta t in the ODE solver
 * @param CGIt the maximum number of Iterations that are executed by the CG/BiCGStab
 * @param CGeps the epsilon used in the CG/BiCGStab
 * @param Solver specifies the sovler that should be used, ExEul, ImEul and CrNic are the
 * possibilities
 * @param t current time
 * @param T time to maturity
 * @param dStrike strike
 */
void testHullWhite(int l, double sigma, double a, std::string fileBound, std::string payoffType,
                   size_t timeSt, double dt, size_t CGIt, double CGeps, std::string Solver,
                   double t, double T, double dStrike) {
  int level = l;
  size_t timesteps = timeSt;
  double stepsize = dt;
  size_t CGiterations = CGIt;
  double CGepsilon = CGeps;

  sgpp::base::BoundingBox1D* myBoundaries = new sgpp::base::BoundingBox1D[1];

  if (readBoudingBoxData(fileBound, 1, myBoundaries) != 0) {
    return;
  }

  sgpp::finance::HullWhiteSolver* myHWSolver = new sgpp::finance::HullWhiteSolver();
  sgpp::base::BoundingBox* myBoundingBox = new sgpp::base::BoundingBox(1, myBoundaries);
  delete[] myBoundaries;

  // init Screen Object
  myHWSolver->initScreen();

  // Construct a grid
  myHWSolver->constructGrid(*myBoundingBox, level);

  // init the basis functions' coefficient vector
  sgpp::base::DataVector* alpha = new sgpp::base::DataVector(myHWSolver->getNumberGridPoints());

  std::cout << "Grid has " << level << " Levels" << std::endl;
  std::cout << "Initial Grid size: " << myHWSolver->getNumberGridPoints() << std::endl;
  std::cout << "Initial Grid size (inner): " << myHWSolver->getNumberInnerGridPoints() << std::endl
            << std::endl
            << std::endl;

  // Init the grid with on payoff function
  myHWSolver->initGridWithPayoff(*alpha, dStrike, payoffType, sigma, a, t, T);

  // Gridpoints @Money
  //  std::cout << "Gridpoints @Money: " << myBSSolver->getGridPointsAtMoney(payoffType, dStrike,
  //  DFLT_EPS_AT_MONEY) << std::endl << std::endl << std::endl;

  // Print the payoff function into a gnuplot file

  myHWSolver->printGrid(*alpha, 20, "payoffHW.gnuplot");
  myHWSolver->printSparseGrid(*alpha, "payoffHW_surplus.grid.gnuplot", true);
  myHWSolver->printSparseGrid(*alpha, "payoffHW_nodal.grid.gnuplot", false);

  // Set stochastic data
  // myHWSolver->setStochasticData(theta, sigma, a);

  // Start solving the Black Scholes Equation
  if (Solver == "ExEul") {
    double theta = 0;
    int count = 0;

    for (int i = 0; i < T / stepsize; i++) {
      theta = calculatetheta(a, sigma, T, count, stepsize);
      myHWSolver->setStochasticData(theta, sigma, a);
      myHWSolver->solveExplicitEuler(timesteps, stepsize, CGiterations, CGepsilon, *alpha, false,
                                     false);
      count = count + 1;
    }
  } else if (Solver == "ImEul") {
    double theta = 0;
    int count = 0;

    for (int i = 0; i < T / stepsize; i++) {
      theta = calculatetheta(a, sigma, T, count, stepsize);
      myHWSolver->setStochasticData(theta, sigma, a);
      myHWSolver->solveImplicitEuler(timesteps, stepsize, CGiterations, CGepsilon, *alpha, false,
                                     false);
      count = count + 1;
    }
  } else if (Solver == "CrNic") {
    myHWSolver->solveCrankNicolson(timesteps, stepsize, CGiterations, CGepsilon, *alpha,
                                   CRNIC_IMEUL_STEPS);
  }
  /*else if (Solver == "AdBas")
  {
    myBSSolver->solveAdamsBashforth(timesteps, stepsize, CGiterations, CGepsilon, *alpha, false,
  false, 20);
  }
  else if (Solver == "VaTim")
  {
    myBSSolver->solveVarTimestep(timesteps, stepsize, CGiterations, CGepsilon, *alpha, false, false,
  20);
  }*/ else {
    std::cout << "!!!! You have chosen an unsupported solver type !!!!" << std::endl;
  }

  if (Solver == "ExEul" || Solver == "ImEul" || Solver == "CrNic") {
    // Print the solved Black Scholes Equation into a gnuplot file
    myHWSolver->printGrid(*alpha, 20, "solvedHW.gnuplot");
    myHWSolver->printSparseGrid(*alpha, "solvedHW_surplus.grid.gnuplot", true);
    myHWSolver->printSparseGrid(*alpha, "solvedHW_nodal.grid.gnuplot", false);
  }

  delete alpha;
  delete myHWSolver;
  delete myBoundingBox;

  // std::cout << "Nothing has been implemented so far " << std::endl;
}

/**
 * Calls the writeHelp method in the HullWhiteSolver Object
 * after creating a screen.
 */

void writeHelp() {
  sgpp::finance::HullWhiteSolver* myHWSolver = new sgpp::finance::HullWhiteSolver();

  myHWSolver->initScreen();

  delete myHWSolver;

  std::stringstream mySStream;

  mySStream << "Some instructions for the use of 1D Hull White Solver:" << std::endl;
  mySStream << "------------------------------------------------------" << std::endl << std::endl;
  mySStream << "Available execution modes are:" << std::endl;
  mySStream << "  solve1DHullWhite" << std::endl;
  mySStream << "Execution modes descriptions:" << std::endl;
  mySStream << "-----------------------------------------------------" << std::endl;
  mySStream << "solve1DHullWhite" << std::endl << "------" << std::endl;
  mySStream << "the following options must be specified:" << std::endl;
  mySStream << "  level: number of levels within the Sparse Grid" << std::endl;
  mySStream << "  value of sigma: sigma value-determine overall level of volatility" << std::endl;
  mySStream << "  value of a: a" << std::endl;
  mySStream << "  file_Boundaries: file that contains the bounding box" << std::endl;
  mySStream << "  payoff_func: function for n-d payoff: std_euro_{call|put}" << std::endl;
  mySStream << "  simeSt: number of time steps of doing HW, default: 1, the total time steps "
               "is defined by T/dT"
            << std::endl;
  mySStream << "  dT: timestep size" << std::endl;
  mySStream << "  CGIterations: Maxmimum number of iterations used in CG mehtod" << std::endl;
  mySStream << "  CGEpsilon: Epsilon used in CG" << std::endl;
  mySStream << "  Solver: the solver to use: ExEul, ImEul or CrNic" << std::endl;
  mySStream << "  t: the current time" << std::endl;
  mySStream << "  T: time to maturity" << std::endl;
  mySStream << "  Strike: the strike" << std::endl;

  mySStream << std::endl;
  mySStream << "Example:" << std::endl;
  mySStream << "5 0.01 0.1 bound.data "
            << " std_euro_call "
            << "1.0 "
            << "0.01 "
            << "400 "
            << "0.000001 "
            << "ImEul "
            << "0.2 "
            << "1.0 "
            << "0.6 " << std::endl;
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

  if (option == "solve1DHullWhite") {
    if (argc != 15) {
      writeHelp();
    } else {
      std::string solver;
      std::string payoff;
      double sigma;
      double a;
      double dStrike;
      std::string fileBound;
      sigma = atof(argv[3]);
      a = atof(argv[4]);
      fileBound.assign(argv[5]);
      payoff.assign(argv[6]);
      solver.assign(argv[11]);
      dStrike = atof(argv[14]);

      testHullWhite(atoi(argv[2]), sigma, a, fileBound, payoff, atoi(argv[7]), atof(argv[8]),
                    atoi(argv[9]), atof(argv[10]), solver, atof(argv[12]), atof(argv[13]), dStrike);
    }
  } else {
    writeHelp();
  }

  return 0;
}
