// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#define CRNIC_IMEUL_STEPS 3
#define GUNPLOT_RESOLUTION 51
#define SOLUTION_FRAMES 100

#define DIV_SIGMA 4.0
#define DISTRI_FACTOR 5.0

#define EXPORT_MATRIX_FILES

#include <sgpp_base.hpp>
#include <sgpp_pde.hpp>
#include <sgpp_finance.hpp>
#include <sgpp_parallel.hpp>
#include <sgpp_solver.hpp>
#include <sgpp_datadriven.hpp>

#include <cstdlib>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <string>
#include <vector>
/**
 * Calls the writeHelp method in the BlackScholesSolver Object
 * after creating a screen.
 */
void writeHelp() {
  SGPP::pde::HeatEquationSolverWithStretching* myHESolver =
      new SGPP::pde::HeatEquationSolverWithStretching();

  myHESolver->initScreen();

  delete myHESolver;

  std::stringstream mySStream;

  mySStream << "Some instructions for the use of Poisson/ Heat Equation Solver:" << std::endl;
  mySStream << "---------------------------------------------------------------" << std::endl
            << std::endl;
  mySStream << "Available execution modes are:" << std::endl;
  mySStream << "  HeatEquation        Solves Heat Equation on a quadratic" << std::endl;
  mySStream << "                      d-dimensional domain" << std::endl << std::endl;
  mySStream << "  PoissonEquation     Solves Poisson Equation on a quadratic" << std::endl;
  mySStream << "                      d-dimensional domain" << std::endl << std::endl << std::endl;

  mySStream << "Execution modes descriptions:" << std::endl;
  mySStream << "-----------------------------------------------------" << std::endl;
  mySStream << "HeatEquation" << std::endl << "------" << std::endl;
  mySStream << "the following options must be specified:" << std::endl;
  mySStream << "  dim: the number of dimensions of Sparse Grid" << std::endl;
  mySStream << "  level: number of levels within the Sparse Grid" << std::endl;
  mySStream << "  left_bound: x_i of left boundary" << std::endl;
  mySStream << "  right_bound: x_i of right boundary" << std::endl;
  mySStream << "  a: thermal diffusivity" << std::endl;
  mySStream << "  initHeat: initial heat distribution" << std::endl;
  mySStream << "  T: time to solve" << std::endl;
  mySStream << "  dT: timestep size" << std::endl;
  mySStream << "  Solver: the solver to use: ExEul, ImEul, CrNic" << std::endl;
  mySStream << "  CGEpsilon: Epsilon used in CG" << std::endl;
  mySStream << "  CGIterations: Maxmimum number of iterations used in CG mehtod" << std::endl;
  mySStream << "  stretchingMode: gives the stretching mode, can be either analytic or discrete"
            << std::endl;
  mySStream << "  file_stretch: file containing the stretching data, file format different "
               "depending on the mode"
            << std::endl;
  mySStream << std::endl;
  mySStream << "Example:" << std::endl;
  mySStream << "HESolver HeatEquation 3 5 0.0 3.0 1.0 smooth 1.0 0.1 ImEul 0.00001 400 analytic "
               "fileStretch.data"
            << std::endl;
  mySStream << std::endl;
  mySStream << "Remark: This test generates following files (gnuplot):" << std::endl;
  mySStream << "  heatStart.gnuplot: the start condition" << std::endl;
  mySStream << "  heatSolved.gnuplot: the numerical solution" << std::endl;
  mySStream << std::endl << std::endl;

  mySStream << "PoissonEquation" << std::endl << "------" << std::endl;
  mySStream << "the following options must be specified:" << std::endl;
  mySStream << "  dim: the number of dimensions of Sparse Grid" << std::endl;
  mySStream << "  level: number of levels within the Sparse Grid" << std::endl;
  mySStream << "  left_bound: x_i of left boundary" << std::endl;
  mySStream << "  right_bound: x_i of right boundary" << std::endl;
  mySStream << "  initHeat: initial heat distribution" << std::endl;
  mySStream << "  CGEpsilon: Epsilon used in CG" << std::endl;
  mySStream << "  CGIterations: Maxmimum number of iterations used in CG mehtod" << std::endl;
  mySStream << "  stretchingMode: gives the stretching mode, can be either analytic or discrete"
            << std::endl;
  mySStream << "  file_stretch: file containing the stretching data, file format different "
               "depending on the mode"
            << std::endl;
  mySStream << std::endl;
  mySStream << "Example:" << std::endl;
  mySStream << "HESolver PoissonEquation 3 5 0.0 3.0 smooth 0.00001 400 analytic fileStretch.data"
            << std::endl;
  mySStream << std::endl;
  mySStream << "Remark: This test generates following files (gnuplot):" << std::endl;
  mySStream << "  poissonStart.gnuplot: the start condition" << std::endl;
  mySStream << "  poissonSolved.gnuplot: the numerical solution" << std::endl;
  mySStream << std::endl << std::endl;

  std::cout << mySStream.str() << std::endl;
}

int readDiscreteStretchingData(std::string tFile, size_t numAssests,
                               std::vector<double>* discreteCoordinates) {
  std::fstream file;
  size_t gridLength = 0;

  file.open(tFile.c_str());

  if (!file.is_open()) {
    std::cout << "Error cannot read file: " << tFile << std::endl;
    return -1;
  }

  // Read the number of points in one dimension, then read the points. Do this for every dimension
  // after.
  for (size_t i = 0; i < numAssests; i++) {
    file >> gridLength;
    discreteCoordinates[i] = std::vector<double>(gridLength, 0);

    for (size_t j = 0; j < gridLength; j++) {
      file >> discreteCoordinates[i][j];
    }
  }

  file.close();

  return 0;
}

int readStretchingData(std::string tFile, size_t numAssests,
                       SGPP::base::Stretching1D* streching1dArray) {
  std::fstream file;
  std::string stretchingType;
  double x_0, xsi;

  file.open(tFile.c_str());

  if (!file.is_open()) {
    std::cout << "Error cannot read file: " << tFile << std::endl;
    return -1;
  }

  for (size_t i = 0; i < numAssests; i++) {
    file >> stretchingType;
    file >> x_0;
    file >> xsi;
    streching1dArray[i].type.assign(stretchingType);
    streching1dArray[i].x_0 = x_0;
    streching1dArray[i].xsi = xsi;
  }

  file.close();

  return 0;
}

void testHeatEquation(size_t dim, size_t level, double bound_left, double bound_right, double a,
                      std::string initFunc, double T, double dt, std::string ODESolver,
                      double cg_eps, size_t cg_its, std::string fileStretch,
                      std::string stretchingMode) {
  size_t timesteps = (size_t)(T / dt);

  SGPP::base::DimensionBoundary* myBoundaries = new SGPP::base::DimensionBoundary[dim];

  // set the bounding box
  for (size_t i = 0; i < dim; i++) {
    myBoundaries[i].leftBoundary = bound_left;
    myBoundaries[i].rightBoundary = bound_right;
    myBoundaries[i].bDirichletLeft = true;
    myBoundaries[i].bDirichletRight = true;
  }

  SGPP::pde::HeatEquationSolverWithStretching* myHESolver =
      new SGPP::pde::HeatEquationSolverWithStretching();
  //  SGPP::BoundingBox* myBoundingBox = new SGPP::BoundingBox(dim, myBoundaries);
  SGPP::base::Stretching* myStretching;

  if (stretchingMode == "analytic") {
    SGPP::base::Stretching1D* stretching1dArray = new SGPP::base::Stretching1D[dim];
    int readStretchData = readStretchingData(fileStretch, dim, stretching1dArray);

    if (readStretchData != 0) {
      std::cout << "Analytic Stretching Data cannot be read, exiting.\n";
      return;
    }

    myStretching = new SGPP::base::Stretching(dim, myBoundaries, stretching1dArray);
    delete[] stretching1dArray;
  } else if (stretchingMode == "discrete") {
    std::vector<double>* discreteCoordinates = new std::vector<double>[dim];
    int readStretchData = readDiscreteStretchingData(fileStretch, dim, discreteCoordinates);

    if (readStretchData != 0) {
      std::cout << "Discrete Stretching Data cannot be read, exiting.\n";
      return;
    }

    myStretching = new SGPP::base::Stretching(dim, discreteCoordinates);
    delete[] discreteCoordinates;
  } else {
    std::cout << "Unsupported Stretching Mode Specified\n";
    return;
  }

  delete[] myBoundaries;

  // init Screen Object
  myHESolver->initScreen();

  // Construct a grid
  myHESolver->constructGrid(*myStretching, level);

  // init the basis functions' coefficient vector (start solution)
  SGPP::base::DataVector* alpha = new SGPP::base::DataVector(myHESolver->getNumberGridPoints());

  if (initFunc == "smooth") {
    myHESolver->initGridWithSmoothHeat(*alpha, bound_right, bound_right / DIV_SIGMA, DISTRI_FACTOR);
  } else {
    writeHelp();
  }

  // Print the initial heat function into a gnuplot file
  if (dim < 3) {
    myHESolver->printGrid(*alpha, GUNPLOT_RESOLUTION, "heatStart.gnuplot");
  }

  // set heat coefficient
  myHESolver->setHeatCoefficient(a);

  // Start solving the Heat Equation
  if (ODESolver == "ExEul") {
    myHESolver->solveExplicitEuler(timesteps, dt, cg_its, cg_eps, *alpha, true, false,
                                   std::max(timesteps / SOLUTION_FRAMES, (size_t)1));
  } else if (ODESolver == "ImEul") {
    myHESolver->solveImplicitEuler(timesteps, dt, cg_its, cg_eps, *alpha, true, false,
                                   std::max(timesteps / SOLUTION_FRAMES, (size_t)1));
  } else if (ODESolver == "CrNic") {
    myHESolver->solveCrankNicolson(timesteps, dt, cg_its, cg_eps, *alpha, CRNIC_IMEUL_STEPS);
  }

  // Print the solved Heat Equation into a gnuplot file
  if (dim < 3) {
    myHESolver->printGrid(*alpha, GUNPLOT_RESOLUTION, "heatSolved.gnuplot");
  }

  delete myHESolver;
  delete myStretching;
  delete alpha;
}

void testPoissonEquation(size_t dim, size_t level, double bound_left, double bound_right,
                         std::string initFunc, double cg_eps, size_t cg_its,
                         std::string fileStretch, std::string stretchingMode) {
  std::cout << "Poisson Equation for Stretching is not implemented yet, please use the "
               "non-stretched version"
            << std::endl;
}

int main(int argc, char* argv[]) {
  std::string option;

  if (argc == 1) {
    //    writeHelp();
    return 0;
  }

  option.assign(argv[1]);

  if (option == "HeatEquation") {
    if (argc != 15) {
      std::cout << "argc:" << argc << std::endl;
      writeHelp();
      return 0;
    }

    size_t dim;
    size_t level;
    double bound_left;
    double bound_right;
    double a;
    std::string initFunc;
    double T;
    double dt;
    std::string ODESolver;
    std::string fileStretch;
    std::string stretchingMode;
    double cg_eps;
    size_t cg_its;

    dim = atoi(argv[2]);
    level = atoi(argv[3]);
    bound_left = atof(argv[4]);
    bound_right = atof(argv[5]);
    a = atof(argv[6]);
    initFunc.assign(argv[7]);
    T = atof(argv[8]);
    dt = atof(argv[9]);
    ODESolver.assign(argv[10]);
    cg_eps = atof(argv[11]);
    cg_its = atoi(argv[12]);
    stretchingMode.assign(argv[13]);
    fileStretch.assign(argv[14]);

    testHeatEquation(dim, level, bound_left, bound_right, a, initFunc, T, dt, ODESolver, cg_eps,
                     cg_its, fileStretch, stretchingMode);
  } else if (option == "PoissonEquation") {
    //    if (argc != 9)
    //    {
    //      writeHelp();
    //      return 0;
    //    }
    //
    //    size_t dim;
    //    size_t level;
    //    double bound_left;
    //    double bound_right;
    //    std::string initFunc;
    //    double cg_eps;
    //    size_t cg_its;
    //
    //    dim = atoi(argv[2]);
    //    level = atoi(argv[3]);
    //    bound_left = atof(argv[4]);
    //    bound_right = atof(argv[5]);
    //    initFunc.assign(argv[6]);
    //    cg_eps = atof(argv[7]);
    //    cg_its = atoi(argv[8]);
    //    testPoissonEquation(dim, level, bound_left, bound_right, initFunc, cg_eps,
    //    cg_its,fileStretch, stretchingMode);
  } else {
    writeHelp();
  }

  return 0;
}
