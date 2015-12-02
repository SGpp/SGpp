// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#define GUNPLOT_RESOLUTION 51
#define SOLUTION_FRAMES 100

#include <cstdlib>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <string>

#include <sgpp_base.hpp>
#include <sgpp_pde.hpp>
#include <sgpp_finance.hpp>
#include <sgpp_parallel.hpp>
#include <sgpp_solver.hpp>
#include <sgpp_datadriven.hpp>

/**
 * Calls the writeHelp method in the BlackScholesSolver Object
 * after creating a screen.
 */
void writeHelp() {
  SGPP::pde::LaserHeatEquationSolver2D* myHESolver = new SGPP::pde::LaserHeatEquationSolver2D(1.0, 0.4, 5, 0.1, 0.001, 4.0);

  myHESolver->initScreen();

  delete myHESolver;

  std::stringstream mySStream;

  mySStream << "Some instructions for the use of 2D Laser Heat Equation Solver:" << std::endl;
  mySStream << "---------------------------------------------------------------" << std::endl << std::endl;
  mySStream << "Parameters are:" << std::endl;
  mySStream << "     level: regular starting level" << std::endl;
  mySStream << "     maxLevel: max. refinement level" << std::endl;
  mySStream << "     initialRefines: number of initial refinements" << std::endl;
  mySStream << "     refine_threshold: refinement threshold" << std::endl;
  mySStream << "     coarsen_threshold: coarsening threshold" << std::endl;
  mySStream << "     v: laser beam velocity" << std::endl;
  mySStream << "     a: material coefficient" << std::endl;
  mySStream << "     heat: initial heating" << std::endl;
  mySStream << "     heat_sigma: heat expansion" << std::endl;
  mySStream << "     T: simulation time" << std::endl;
  mySStream << "     dt: timestepsize" << std::endl;
  mySStream << "     cg_its: max. CG iterations" << std::endl;
  mySStream << "     cg_eps: CG epsilon" << std::endl;
  mySStream << "     animation: 0 no plots, 1 generate plots for animation" << std::endl << std::endl << std::endl;
  mySStream << "Example:" << std::endl;
  mySStream << "     LaserHESolver2D 5 8 5 0.01 0.001 1.0 10.0 4.0 0.04 2.0 0.001 400 0.00001 0" << std::endl;
  mySStream << std::endl << std::endl;

  std::cout << mySStream.str() << std::endl;
}

void testLaserHeatEquation( size_t level, size_t maxLevel, size_t initialRefines, double refine_threshold,
                            double coarsen_threshold, double v, double a, double heat, double heat_sigma,
                            double T, double dt, size_t cg_its, double cg_eps, size_t animation) {
  size_t timesteps = (size_t)(T / dt);

  SGPP::base::DimensionBoundary* myBoundaries = new SGPP::base::DimensionBoundary[2];

  // set the bounding box
  for (size_t i = 0; i < 2; i++) {
    myBoundaries[i].leftBoundary = 0.0;
    myBoundaries[i].rightBoundary = 1.0;
    myBoundaries[i].bDirichletLeft = true;
    myBoundaries[i].bDirichletRight = true;
  }

  SGPP::pde::LaserHeatEquationSolver2D* myHESolver = new SGPP::pde::LaserHeatEquationSolver2D(v, heat_sigma, maxLevel, refine_threshold, coarsen_threshold, heat);
  SGPP::base::BoundingBox* myBoundingBox = new SGPP::base::BoundingBox(2, myBoundaries);
  delete[] myBoundaries;

  // init Screen Object
  myHESolver->initScreen();

  // Construct a grid
  myHESolver->constructGrid(*myBoundingBox, level);

  // init the basis functions' coefficient vector (start solution)
  SGPP::base::DataVector* alpha = new SGPP::base::DataVector(myHESolver->getNumberGridPoints());
  myHESolver->refineInitialGridWithLaserHeat(*alpha, initialRefines);

  // Print the initial heat function into a gnuplot file
  myHESolver->printGrid(*alpha, GUNPLOT_RESOLUTION, "laserheatStart.gnuplot");

  // set heat coefficient
  myHESolver->setHeatCoefficient(a);

  // Start solving the Heat Equation
  if (animation == 0) {
    myHESolver->solveImplicitEuler(timesteps, dt, cg_its, cg_eps, *alpha, true, false, std::max(timesteps / SOLUTION_FRAMES, (size_t)1));
  } else {
    myHESolver->solveImplicitEuler(timesteps, dt, cg_its, cg_eps, *alpha, true, true, std::max(timesteps / SOLUTION_FRAMES, (size_t)1));
  }

  // Print the solved Heat Equation into a gnuplot file
  myHESolver->printGrid(*alpha, GUNPLOT_RESOLUTION, "laserheatSolved.gnuplot");

  delete myHESolver;
  delete myBoundingBox;
  delete alpha;
}

int main(int argc, char* argv[]) {
  if (argc != 15) {
    writeHelp();
    return 0;
  }

  size_t level;
  size_t maxLevel;
  size_t initialRefines;
  double refine_threshold;
  double coarsen_threshold;

  double v;
  double a;
  double heat;
  double heat_sigma;

  double T;
  double dt;
  size_t cg_its;
  double cg_eps;

  size_t animation;

  level = atoi(argv[1]);
  maxLevel = atoi(argv[2]);
  initialRefines = atoi(argv[3]);
  refine_threshold = atof(argv[4]);
  coarsen_threshold = atof(argv[5]);

  v = atof(argv[6]);
  a = atof(argv[7]);
  heat = atof(argv[8]);
  heat_sigma = atof(argv[9]);

  T = atof(argv[10]);
  dt = atof(argv[11]);
  cg_its = atoi(argv[12]);
  cg_eps = atof(argv[13]);

  animation = atoi(argv[14]);

  testLaserHeatEquation(level, maxLevel, initialRefines, refine_threshold, coarsen_threshold, v, a, heat, heat_sigma, T, dt, cg_its, cg_eps, animation);

  return 0;
}