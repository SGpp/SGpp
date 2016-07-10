// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#define CRNIC_IMEUL_STEPS 3
#define GUNPLOT_RESOLUTION 51
#define SOLUTION_FRAMES 100

#define DIV_SIGMA 4.0
#define DISTRI_FACTOR 5.0

#define EXP_FACTOR 2.0

#define NUMEVALPOINTS 11

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wcast-qual"
#include <mpi.h>
#pragma GCC diagnostic pop

#include <sgpp_base.hpp>
#include <sgpp_pde.hpp>
#include <sgpp_finance.hpp>
#include <sgpp_parallel.hpp>
#include <sgpp_solver.hpp>
#include <sgpp_datadriven.hpp>
#include <sgpp_mpi.hpp>

#include <cstdlib>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <string>
#include <iomanip>
#include <vector>

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
 * reads function values from a file
 *
 * @param values DataVector into which the values will be stored
 * @param tFile file from which the values are read
 */
int readCuboidValues(sgpp::base::DataVector& values, std::string tFile) {
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
 * Writes a DataMatrix into a file
 *
 * @param data the DataMatrix that should be written into a file
 * @param tFile the file into which the data is written
 *
 * @return error code
 */
int writeDataMatrix(sgpp::base::DataMatrix& data, std::string tFile) {
  std::ofstream file;
  file.open(tFile.c_str());

  if (!file.is_open()) {
    std::cout << "Error cannot write file: " << tFile << std::endl;
    return -1;
  }

  for (size_t i = 0; i < data.getNrows(); i++) {
    for (size_t j = 0; j < data.getNcols(); j++) {
      file << std::scientific << std::setprecision(16) << data.get(i, j) << " ";
    }

    file << std::endl;
  }

  file.close();

  return 0;
}

/**
 * Writes a DataVector into a file
 *
 * @param data the DataVector that should be written into a file
 * @param tFile the file into which the data is written
 *
 * @return error code
 */
int writeDataVector(sgpp::base::DataVector& data, std::string tFile) {
  std::ofstream file;
  file.open(tFile.c_str());

  if (!file.is_open()) {
    std::cout << "Error cannot write file: " << tFile << std::endl;
    return -1;
  }

  for (size_t i = 0; i < data.getSize(); i++) {
    file << std::scientific << std::setprecision(16) << data.get(i) << " " << std::endl;
  }

  file.close();

  return 0;
}

/**
 * Calls the writeHelp method in the BlackScholesSolver Object
 * after creating a screen.
 */
void writeHelp() {
  sgpp::pde::HeatEquationSolver* myHESolver = new sgpp::pde::HeatEquationSolver();

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
  mySStream << "  start_level: number of start levels of the Sparse Grid" << std::endl;
  mySStream << "  end_level: number of max. levels of the Sparse Grid" << std::endl;
  mySStream << "  left_bound: x_i of left boundary" << std::endl;
  mySStream << "  right_bound: x_i of right boundary" << std::endl;
  mySStream << "  a: thermal diffusivity" << std::endl;
  mySStream << "  initHeat: initial heat distribution" << std::endl;
  mySStream << "  T: time to solve" << std::endl;
  mySStream << "  dT: timestep size" << std::endl;
  mySStream << "  Solver: the solver to use: ExEul, ImEul, CrNic" << std::endl;
  mySStream << "  CGEpsilon: Epsilon used in CG" << std::endl;
  mySStream << "  CGIterations: Maxmimum number of iterations used in CG mehtod" << std::endl;
  mySStream << std::endl;
  mySStream << "Example:" << std::endl;
  mySStream << "HESolverMPI HeatEquation 3 2 5 0.0 3.0 1.0 smooth 1.0 0.1 ImEul 0.00001 400"
            << std::endl;
  mySStream << std::endl;
  mySStream << "Remark: This test generates following files (gnuplot):" << std::endl;
  mySStream << "  heatStart.gnuplot: the start condition" << std::endl;
  mySStream << "  heatSolved.gnuplot: the numerical solution" << std::endl;
  mySStream << std::endl << std::endl;

  mySStream << "PoissonEquation" << std::endl << "------" << std::endl;
  mySStream << "the following options must be specified:" << std::endl;
  mySStream << "  dim: the number of dimensions of Sparse Grid" << std::endl;
  mySStream << "  start_level: number of start levels of the Sparse Grid" << std::endl;
  mySStream << "  end_level: number of max levels of the Sparse Grid" << std::endl;
  mySStream << "  left_bound: x_i of left boundary" << std::endl;
  mySStream << "  right_bound: x_i of right boundary" << std::endl;
  mySStream << "  initHeat: initial heat distribution" << std::endl;
  mySStream << "  CGEpsilon: Epsilon used in CG" << std::endl;
  mySStream << "  CGIterations: Maxmimum number of iterations used in CG mehtod" << std::endl;
  mySStream << std::endl;
  mySStream << "Example:" << std::endl;
  mySStream << "HESolverMPI PoissonEquation 3 2 5 0.0 3.0 smooth 0.00001 400" << std::endl;
  mySStream << std::endl;
  mySStream << "Remark: This test generates following files (gnuplot):" << std::endl;
  mySStream << "  poissonStart.gnuplot: the start condition" << std::endl;
  mySStream << "  poissonSolved.gnuplot: the numerical solution" << std::endl;
  mySStream << std::endl << std::endl;

  mySStream << "PoissonEquationAdapt" << std::endl << "------" << std::endl;
  mySStream << "the following options must be specified:" << std::endl;
  mySStream << "  dim: the number of dimensions of Sparse Grid" << std::endl;
  mySStream << "  start_level: number of start-levels of the Sparse Grid" << std::endl;
  mySStream << "  refinemode: method used for initial refinement: classic, maxLevel" << std::endl;
  mySStream << "  max_ref_level: maxLevel used" << std::endl;
  mySStream << "  num_refines: Number of initial refinements" << std::endl;
  mySStream << "  refine_thres: Refinement threshold" << std::endl;
  mySStream << "  left_bound: x_i of left boundary" << std::endl;
  mySStream << "  right_bound: x_i of right boundary" << std::endl;
  mySStream << "  initHeat: initial heat distribution" << std::endl;
  mySStream << "  CGEpsilon: Epsilon used in CG" << std::endl;
  mySStream << "  CGIterations: Maxmimum number of iterations used in CG mehtod" << std::endl;
  mySStream << std::endl;
  mySStream << "Example:" << std::endl;
  mySStream << "HESolverMPI PoissonEquationAdapt 3 2 maxLevel 5 5 1e-7 0.0 3.0 smooth 0.00001 400"
            << std::endl;
  mySStream << std::endl;
  mySStream << "Remark: This test generates following files (gnuplot):" << std::endl;
  mySStream << "  poissonStartAdapt.gnuplot: the start condition" << std::endl;
  mySStream << "  poissonSolvedAdapt.gnuplot: the numerical solution" << std::endl;
  mySStream << std::endl << std::endl;

  std::cout << mySStream.str() << std::endl;
}

void testHeatEquation(size_t dim, size_t start_level, size_t end_level, double bound_left,
                      double bound_right, double a, std::string initFunc, double T, double dt,
                      std::string ODESolver, double cg_eps, size_t cg_its) {
  size_t timesteps = (size_t)(T / dt);
  sgpp::parallel::HeatEquationSolverMPI* myHESolver = new sgpp::parallel::HeatEquationSolverMPI();
  sgpp::base::DataVector* alpha = NULL;
  sgpp::base::DataMatrix EvalPoints(1, dim);
  std::string tFileEvalCuboid = "EvalPointsHeatEquationMPI.data";
  std::string tFileEvalCuboidValues = "EvalValuesHeatEquationMPI.data";
  size_t evalPoints = NUMEVALPOINTS;
  std::vector<sgpp::base::DataVector> results;

  for (size_t l = start_level; l <= end_level; l++) {
    if (sgpp::parallel::myGlobalMPIComm->getMyRank() == 0) {
      std::vector<sgpp::base::BoundingBox1D> myBoundaries(dim, sgpp::base::BoundingBox1D());

      // set the bounding box
      for (size_t i = 0; i < dim; i++) {
        myBoundaries[i].leftBoundary = bound_left;
        myBoundaries[i].rightBoundary = bound_right;
        myBoundaries[i].bDirichletLeft = true;
        myBoundaries[i].bDirichletRight = true;
      }

      sgpp::base::BoundingBox* myBoundingBox = new sgpp::base::BoundingBox(myBoundaries);

      // in first iteration -> calculate the evaluation points
      if (l == start_level) {
        sgpp::base::EvalCuboidGenerator* myEvalCuboidGen = new sgpp::base::EvalCuboidGenerator();
        myEvalCuboidGen->getEvaluationCuboid(EvalPoints, *myBoundingBox, evalPoints);
        writeDataMatrix(EvalPoints, tFileEvalCuboid);
        delete myEvalCuboidGen;
      }

      // init Screen Object
      myHESolver->initScreen();

      // Construct a grid
      myHESolver->constructGrid(*myBoundingBox, static_cast<int>(l));

      // init the basis functions' coefficient vector (start solution)
      alpha = new sgpp::base::DataVector(myHESolver->getNumberGridPoints());

      if (initFunc == "smooth") {
        myHESolver->initGridWithSmoothHeat(*alpha, bound_right, bound_right / DIV_SIGMA,
                                           DISTRI_FACTOR);
      } else {
        writeHelp();
        sgpp::parallel::myGlobalMPIComm->Abort();
      }

      delete myBoundingBox;
    }

    // Communicate grid
    if (sgpp::parallel::myGlobalMPIComm->getMyRank() == 0) {
      std::string serialized_grid = myHESolver->getGrid();

      sgpp::parallel::myGlobalMPIComm->broadcastGrid(serialized_grid);
    } else {
      // Now receive the grid
      std::string serialized_grid = "";

      sgpp::parallel::myGlobalMPIComm->receiveGrid(serialized_grid);
      myHESolver->setGrid(serialized_grid);

      alpha = new sgpp::base::DataVector(myHESolver->getNumberGridPoints());
    }

    // Communicate coefficients
    sgpp::parallel::myGlobalMPIComm->broadcastGridCoefficientsFromRank0(*alpha);

    // Print initial grid only on rank 0
    if (sgpp::parallel::myGlobalMPIComm->getMyRank() == 0) {
      // Print the initial heat function into a gnuplot file
      if (dim < 3) {
        myHESolver->printGrid(*alpha, GUNPLOT_RESOLUTION, "heatStartMPI.gnuplot");
      }
    }

    // set heat coefficient
    myHESolver->setHeatCoefficient(a);

    // Start solving the Heat Equation
    if (ODESolver == "ExEul") {
      myHESolver->solveExplicitEuler(timesteps, dt, cg_its, cg_eps, *alpha, true, false);
    } else if (ODESolver == "ImEul") {
      myHESolver->solveImplicitEuler(timesteps, dt, cg_its, cg_eps, *alpha, true, false);
    } else if (ODESolver == "CrNic") {
      myHESolver->solveCrankNicolson(timesteps, dt, cg_its, cg_eps, *alpha, CRNIC_IMEUL_STEPS);
    }

    // print solved grid only on rank 0
    if (sgpp::parallel::myGlobalMPIComm->getMyRank() == 0) {
      // Print the solved Heat Equation into a gnuplot file
      if (dim < 3) {
        myHESolver->printGrid(*alpha, GUNPLOT_RESOLUTION, "heatSolvedMPI.gnuplot");
      }

      // Calculate Norms
      // Evaluate Cuboid
      sgpp::base::DataVector PoisEvals(EvalPoints.getNrows());
      myHESolver->evaluateCuboid(*alpha, PoisEvals, EvalPoints);
      results.push_back(PoisEvals);

      // write solution in a additional file
      std::stringstream level_string;
      level_string << l;
      writeDataVector(PoisEvals, tFileEvalCuboidValues + ".level_" + level_string.str());
      writeDataVector(PoisEvals, tFileEvalCuboidValues);

      if (l > start_level) {
        std::cout << "====================================================================="
                  << std::endl;
        std::cout << "====================================================================="
                  << std::endl
                  << std::endl;
        std::cout << "Calculating norms of relative errors to a grid" << std::endl;
        std::cout << "with " << l << " levels and testing-coboid" << std::endl;

        double oldMaxNorm = 0.0;
        double oldTwoNorm = 0.0;

        // Calculate relative errors and some norms
        for (size_t j = 0; j < l - start_level; j++) {
          sgpp::base::DataVector maxLevel(results[l - start_level]);
          sgpp::base::DataVector relError(results[j]);
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
          std::cout << "Level " << j + start_level << ": max-norm(rel-error)=" << maxNorm
                    << "; two-norm(rel-error)=" << l2Norm
                    << "; rate max-norm: " << log(oldMaxNorm / maxNorm)
                    << "; rate two-norm: " << log(oldTwoNorm / l2Norm) << std::endl;

          oldMaxNorm = maxNorm;
          oldTwoNorm = l2Norm;
        }
      }

      std::cout << std::endl << std::endl;
    }

    myHESolver->deleteGrid();

    delete alpha;
    alpha = NULL;
  }

  delete myHESolver;
}

void testPoissonEquation(size_t dim, size_t start_level, size_t end_level, double bound_left,
                         double bound_right, std::string initFunc, double cg_eps, size_t cg_its) {
  sgpp::parallel::PoissonEquationSolverMPI* myPoisSolver =
      new sgpp::parallel::PoissonEquationSolverMPI();
  sgpp::base::DataVector* alpha = NULL;
  sgpp::base::DataMatrix EvalPoints(1, dim);
  std::string tFileEvalCuboid = "EvalPointsPoissonMPI.data";
  std::string tFileEvalCuboidValues = "EvalValuesPoissonMPI.data";
  size_t evalPoints = NUMEVALPOINTS;
  std::vector<sgpp::base::DataVector> results;

  for (size_t l = start_level; l <= end_level; l++) {
    if (sgpp::parallel::myGlobalMPIComm->getMyRank() == 0) {
      std::vector<sgpp::base::BoundingBox1D> myBoundaries(dim, sgpp::base::BoundingBox1D());

      // set the bounding box
      for (size_t i = 0; i < dim; i++) {
        myBoundaries[i].leftBoundary = bound_left;
        myBoundaries[i].rightBoundary = bound_right;
        myBoundaries[i].bDirichletLeft = true;
        myBoundaries[i].bDirichletRight = true;
      }

      sgpp::base::BoundingBox* myBoundingBox = new sgpp::base::BoundingBox(myBoundaries);

#ifndef __MIC__

      // in first iteration -> calculate the evaluation points
      if (l == start_level) {
        sgpp::base::EvalCuboidGenerator* myEvalCuboidGen = new sgpp::base::EvalCuboidGenerator();
        myEvalCuboidGen->getEvaluationCuboid(EvalPoints, *myBoundingBox, evalPoints);
        writeDataMatrix(EvalPoints, tFileEvalCuboid);
        delete myEvalCuboidGen;
      }

#endif
      // init Screen Object
      myPoisSolver->initScreen();

      // Construct a grid
      myPoisSolver->constructGrid(*myBoundingBox, static_cast<int>(l));

      // init the basis functions' coefficient vector (start solution)
      alpha = new sgpp::base::DataVector(myPoisSolver->getNumberGridPoints());

      if (initFunc == "smooth") {
        myPoisSolver->initGridWithSmoothHeat(*alpha, bound_right, bound_right / DIV_SIGMA,
                                             DISTRI_FACTOR);
      } else if (initFunc == "exp") {
        myPoisSolver->initGridWithExpHeat(*alpha, EXP_FACTOR);
      } else {
        writeHelp();
        sgpp::parallel::myGlobalMPIComm->Abort();
        delete alpha, delete myPoisSolver;
        delete myBoundingBox;
        return;
      }

      delete myBoundingBox;
    }

    // Communicate grid
    if (sgpp::parallel::myGlobalMPIComm->getMyRank() == 0) {
      std::string serialized_grid = myPoisSolver->getGrid();

      sgpp::parallel::myGlobalMPIComm->broadcastGrid(serialized_grid);
    } else {
      // Now receive the grid
      std::string serialized_grid = "";

      sgpp::parallel::myGlobalMPIComm->receiveGrid(serialized_grid);
      myPoisSolver->setGrid(serialized_grid);

      alpha = new sgpp::base::DataVector(myPoisSolver->getNumberGridPoints());
    }

    // Communicate coefficients
    sgpp::parallel::myGlobalMPIComm->broadcastGridCoefficientsFromRank0(*alpha);

    if (sgpp::parallel::myGlobalMPIComm->getMyRank() == 0) {
      // Print the initial heat function into a gnuplot file
      if (dim < 3) {
        myPoisSolver->printGrid(*alpha, GUNPLOT_RESOLUTION, "poissonStartMPI.gnuplot");
      }
    }

    // solve Poisson Equation
    myPoisSolver->solvePDE(*alpha, *alpha, cg_its, cg_eps, true);

#ifndef __MIC__

    if (sgpp::parallel::myGlobalMPIComm->getMyRank() == 0) {
      // Print the solved Heat Equation into a gnuplot file
      if (dim < 3) {
        myPoisSolver->printGrid(*alpha, GUNPLOT_RESOLUTION, "poissonSolvedMPI.gnuplot");
      }

      // Calculate Norms
      // Evaluate Cuboid
      sgpp::base::DataVector PoisEvals(EvalPoints.getNrows());
      myPoisSolver->evaluateCuboid(*alpha, PoisEvals, EvalPoints);
      results.push_back(PoisEvals);

      // write solution in a additional file
      std::stringstream level_string;
      level_string << l;
      writeDataVector(PoisEvals, tFileEvalCuboidValues + ".level_" + level_string.str());
      writeDataVector(PoisEvals, tFileEvalCuboidValues);

      if (l > start_level) {
        std::cout << "====================================================================="
                  << std::endl;
        std::cout << "====================================================================="
                  << std::endl
                  << std::endl;
        std::cout << "Calculating norms of relative errors to a grid" << std::endl;
        std::cout << "with " << l << " levels and testing-coboid" << std::endl;

        double oldMaxNorm = 0.0;
        double oldTwoNorm = 0.0;

        // Calculate relative errors and some norms
        for (size_t j = 0; j < l - start_level; j++) {
          sgpp::base::DataVector maxLevel(results[l - start_level]);
          sgpp::base::DataVector relError(results[j]);
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
          std::cout << "Level " << j + start_level << ": max-norm(rel-error)=" << maxNorm
                    << "; two-norm(rel-error)=" << l2Norm
                    << "; rate max-norm: " << log(oldMaxNorm / maxNorm)
                    << "; rate two-norm: " << log(oldTwoNorm / l2Norm) << std::endl;

          oldMaxNorm = maxNorm;
          oldTwoNorm = l2Norm;
        }
      }

      std::cout << std::endl << std::endl;
    }

#endif

    myPoisSolver->deleteGrid();

    delete alpha;
    alpha = NULL;
  }

  delete myPoisSolver;
}

void testPoissonEquationAdapt(size_t dim, size_t start_level, std::string refine,
                              size_t max_ref_level, size_t num_refines, double refine_thres,
                              double bound_left, double bound_right, std::string initFunc,
                              double cg_eps, size_t cg_its) {
  sgpp::parallel::PoissonEquationSolverMPI* myPoisSolver =
      new sgpp::parallel::PoissonEquationSolverMPI();
  sgpp::base::DataVector* alpha = NULL;
  std::string tFileEvalCuboid = "EvalPointsPoissonMPI.data";
  std::string tFileEvalCuboidValues = "EvalValuesPoissonMPI.data";

  if (sgpp::parallel::myGlobalMPIComm->getMyRank() == 0) {
    std::vector<sgpp::base::BoundingBox1D> myBoundaries(dim, sgpp::base::BoundingBox1D());

    // set the bounding box
    for (size_t i = 0; i < dim; i++) {
      myBoundaries[i].leftBoundary = bound_left;
      myBoundaries[i].rightBoundary = bound_right;
      myBoundaries[i].bDirichletLeft = true;
      myBoundaries[i].bDirichletRight = true;
    }

    sgpp::base::BoundingBox* myBoundingBox = new sgpp::base::BoundingBox(myBoundaries);

    // init Screen Object
    myPoisSolver->initScreen();

    // Construct a grid
    myPoisSolver->constructGrid(*myBoundingBox, static_cast<int>(start_level));

    // init the basis functions' coefficient vector (start solution)
    // init the basis functions' coefficient vector (start solution)
    alpha = new sgpp::base::DataVector(myPoisSolver->getNumberGridPoints());

    if (initFunc == "smooth") {
      std::cout << "Starting Grid size: " << myPoisSolver->getNumberGridPoints() << std::endl;
      std::cout << "Starting Grid size (inner): " << myPoisSolver->getNumberInnerGridPoints()
                << std::endl
                << std::endl;

      for (size_t i = 0; i < num_refines; i++) {
        std::cout << "Refining Grid..." << std::endl;
        myPoisSolver->initGridWithSmoothHeatFullDomain(*alpha, bound_right, bound_right / DIV_SIGMA,
                                                       DISTRI_FACTOR);

        if (refine == "classic") {
          myPoisSolver->refineInitialGridSurplus(*alpha, -1, refine_thres);
        } else {
          myPoisSolver->refineInitialGridSurplusToMaxLevel(*alpha, refine_thres,
                                                           static_cast<int>(max_ref_level));
        }

        std::cout << "Refined Grid size: " << myPoisSolver->getNumberGridPoints() << std::endl;
        std::cout << "Refined Grid size (inner): " << myPoisSolver->getNumberInnerGridPoints()
                  << std::endl;
      }

      std::cout << std::endl << std::endl;

      myPoisSolver->initGridWithSmoothHeat(*alpha, bound_right, bound_right / DIV_SIGMA,
                                           DISTRI_FACTOR);
    } else if (initFunc == "exp") {
      std::cout << "Starting Grid size: " << myPoisSolver->getNumberGridPoints() << std::endl;
      std::cout << "Starting Grid size (inner): " << myPoisSolver->getNumberInnerGridPoints()
                << std::endl
                << std::endl;

      for (size_t i = 0; i < num_refines; i++) {
        std::cout << "Refining Grid..." << std::endl;
        myPoisSolver->initGridWithExpHeatFullDomain(*alpha, EXP_FACTOR);

        if (refine == "classic") {
          myPoisSolver->refineInitialGridSurplus(*alpha, -1, refine_thres);
        } else {
          myPoisSolver->refineInitialGridSurplusToMaxLevel(*alpha, refine_thres,
                                                           static_cast<int>(max_ref_level));
        }

        std::cout << "Refined Grid size: " << myPoisSolver->getNumberGridPoints() << std::endl;
        std::cout << "Refined Grid size (inner): " << myPoisSolver->getNumberInnerGridPoints()
                  << std::endl;
      }

      std::cout << std::endl << std::endl;

      myPoisSolver->initGridWithExpHeat(*alpha, EXP_FACTOR);
    } else {
      writeHelp();
      sgpp::parallel::myGlobalMPIComm->Abort();
      delete alpha, delete myPoisSolver;
      delete myBoundingBox;
      return;
    }

    delete myBoundingBox;
  }

  // Communicate grid
  if (sgpp::parallel::myGlobalMPIComm->getMyRank() == 0) {
    std::string serialized_grid = myPoisSolver->getGrid();

    sgpp::parallel::myGlobalMPIComm->broadcastGrid(serialized_grid);
  } else {
    // Now receive the grid
    std::string serialized_grid = "";

    sgpp::parallel::myGlobalMPIComm->receiveGrid(serialized_grid);
    myPoisSolver->setGrid(serialized_grid);

    alpha = new sgpp::base::DataVector(myPoisSolver->getNumberGridPoints());
  }

  // Communicate coefficients
  sgpp::parallel::myGlobalMPIComm->broadcastGridCoefficientsFromRank0(*alpha);

  if (sgpp::parallel::myGlobalMPIComm->getMyRank() == 0) {
    // Print the initial heat function into a gnuplot file
    if (dim < 3) {
      myPoisSolver->printSparseGrid(*alpha, "poissonStartMPIAdapt.gnuplot", false);
    }
  }

  // solve Poisson Equation
  myPoisSolver->solvePDE(*alpha, *alpha, cg_its, cg_eps, true);

  if (sgpp::parallel::myGlobalMPIComm->getMyRank() == 0) {
    // Print the solved Heat Equation into a gnuplot file
    if (dim < 3) {
      myPoisSolver->printSparseGrid(*alpha, "poissonSolvedMPIAdapt.gnuplot", false);
    }

    // calculate relative errors
    ////////////////////////////

    // read Evaluation cuboid
    sgpp::base::DataMatrix EvalCuboid(1, dim);
    int retCuboid = readEvalutionCuboid(EvalCuboid, tFileEvalCuboid, dim);

    // read reference values for evaluation cuboid
    sgpp::base::DataVector EvalCuboidValues(1);
    int retCuboidValues = readCuboidValues(EvalCuboidValues, tFileEvalCuboidValues);

    double maxNorm = 0.0;
    double l2Norm = 0.0;

    if (retCuboid == 0 && retCuboidValues == 0) {
      std::cout << "Calculating relative errors..." << std::endl;
      // Evaluate Cuboid
      sgpp::base::DataVector curVals(EvalCuboid.getNrows());
      myPoisSolver->evaluateCuboid(*alpha, curVals, EvalCuboid);

      sgpp::base::DataVector relError(curVals);

      // calculate relative error
      relError.sub(EvalCuboidValues);
      relError.componentwise_div(EvalCuboidValues);

      // calculate max. norm of relative error
      maxNorm = relError.maxNorm();

      // calculate two norm of relative error
      l2Norm = relError.RMSNorm();

      // Printing norms
      std::cout << "Results: max-norm(rel-error)=" << maxNorm << "; two-norm(rel-error)=" << l2Norm
                << std::endl;

      // reprint data with prefix -> can be easily grep-ed
      std::cout << std::endl << std::endl;
    } else {
      std::cout << "Couldn't open evaluation cuboid data -> skipping tests!" << std::endl
                << std::endl;
    }
  }

  delete alpha;
  delete myPoisSolver;
}

int main(int argc, char* argv[]) {
  std::string option;
  int mpi_myid;
  int mpi_ranks;

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_ranks);
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_myid);
  sgpp::parallel::myGlobalMPIComm = new sgpp::parallel::MPICommunicator(mpi_myid, mpi_ranks);

  std::streambuf* stdoutBuf = std::cout.rdbuf();
  std::ofstream dummy_out("/dev/null");

  if (mpi_myid != 0) {  // disable output for all processes but proc 0
    std::cout.rdbuf(dummy_out.rdbuf());
  }

  if (argc == 1) {
    if (mpi_myid == 0) {
      writeHelp();
    }

    sgpp::parallel::myGlobalMPIComm->Abort();
    return 0;
  }

  option.assign(argv[1]);

  if (option == "HeatEquation") {
    if (argc != 14) {
      if (mpi_myid == 0) {
        writeHelp();
      }

      sgpp::parallel::myGlobalMPIComm->Abort();
      return 0;
    }

    size_t dim;
    size_t start_level;
    size_t end_level;
    double bound_left;
    double bound_right;
    double a;
    std::string initFunc;
    double T;
    double dt;
    std::string ODESolver;
    double cg_eps;
    size_t cg_its;

    dim = atoi(argv[2]);
    start_level = atoi(argv[3]);
    end_level = atoi(argv[4]);
    bound_left = atof(argv[5]);
    bound_right = atof(argv[6]);
    a = atof(argv[7]);
    initFunc.assign(argv[8]);
    T = atof(argv[9]);
    dt = atof(argv[10]);
    ODESolver.assign(argv[11]);
    cg_eps = atof(argv[12]);
    cg_its = atoi(argv[13]);

    testHeatEquation(dim, start_level, end_level, bound_left, bound_right, a, initFunc, T, dt,
                     ODESolver, cg_eps, cg_its);
  } else if (option == "PoissonEquation") {
    if (argc != 10) {
      if (mpi_myid == 0) {
        writeHelp();
      }

      sgpp::parallel::myGlobalMPIComm->Abort();
      return 0;
    }

    size_t dim;
    size_t start_level;
    size_t end_level;
    double bound_left;
    double bound_right;
    std::string initFunc;
    double cg_eps;
    size_t cg_its;

    dim = atoi(argv[2]);
    start_level = atoi(argv[3]);
    end_level = atoi(argv[4]);
    bound_left = atof(argv[5]);
    bound_right = atof(argv[6]);
    initFunc.assign(argv[7]);
    cg_eps = atof(argv[8]);
    cg_its = atoi(argv[9]);
    testPoissonEquation(dim, start_level, end_level, bound_left, bound_right, initFunc, cg_eps,
                        cg_its);
  } else if (option == "PoissonEquationAdapt") {
    if (argc != 13) {
      if (mpi_myid == 0) {
        writeHelp();
      }

      sgpp::parallel::myGlobalMPIComm->Abort();
      return 0;
    }

    size_t dim;
    size_t start_level;
    std::string refine;
    size_t max_ref_level;
    size_t num_refines;
    double refine_thres;

    double bound_left;
    double bound_right;
    std::string initFunc;
    double cg_eps;
    size_t cg_its;

    dim = atoi(argv[2]);
    start_level = atoi(argv[3]);
    refine.assign(argv[4]);
    max_ref_level = atoi(argv[5]);
    num_refines = atoi(argv[6]);
    refine_thres = atof(argv[7]);
    bound_left = atof(argv[8]);
    bound_right = atof(argv[9]);
    initFunc.assign(argv[10]);
    cg_eps = atof(argv[11]);
    cg_its = atoi(argv[12]);

    testPoissonEquationAdapt(dim, start_level, refine, max_ref_level, num_refines, refine_thres,
                             bound_left, bound_right, initFunc, cg_eps, cg_its);
  } else {
    if (mpi_myid == 0) {
      writeHelp();
    }
  }

  delete sgpp::parallel::myGlobalMPIComm;

  if (mpi_myid != 0) {  // restore stdout buffer
    std::cout.rdbuf(stdoutBuf);
  }

  MPI_Finalize();

  return 0;
}
