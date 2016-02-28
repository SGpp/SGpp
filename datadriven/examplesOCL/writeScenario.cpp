// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <iostream>
#include <string>

#include "sgpp/base/grid/Grid.hpp"
#include "sgpp/solver/TypesSolver.hpp"

#include "sgpp/base/opencl/OCLManagerMultiPlatform.hpp"
#include "sgpp/base/opencl/OCLOperationConfiguration.hpp"
#include "sgpp/datadriven/application/LearnerScenario.hpp"

int main(int argc, char **argv) {
  int maxLevel = 10;
  //    std::string fileName = "friedman_4d_small.arff";
  std::string fileName = "DR5_train.arff";
  SGPP::base::RegularGridConfiguration gridConfig;
  SGPP::solver::SLESolverConfiguration SLESolverConfigRefine;
  SGPP::solver::SLESolverConfiguration SLESolverConfigFinal;
  SGPP::base::AdpativityConfiguration adaptConfig;

  // setup grid
  gridConfig.dim_ = 0;  // dim is inferred from the data
  gridConfig.level_ = maxLevel;
  gridConfig.type_ = SGPP::base::GridType::ModLinear;

  // Set Adaptivity
  adaptConfig.maxLevelType_ = false;
  adaptConfig.noPoints_ = 80;
  adaptConfig.numRefinements_ = 0;
  adaptConfig.percent_ = 200.0;
  adaptConfig.threshold_ = 0.0;

  // Set solver during refinement
  SLESolverConfigRefine.eps_ = 0;
  SLESolverConfigRefine.maxIterations_ = 1;
  SLESolverConfigRefine.threshold_ = -1.0;
  SLESolverConfigRefine.type_ = SGPP::solver::SLESolverType::CG;

  // Set solver for final step
  SLESolverConfigFinal.eps_ = 0;
  SLESolverConfigFinal.maxIterations_ = 1;
  SLESolverConfigFinal.threshold_ = -1.0;
  SLESolverConfigFinal.type_ = SGPP::solver::SLESolverType::CG;

  double lambda = 0.000001;

  SGPP::datadriven::LearnerScenario scenario(fileName, lambda, gridConfig, SLESolverConfigRefine,
                                             SLESolverConfigFinal, adaptConfig);

  //  scenario.writeToFile("DR5_modlinear.scenario");
  scenario.serialize("DR5_modlinear.scenario");

  std::cout << "done" << std::endl;

  return 0;
}
