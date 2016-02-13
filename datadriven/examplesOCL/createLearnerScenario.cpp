// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <iostream>
#include <string>

#include "sgpp/globaldef.hpp"

#include "sgpp/datadriven/application/LearnerScenario.hpp"

int main(int argc, char** argv) {
  int maxLevel = 1;
  //    std::string fileName = "friedman_4d_small.arff";
  std::string fileName = "friedman_4d.arff";
  sg::base::RegularGridConfiguration gridConfig;
  sg::solver::SLESolverConfiguration SLESolverConfigRefine;
  sg::solver::SLESolverConfiguration SLESolverConfigFinal;
  sg::base::AdpativityConfiguration adaptConfig;

  // setup grid
  gridConfig.dim_ = 0;  // dim is inferred from the data
  gridConfig.level_ = maxLevel;
  gridConfig.type_ = SGPP::base::GridType::Linear;

  // Set Adaptivity
  adaptConfig.maxLevelType_ = false;
  adaptConfig.noPoints_ = 80;
  adaptConfig.numRefinements_ = 0;
  adaptConfig.percent_ = 200.0;
  adaptConfig.threshold_ = 0.0;

  // Set solver during refinement
  SLESolverConfigRefine.eps_ = 0;
  SLESolverConfigRefine.maxIterations_ = 5;
  SLESolverConfigRefine.threshold_ = -1.0;
  SLESolverConfigRefine.type_ = SGPP::solver::SLESolverType::CG;

  // Set solver for final step
  SLESolverConfigFinal.eps_ = 0;
  SLESolverConfigFinal.maxIterations_ = 20;
  SLESolverConfigFinal.threshold_ = -1.0;
  SLESolverConfigFinal.type_ = SGPP::solver::SLESolverType::CG;

  double lambda = 0.000001;

  SGPP::datadriven::LearnerScenario scenario(fileName, lambda, gridConfig,
                                             SLESolverConfigRefine,
                                             SLESolverConfigFinal, adaptConfig);

  scenario.writeToFile("learnerSimple.scenario");

  SGPP::datadriven::LearnerScenario readScenario("learnerSimple.scenario");

  readScenario.writeToFile("test.scenario");
}
