// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/datadriven/application/MetaLearner.hpp>
#include <sgpp/datadriven/operation/hash/DatadrivenOperationCommon.hpp>

#include <iostream>
#include <string>

int main(int argc, char** argv) {
  int baseLevel = 7;

  // std::string fileName = "debugging.arff";
  std::string fileName = "DR5_train.arff";
  //  std::string fileName = "friedman2_90000.arff";
  //  std::string fileName = "bigger.arff";

  sgpp::base::RegularGridConfiguration gridConfig;
  sgpp::solver::SLESolverConfiguration SLESolverConfigRefine;
  sgpp::solver::SLESolverConfiguration SLESolverConfigFinal;
  sgpp::base::AdaptivityConfiguration adaptivityConfig;

  // setup grid
  gridConfig.dim_ = 0;  // dim is inferred from the data
  gridConfig.level_ = baseLevel;
  gridConfig.type_ = sgpp::base::GridType::Linear;

  // setup adaptivity
  adaptivityConfig.maxLevelType_ = false;
  adaptivityConfig.numRefinementPoints_ = 80;
  adaptivityConfig.numRefinements_ = 0;
  adaptivityConfig.percent_ = 200.0;
  adaptivityConfig.refinementThreshold_ = 0.0;

  // setup solver during refinement
  SLESolverConfigRefine.eps_ = 0;
  SLESolverConfigRefine.maxIterations_ = 5;
  SLESolverConfigRefine.threshold_ = -1.0;
  SLESolverConfigRefine.type_ = sgpp::solver::SLESolverType::CG;

  // setup solver for final step
  SLESolverConfigFinal.eps_ = 0;
  SLESolverConfigFinal.maxIterations_ = 5;
  SLESolverConfigFinal.threshold_ = -1.0;
  SLESolverConfigFinal.type_ = sgpp::solver::SLESolverType::CG;

  double lambda = 0.000001;

  bool verbose = true;
  sgpp::datadriven::MetaLearner learner(gridConfig, SLESolverConfigRefine, SLESolverConfigFinal,
                                        adaptivityConfig, lambda, verbose);

  // Configuration for intrisics-based streaming implementation
  sgpp::datadriven::OperationMultipleEvalConfiguration configuration(
      sgpp::datadriven::OperationMultipleEvalType::STREAMING,
      sgpp::datadriven::OperationMultipleEvalSubType::DEFAULT);

  // Configuration for subspace-based evaluation (with subspace-skipping)
  // sgpp::datadriven::OperationMultipleEvalConfiguration configuration(
  // sgpp::datadriven::OperationMultipleEvalType::SUBSPACELINEAR,
  // sgpp::datadriven::OperationMultipleEvalSubType::COMBINED);

  // only execute learning (no comparisons or tests, for performance
  // measurements)
  learner.learn(configuration, fileName);

  // execute learning with the specified configuration and use the
  // implementation from base as comparison
  // result grids are compared by sampling the domain (again with a grid) and
  // comparing the evaluated values
  //    learner.learnAndCompare(configuration, fileName, 8);

  return EXIT_SUCCESS;
}
