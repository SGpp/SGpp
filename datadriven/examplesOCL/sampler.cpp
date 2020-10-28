// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/opencl/OCLOperationConfiguration.hpp>
#include <sgpp/datadriven/application/MetaLearner.hpp>
#include <sgpp/datadriven/operation/hash/DatadrivenOperationCommon.hpp>

#include <iostream>
#include <string>

int main(int argc, char** argv) {
  //  int maxLevel = 9;
  int maxLevel = 12;

  // std::string fileName = "debugging.arff";
  //  std::string fileName = "DR5_train_larger.arff";
  std::string fileName = "friedman2_4d_1000000.arff";
  //  std::string fileName = "friedman_4d_large.arff";
  //  std::string fileName = "friedman2_90000.arff";
  //  std::string fileName = "bigger.arff";

  sgpp::base::RegularGridConfiguration gridConfig;
  sgpp::solver::SLESolverConfiguration SLESolverConfigRefine;
  sgpp::solver::SLESolverConfiguration SLESolverConfigFinal;
  sgpp::base::AdaptivityConfiguration adaptivityConfig;

  // setup grid
  gridConfig.dim_ = 0;  // dim is inferred from the data
  gridConfig.level_ = maxLevel;
  gridConfig.type_ = sgpp::base::GridType::Linear;

  // Set Adaptivity
  adaptivityConfig.maxLevelType_ = false;
  adaptivityConfig.numRefinementPoints_ = 80;
  adaptivityConfig.numRefinements_ = 0;
  adaptivityConfig.percent_ = 200.0;
  adaptivityConfig.refinementThreshold_ = 0.0;

  // Set solver during refinement
  SLESolverConfigRefine.eps_ = 0;
  SLESolverConfigRefine.maxIterations_ = 10;
  SLESolverConfigRefine.threshold_ = -1.0;
  SLESolverConfigRefine.type_ = sgpp::solver::SLESolverType::CG;

  // Set solver for final step
  SLESolverConfigFinal.eps_ = 0;
  SLESolverConfigFinal.maxIterations_ = 3;
  SLESolverConfigFinal.threshold_ = -1.0;
  SLESolverConfigFinal.type_ = sgpp::solver::SLESolverType::CG;

  std::string metaInformation =
      "refine: " + std::to_string(adaptivityConfig.numRefinements_) + " points: " +
      std::to_string(adaptivityConfig.numRefinementPoints_) + " iterations: " +
      std::to_string(SLESolverConfigRefine.maxIterations_);

  double lambda = 0.000001;

  bool verbose = true;
  sgpp::datadriven::MetaLearner learner(gridConfig, SLESolverConfigRefine, SLESolverConfigFinal,
                                        adaptivityConfig, lambda, verbose);

  // learner.learn(kernelType, fileName);
  // learner.learnReference(fileName);

  // buggy are:
  // subspace simple - 0
  // subspacelinear combined - 60
  // streaming default - 1600 (13 without avx)
  // streaming ocl - 13

  //    sgpp::base::OCLOperationConfiguration parameters("tunedParameters.cfg");
  sgpp::base::OCLOperationConfiguration parameters("vgpu2_whole.cfg");

  sgpp::datadriven::OperationMultipleEvalConfiguration configuration(
      sgpp::datadriven::OperationMultipleEvalType::STREAMING,
      sgpp::datadriven::OperationMultipleEvalSubType::OCLMP, parameters);

  learner.learn(configuration, fileName);
  // learner.learnReference(fileName);

  // learner.learnAndTest(fileName, testFileName,
  // isBinaryClassificationProblem);
  // learner.learnAndCompare(configuration, fileName, 5);

  // learner.writeStatisticsFile("statistics.csv", "test");

  return EXIT_SUCCESS;
}
