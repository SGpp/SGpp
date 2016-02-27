// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <iostream>
#include <string>

#include "sgpp/globaldef.hpp"
#include "sgpp/datadriven/application/LearnerScenario.hpp"
#include "sgpp/datadriven/application/MetaLearner.hpp"
#include "sgpp/datadriven/operation/hash/simple/DatadrivenOperationCommon.hpp"
#include "sgpp/base/operation/hash/OperationMultipleEval.hpp"
#include "sgpp/datadriven/DatadrivenOpFactory.hpp"
#include "sgpp/datadriven/tools/ARFFTools.hpp"
#include "sgpp/base/opencl/OCLOperationConfiguration.hpp"

int main(int argc, char** argv) {
  //    std::string fileName = "friedman_4d_small.arff";

  //  std::string scenarioFileName("DR5_Linear.scenario");
  //  std::string fileName = "DR5_train.arff";
  //  std::string alphaReferenceFileName = "DR5_Linear_AlphaReference.vec";

  std::string scenarioFileName("friedman2_4d_Linear.scenario");
  std::string fileName = "friedman2_4d_300000.arff";
  std::string alphaReferenceFileName = "friedman2_4d_Linear_AlphaReference.vec";

  std::string parameterFileName("allDevices.cfg");
  double lambda = 0.0000001;

  SGPP::base::RegularGridConfiguration gridConfig;
  SGPP::solver::SLESolverConfiguration SLESolverConfigRefine;
  SGPP::solver::SLESolverConfiguration SLESolverConfigFinal;
  SGPP::base::AdpativityConfiguration adaptConfig;
  SGPP::datadriven::TestsetConfiguration testsetConfig;

  // setup grid
  gridConfig.dim_ = 0;  // dim is inferred from the data
  gridConfig.level_ = 10;
  gridConfig.type_ = SGPP::base::GridType::Linear;
  // dummy values
  gridConfig.boundaryLevel_ = 0;
  gridConfig.maxDegree_ = 30;

  // Set Adaptivity
  adaptConfig.maxLevelType_ = false;
  adaptConfig.noPoints_ = 80;
  adaptConfig.numRefinements_ = 0;
  adaptConfig.percent_ = 200.0;
  adaptConfig.threshold_ = 0.0;

  // Set solver during refinement
  SLESolverConfigRefine.eps_ = 0;
  SLESolverConfigRefine.maxIterations_ = 200;
  SLESolverConfigRefine.threshold_ = -1.0;
  SLESolverConfigRefine.type_ = SGPP::solver::SLESolverType::CG;

  // Set solver for final step
  SLESolverConfigFinal.eps_ = 0;
  SLESolverConfigFinal.maxIterations_ = 1;
  SLESolverConfigFinal.threshold_ = -1.0;
  SLESolverConfigFinal.type_ = SGPP::solver::SLESolverType::CG;

  // create scenario to learn without testset configuration to get mse and largestDifference for
  // final scenario configuration
  SGPP::datadriven::LearnerScenario scenario(fileName, lambda, gridConfig, SLESolverConfigRefine,
                                             SLESolverConfigFinal, adaptConfig);

  SGPP::datadriven::MetaLearner metaLearner(
      scenario.getGridConfig(), scenario.getSolverConfigurationRefine(),
      scenario.getSolverConfigurationFinal(), scenario.getAdaptivityConfiguration(),
      scenario.getLambda(), true);

  //  SGPP::datadriven::MetaLearner metaLearner(gridConfig, SLESolverConfigRefine,
  //  SLESolverConfigFinal,
  //                                            adaptConfig, lambda, true);

  SGPP::base::OCLOperationConfiguration parameters(parameterFileName);

  SGPP::datadriven::OperationMultipleEvalConfiguration configuration(
      SGPP::datadriven::OperationMultipleEvalType::STREAMING,
      SGPP::datadriven::OperationMultipleEvalSubType::OCLMP, parameters);

  std::string datasetFileNameFetched = scenario.getDatasetFileName();
  metaLearner.learn(configuration, datasetFileNameFetched, true);

  //  SGPP::base::Grid& grid = metaLearner.getLearnedGrid();
  SGPP::base::DataVector& alpha = metaLearner.getLearnedAlpha();

  alpha.toFile(alphaReferenceFileName);

  testsetConfig.hasTestDataset = true;
  testsetConfig.alphaReferenceFileName = alphaReferenceFileName;
  testsetConfig.expectedMSE = 0;
  testsetConfig.expectedLargestDifference = 0;

  SGPP::datadriven::LearnerScenario scenarioWithTest(fileName, lambda, gridConfig,
                                                     SLESolverConfigRefine, SLESolverConfigFinal,
                                                     adaptConfig, testsetConfig);

  scenarioWithTest.serialize(scenarioFileName);

  std::cout << "all done!" << std::endl;
}
