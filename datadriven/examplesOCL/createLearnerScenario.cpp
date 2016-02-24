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

  double lambda = 0.0000001;
  std::string fileName = "DR5_train.arff";
  //  std::string fileName = "friedman2_4d_500000.arff";
  //  std::string fileName = "friedman_4d.arff";
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
  SLESolverConfigFinal.maxIterations_ = 20;
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

  SGPP::base::OCLOperationConfiguration parameters("allDevices.cfg");

  SGPP::datadriven::OperationMultipleEvalConfiguration configuration(
      SGPP::datadriven::OperationMultipleEvalType::STREAMING,
      SGPP::datadriven::OperationMultipleEvalSubType::OCLMP, parameters);

  std::string datasetFileNameFetched = scenario.getDatasetFileName();
  metaLearner.learn(configuration, datasetFileNameFetched, true);

  SGPP::base::Grid& grid = metaLearner.getLearnedGrid();
  SGPP::base::DataVector& alpha = metaLearner.getLearnedAlpha();

  alpha.toString();
  //  SGPP::datadriven::Dataset dataset =
  // SGPP::datadriven::ARFFTools::readARFF(datasetFileNameFetched);
  //  std::unique_ptr<SGPP::base::OperationMultipleEval> eval =
  //      SGPP::op_factory::createOperationMultipleEval(grid, dataset.getData(), configuration);
  //
  //  SGPP::base::DataVector referenceValues(dataset.getNumberInstances());
  //
  //  eval->mult(alpha, referenceValues);
  //
  //  double mse = 0.0;
  //  double largestDifference = 0.0;
  //  for (size_t i = 0; i < referenceValues.getSize(); i++) {
  //    double difference = fabs(referenceValues[i] - dataset.getTargets()[i]);
  //
  //    if (difference > largestDifference) {
  //      largestDifference = difference;
  //    }
  //
  //    mse += difference * difference;
  //  }
  //
  //  mse /= static_cast<double>(referenceValues.getSize());

  //  testsetConfig.hasTestDataset = true;
  //  testsetConfig.datasetFileName = "DR5_train.arff";
  //  //  testsetConfig.datasetFileName = "friedman2_4d_500000.arff";
  //  //  testsetConfig.datasetFileName = "friedman_4d.arff";
  //  testsetConfig.expectedMSE = mse;
  //  testsetConfig.expectedLargestDifference = largestDifference;
  //
  //  SGPP::datadriven::LearnerScenario scenarioWithTest(fileName, lambda, gridConfig,
  //                                                     SLESolverConfigRefine,
  //                                                     SLESolverConfigFinal,
  //                                                     adaptConfig, testsetConfig);
  //
  //  scenarioWithTest.serialize("learnerSimple.scenario");
  //  //  scenario.writeToFile("learnerSimple.scenario");
  //
  //  SGPP::datadriven::LearnerScenario readScenario("learnerSimple.scenario");
  //
  //  //  readScenario.writeToFile("test.scenario");
  //  readScenario.serialize("test.scenario");

  std::cout << "all done!" << std::endl;
}
