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

using sgpp::datadriven::OperationMultipleEvalType;
using sgpp::datadriven::OperationMultipleEvalSubType;
using sgpp::base::OCLOperationConfiguration;

int main(int argc, char** argv) {
  ///////////////////////////////// Configuration start ////////////////////////////////////

  // to be set manually in the generated scenario file:
  // testsetConfig.expectedMSE
  // testsetConfig.expectedLargestDifference

  sgpp::datadriven::InternalPrecision internalPrecision =
      sgpp::datadriven::InternalPrecision::Float;
  //  std::string kernelName = "StreamingOCLMultiPlatform";
  std::string datasetFileName = "friedman2_4d_300000.arff";

  OperationMultipleEvalType operationType = OperationMultipleEvalType::STREAMING;
  OperationMultipleEvalSubType operationSubType = OperationMultipleEvalSubType::OCLMASKMP;

  std::string kernelName;
  if (operationType == OperationMultipleEvalType::STREAMING) {
    if (operationSubType == OperationMultipleEvalSubType::OCLMASKMP) {
      kernelName = "StreamingModOCLMaskMultiPlatform";
    } else if (operationSubType == OperationMultipleEvalSubType::OCLFASTMP) {
      kernelName = "StreamingModOCLFastMultiPlatform";
    } else {
      throw;
    }
  } else {
    throw;
  }

  std::string parameterFileName("allDevices.cfg");
  double lambda = 0.0000001;

  sgpp::base::RegularGridConfiguration gridConfig;
  sgpp::solver::SLESolverConfiguration SLESolverConfigRefine;
  sgpp::solver::SLESolverConfiguration SLESolverConfigFinal;
  sgpp::base::AdpativityConfiguration adaptConfig;
  sgpp::datadriven::TestsetConfiguration testsetConfig;

  // setup grid
  gridConfig.dim_ = 0;  // dim is inferred from the data
  gridConfig.level_ = 10;
  gridConfig.type_ = sgpp::base::GridType::ModLinear;
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
  SLESolverConfigRefine.type_ = sgpp::solver::SLESolverType::CG;

  // Set solver for final step
  SLESolverConfigFinal.eps_ = 0;
  SLESolverConfigFinal.maxIterations_ = 1;
  SLESolverConfigFinal.threshold_ = -1.0;
  SLESolverConfigFinal.type_ = sgpp::solver::SLESolverType::CG;

  ///////////////////////////////// Configuration end ////////////////////////////////////

  size_t dotPosition = datasetFileName.find('.');
  std::string datasetName = datasetFileName.substr(0, dotPosition);

  std::string precisionString;
  if (internalPrecision == sgpp::datadriven::InternalPrecision::Float) {
    precisionString = "float";
  } else {
    precisionString = "double";
  }
  std::string scenarioFileName(datasetName + "_" + kernelName + "_" + precisionString +
                               ".scenario");

  std::string alphaReferenceFileName =
      datasetName + "_" + kernelName + "_" + precisionString + "_AlphaReference.vec";

  // create scenario to learn without testset configuration to get mse and largestDifference for
  // final scenario configuration
  sgpp::datadriven::LearnerScenario scenario(datasetFileName, lambda, internalPrecision, gridConfig,
                                             SLESolverConfigRefine, SLESolverConfigFinal,
                                             adaptConfig);

  sgpp::datadriven::MetaLearner metaLearner(
      scenario.getGridConfig(), scenario.getSolverConfigurationRefine(),
      scenario.getSolverConfigurationFinal(), scenario.getAdaptivityConfiguration(),
      scenario.getLambda(), true);

  //  sgpp::datadriven::MetaLearner metaLearner(gridConfig, SLESolverConfigRefine,
  //  SLESolverConfigFinal,
  //                                            adaptConfig, lambda, true);

  OCLOperationConfiguration parameters(parameterFileName);

  // set precision in configuration as specified
  parameters.replaceIDAttr("INTERNAL_PRECISION", precisionString);

  sgpp::datadriven::OperationMultipleEvalConfiguration configuration(operationType,
                                                                     operationSubType, parameters);

  std::string datasetFileNameFetched = scenario.getDatasetFileName();
  metaLearner.learn(configuration, datasetFileNameFetched, true);

  sgpp::base::DataVector& alpha = metaLearner.getLearnedAlpha();

  alpha.toFile(alphaReferenceFileName);

  testsetConfig.hasTestDataset = true;
  testsetConfig.alphaReferenceFileName = alphaReferenceFileName;
  testsetConfig.expectedMSE = 0;
  testsetConfig.expectedLargestDifference = 0;

  sgpp::datadriven::LearnerScenario scenarioWithTest(
      datasetFileName, lambda, internalPrecision, gridConfig, SLESolverConfigRefine,
      SLESolverConfigFinal, adaptConfig, testsetConfig);

  scenarioWithTest.serialize(scenarioFileName);

  std::cout << "all done!" << std::endl;
}
