// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/opencl/OCLOperationConfiguration.hpp>
#include <sgpp/base/operation/hash/OperationMultipleEval.hpp>
#include <sgpp/datadriven/DatadrivenOpFactory.hpp>
#include <sgpp/datadriven/application/LearnerScenario.hpp>
#include <sgpp/datadriven/application/MetaLearner.hpp>
#include <sgpp/datadriven/operation/hash/DatadrivenOperationCommon.hpp>
#include <sgpp/datadriven/tools/ARFFTools.hpp>
#include <sgpp/globaldef.hpp>

#include <iostream>
#include <string>

using sgpp::base::OCLOperationConfiguration;
using sgpp::datadriven::OperationMultipleEvalSubType;
using sgpp::datadriven::OperationMultipleEvalType;

int main(int argc, char** argv) {
  sgpp::base::RegularGridConfiguration gridConfig;
  sgpp::solver::SLESolverConfiguration SLESolverConfigRefine;
  sgpp::solver::SLESolverConfiguration SLESolverConfigFinal;
  sgpp::base::AdaptivityConfiguration adaptivityConfig;
  sgpp::datadriven::TestsetConfiguration testsetConfig;

  ///////////////////////////////// Configuration start ////////////////////////////////////

  // to be set manually in the generated scenario file:
  // testsetConfig.expectedMSE
  // testsetConfig.expectedLargestDifference

  sgpp::datadriven::InternalPrecision internalPrecision =
      sgpp::datadriven::InternalPrecision::Float;
  std::string datasetFileName = "DR5_train.arff";

  double lambda = 0.0000001;

  // setup grid
  gridConfig.dim_ = 0;  // dim is inferred from the data
  gridConfig.level_ = 10;
  gridConfig.type_ = sgpp::base::GridType::ModLinear;
  // dummy values
  gridConfig.boundaryLevel_ = 0;
  gridConfig.maxDegree_ = 30;

  // Set Adaptivity
  adaptivityConfig.maxLevelType_ = false;
  adaptivityConfig.numRefinementPoints_ = 80;
  adaptivityConfig.numRefinements_ = 0;
  adaptivityConfig.percent_ = 200.0;
  adaptivityConfig.refinementThreshold_ = 0.0;

  // Set solver during refinement
  SLESolverConfigRefine.eps_ = 0;
  SLESolverConfigRefine.maxIterations_ = 200;
  SLESolverConfigRefine.threshold_ = -1.0;
  SLESolverConfigRefine.type_ = sgpp::solver::SLESolverType::CG;

  // Set solver for final step
  SLESolverConfigFinal.eps_ = 0;
  SLESolverConfigFinal.maxIterations_ = 2;
  SLESolverConfigFinal.threshold_ = -1.0;
  SLESolverConfigFinal.type_ = sgpp::solver::SLESolverType::CG;

  ///////////////////////////////// Configuration end ////////////////////////////////////

  std::string parameterFileName;
  if (internalPrecision == sgpp::datadriven::InternalPrecision::Double) {
    parameterFileName = "createPlatformDouble.cfg";
  } else {
    parameterFileName = "createPlatformFloat.cfg";
  }

  std::string basisFunctionName;
  if (gridConfig.type_ == sgpp::base::GridType::Linear) {
    basisFunctionName = "Linear";
  } else if (gridConfig.type_ == sgpp::base::GridType::ModLinear) {
    basisFunctionName = "ModLinear";
  } else {
    throw;
  }

  // use the implementation of base for double!
  // use the default streaming implementation for float
  OperationMultipleEvalType operationType;
  OperationMultipleEvalSubType operationSubType;

  if (internalPrecision == sgpp::datadriven::InternalPrecision::Double) {
    operationType = OperationMultipleEvalType::DEFAULT;
    operationSubType = OperationMultipleEvalSubType::DEFAULT;

  } else {
    if (gridConfig.type_ == sgpp::base::GridType::Linear) {
      operationType = OperationMultipleEvalType::STREAMING;
      operationSubType = OperationMultipleEvalSubType::OCLMP;
    } else if (gridConfig.type_ == sgpp::base::GridType::ModLinear) {
      operationType = OperationMultipleEvalType::STREAMING;
      operationSubType = OperationMultipleEvalSubType::OCLMASKMP;
    } else {
      throw;
    }
  }

  size_t dotPosition = datasetFileName.find('.');
  std::string datasetName = datasetFileName.substr(0, dotPosition);

  std::string precisionString;
  if (internalPrecision == sgpp::datadriven::InternalPrecision::Float) {
    precisionString = "float";
  } else {
    precisionString = "double";
  }
  std::string scenarioFileName(datasetName + "_" + basisFunctionName + "_" + precisionString +
                               ".scenario");

  std::string alphaReferenceFileName =
      datasetName + "_" + basisFunctionName + "_" + precisionString + "_AlphaReference.vec";

  // create scenario to learn without testset configuration to get mse and largestDifference for
  // final scenario configuration
  sgpp::datadriven::LearnerScenario scenario(datasetFileName, lambda, internalPrecision, gridConfig,
                                             SLESolverConfigRefine, SLESolverConfigFinal,
                                             adaptivityConfig);

  sgpp::datadriven::MetaLearner metaLearner(
      scenario.getGridConfig(), scenario.getSolverConfigurationRefine(),
      scenario.getSolverConfigurationFinal(), scenario.getAdaptivityConfiguration(),
      scenario.getLambda(), true);

  //  sgpp::datadriven::MetaLearner metaLearner(gridConfig, SLESolverConfigRefine,
  //  SLESolverConfigFinal,
  //                                            adaptivityConfig, lambda, true);

  OCLOperationConfiguration parameters(parameterFileName);

  // set precision in configuration as specified
  parameters.replaceIDAttr("INTERNAL_PRECISION", precisionString);

  sgpp::datadriven::OperationMultipleEvalConfiguration configuration(operationType,
                                                                     operationSubType, parameters);

  std::string datasetFileNameFetched = scenario.getDatasetFileName();
  //  metaLearner.learn(configuration, datasetFileNameFetched, true);
  metaLearner.learn(configuration, datasetFileNameFetched, true);

  sgpp::base::DataVector& alpha = metaLearner.getLearnedAlpha();

  alpha.toFile(alphaReferenceFileName);

  testsetConfig.hasTestDataset = true;
  testsetConfig.alphaReferenceFileName = alphaReferenceFileName;
  testsetConfig.expectedMSE = 0;
  testsetConfig.expectedLargestDifference = 0;

  sgpp::datadriven::LearnerScenario scenarioWithTest(
      datasetFileName, lambda, internalPrecision, gridConfig, SLESolverConfigRefine,
      SLESolverConfigFinal, adaptivityConfig, testsetConfig);

  scenarioWithTest.serialize(scenarioFileName);

  std::cout << "all done!" << std::endl;
}
