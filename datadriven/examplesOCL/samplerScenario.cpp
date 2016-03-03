// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <iostream>
#include <string>

#include "sgpp/datadriven/application/MetaLearner.hpp"
#include "sgpp/datadriven/operation/hash/simple/DatadrivenOperationCommon.hpp"
#include "sgpp/base/opencl/OCLOperationConfiguration.hpp"
#include "sgpp/datadriven/application/LearnerScenario.hpp"

int main(int argc, char** argv) {
  std::string scenarioFileName(
      "friedman2_4d_300000_StreamingModOCLMaskMultiPlatform_double.scenario");
  std::string parameterFile("reproduce.cfg");

  SGPP::datadriven::LearnerScenario scenario(scenarioFileName);

  bool verbose = true;
  SGPP::datadriven::MetaLearner learner(
      scenario.getGridConfig(), scenario.getSolverConfigurationRefine(),
      scenario.getSolverConfigurationFinal(), scenario.getAdaptivityConfiguration(),
      scenario.getLambda(), verbose);

  SGPP::base::OCLOperationConfiguration parameters(parameterFile);

  SGPP::datadriven::OperationMultipleEvalConfiguration configuration(
      SGPP::datadriven::OperationMultipleEvalType::STREAMING,
      SGPP::datadriven::OperationMultipleEvalSubType::OCLMASKMP, parameters);

  std::string datasetFile = scenario.getDatasetFileName();
  learner.learn(configuration, datasetFile);

  return EXIT_SUCCESS;
}
