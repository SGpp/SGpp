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
  std::string scenarioFileName("friedman2_4d_300000_StreamingOCLMultiPlatform_float.scenario");
  std::string parameterFile("reproduce.cfg");

  sgpp::datadriven::LearnerScenario scenario(scenarioFileName);

  bool verbose = true;
  sgpp::datadriven::MetaLearner learner(
      scenario.getGridConfig(), scenario.getSolverConfigurationRefine(),
      scenario.getSolverConfigurationFinal(), scenario.getAdaptivityConfiguration(),
      scenario.getLambda(), verbose);

  sgpp::base::OCLOperationConfiguration parameters(parameterFile);

  sgpp::datadriven::OperationMultipleEvalConfiguration configuration(
      sgpp::datadriven::OperationMultipleEvalType::STREAMING,
      sgpp::datadriven::OperationMultipleEvalSubType::OCLMP, parameters);

  std::string datasetFile = scenario.getDatasetFileName();
  learner.learn(configuration, datasetFile);

  return EXIT_SUCCESS;
}
