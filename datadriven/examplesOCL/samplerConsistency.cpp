// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "sgpp/datadriven/application/MetaLearner.hpp"
#include "sgpp/datadriven/operation/hash/DatadrivenOperationCommon.hpp"
#include "sgpp/base/opencl/OCLOperationConfiguration.hpp"
#include "sgpp/datadriven/application/LearnerScenario.hpp"
#include "sgpp/base/datatypes/DataMatrix.hpp"
#include "sgpp/base/datatypes/DataVector.hpp"
#include "sgpp/base/exception/application_exception.hpp"

using sgpp::base::DataVector;
using sgpp::base::DataMatrix;

std::string baseFolder = "datadriven/performanceTests/scenarios/";

void verifyLearned(sgpp::datadriven::TestsetConfiguration &testsetConfiguration,
                   DataVector &alpha) {
  DataVector alphaReference = DataVector::fromFile(testsetConfiguration.alphaReferenceFileName);

  if (alphaReference.getSize() != alpha.getSize()) {
    throw sgpp::base::application_exception("error: size of reference vector doesn't match");
  }

  //  std::ofstream outFile("differences.log");

  double mse = 0.0;
  double largestDifference = 0.0;
  for (size_t i = 0; i < alpha.getSize(); i++) {
    double difference = fabs(alphaReference[i] - alpha[i]);

    //    outFile << "alpha[" << i << "] = " << alpha[i] << ", alphaReference[" << i
    //            << "] = " << alphaReference[i] << " difference: " << difference << std::endl;
    //    outFile << alpha[i] << std::endl;

    if (difference > largestDifference) {
      largestDifference = difference;
    }

    mse += difference * difference;
  }
  mse /= static_cast<double>(alpha.getSize());

  //  outFile.close();

  if (mse > testsetConfiguration.expectedMSE ||
      largestDifference > testsetConfiguration.expectedLargestDifference) {
    std::string message("error: violated the expected error, mse: " + std::to_string(mse) +
                        " (excepted: " + std::to_string(testsetConfiguration.expectedMSE) +
                        ") largestDifference: " + std::to_string(largestDifference) +
                        " (excepted: " +
                        std::to_string(testsetConfiguration.expectedLargestDifference) + ")");
    throw sgpp::base::application_exception(message.c_str());
  } else {
    std::cout << "mse: " << mse << " ok!" << std::endl;
    std::cout << "largestDifference: " << largestDifference << " ok!" << std::endl;
  }
}

int main(int argc, char **argv) {
  std::vector<std::string> scenarios = {
      baseFolder + "friedman2_4d_300000_Linear_double.scenario",
      baseFolder + "friedman2_4d_300000_Linear_float.scenario",
      baseFolder + "friedman2_4d_300000_ModLinear_double.scenario",
      baseFolder + "friedman2_4d_300000_ModLinear_float.scenario",
      baseFolder + "friedman1_10d_150000_Linear_double.scenario",
      baseFolder + "friedman1_10d_150000_Linear_float.scenario",
      baseFolder + "friedman1_10d_150000_ModLinear_double.scenario",
      baseFolder + "friedman1_10d_150000_ModLinear_float.scenario",
      baseFolder + "DR5_train_Linear_double.scenario",
      baseFolder + "DR5_train_Linear_float.scenario",
      baseFolder + "DR5_train_ModLinear_double.scenario",
      baseFolder + "DR5_train_ModLinear_float.scenario"};

  std::vector<std::string> parameters = {
      "friedman2_4d_300000_Linear_double_StreamingOCLMultiPlatform_tuned.cfg",
      "friedman2_4d_300000_Linear_float_StreamingOCLMultiPlatform_tuned.cfg",
      "friedman2_4d_300000_ModLinear_double_StreamingModOCLMaskMultiPlatform_tuned.cfg",
      "friedman2_4d_300000_ModLinear_float_StreamingModOCLMaskMultiPlatform_tuned.cfg",
      "friedman1_10d_150000_Linear_double_StreamingOCLMultiPlatform_tuned.cfg",
      "friedman1_10d_150000_Linear_float_StreamingOCLMultiPlatform_tuned.cfg",
      "friedman1_10d_150000_ModLinear_double_StreamingModOCLMaskMultiPlatform_tuned.cfg",
      "friedman1_10d_150000_ModLinear_float_StreamingModOCLMaskMultiPlatform_tuned.cfg",
      "DR5_train_Linear_double_StreamingOCLMultiPlatform_tuned.cfg",
      "DR5_train_Linear_float_StreamingOCLMultiPlatform_tuned.cfg",
      "DR5_train_ModLinear_double_StreamingModOCLMaskMultiPlatform_tuned.cfg",
      "DR5_train_ModLinear_float_StreamingModOCLMaskMultiPlatform_tuned.cfg"};

  std::ofstream outFile("consistency.log");
  for (size_t i = 0; i < scenarios.size(); i++) {
    std::string &scenarioFileName = scenarios[i];
    std::string &parameterFile = parameters[i];

    outFile << scenarioFileName << std::endl;

    sgpp::datadriven::LearnerScenario scenario(scenarioFileName);

    sgpp::datadriven::OperationMultipleEvalSubType subType;
    if (scenario["grid"]["type"].get().compare("Linear") == 0) {
      subType = sgpp::datadriven::OperationMultipleEvalSubType::OCLMP;
    } else if (scenario["grid"]["type"].get().compare("ModLinear") == 0) {
      subType = sgpp::datadriven::OperationMultipleEvalSubType::OCLMASKMP;
    } else {
      throw;
    }

    sgpp::base::OCLOperationConfiguration parameters(parameterFile);
    sgpp::datadriven::OperationMultipleEvalConfiguration configuration(
        sgpp::datadriven::OperationMultipleEvalType::STREAMING, subType, parameters);

    for (size_t repeat = 0; repeat < 10; repeat++) {
      std::cout << "repeat: " << repeat << std::endl;
      bool verbose = true;
      sgpp::datadriven::MetaLearner learner(
          scenario.getGridConfig(), scenario.getSolverConfigurationRefine(),
          scenario.getSolverConfigurationFinal(), scenario.getAdaptivityConfiguration(),
          scenario.getLambda(), verbose);

      std::string datasetFile = scenario.getDatasetFileName();
      try {
        learner.learn(configuration, datasetFile);
        //  learner.learnAndCompare(configuration, datasetFile, 4);

        DataVector &alpha = learner.getLearnedAlpha();
        sgpp::datadriven::TestsetConfiguration testsetConfiguration =
            scenario.getTestsetConfiguration();
        verifyLearned(testsetConfiguration, alpha);
        sgpp::datadriven::LearnerTiming timing = learner.getLearnerTiming();
        std::cout << "time complete: " << timing.timeComplete_ << std::endl;
        outFile << timing.timeComplete_ << std::endl;
      } catch (sgpp::base::operation_exception &e) {
        std::cout << "exception caught: " << e.what() << std::endl;
      }
    }
  }
  outFile.close();
  return EXIT_SUCCESS;
}
