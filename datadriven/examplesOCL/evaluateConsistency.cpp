// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/exception/application_exception.hpp>
#include <sgpp/base/opencl/OCLOperationConfiguration.hpp>
#include <sgpp/datadriven/application/LearnerScenario.hpp>
#include <sgpp/datadriven/application/MetaLearner.hpp>
#include <sgpp/datadriven/operation/hash/DatadrivenOperationCommon.hpp>

#include <cstring>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>


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
    std::stringstream errorStream;
    errorStream << "error: violated the expected error, mse: " << mse
                << " (excepted: " << testsetConfiguration.expectedMSE
                << ") largestDifference: " << largestDifference
                << " (excepted: " << testsetConfiguration.expectedLargestDifference << ")";
    std::string message = errorStream.str();
    throw sgpp::base::application_exception(message.c_str());
  } else {
    std::stringstream messageStream;
    messageStream << std::scientific;
    messageStream << "mse: " << mse << " ok! (excepted: " << testsetConfiguration.expectedMSE
                  << ")";
    messageStream << " largestDifference: " << largestDifference
                  << " ok! (expected: " << testsetConfiguration.expectedLargestDifference << ")"
                  << std::endl;
    std::cout << messageStream.str();
  }
}

int main(int argc, char **argv) {
  if (argc < 2 || argc > 3) {
    throw;
  }

  bool floatOnly = false;
  if (argc == 3) {
    if (strcmp(argv[2], "floatOnly") == 0) {
      std::cout << "float-only enabled!" << std::endl;
      floatOnly = true;
    } else {
      throw;
    }
  }

  std::vector<std::string> scenarios = {
      baseFolder + "chess_4d_500000_Linear_double.scenario",
      baseFolder + "chess_4d_500000_Linear_float.scenario",
      baseFolder + "chess_4d_500000_ModLinear_double.scenario",
      baseFolder + "chess_4d_500000_ModLinear_float.scenario",
      baseFolder + "friedman1_10d_500000_Linear_double.scenario",
      baseFolder + "friedman1_10d_500000_Linear_float.scenario",
      baseFolder + "friedman1_10d_500000_ModLinear_double.scenario",
      baseFolder + "friedman1_10d_500000_ModLinear_float.scenario",
      baseFolder + "DR5_train_Linear_double.scenario",
      baseFolder + "DR5_train_Linear_float.scenario",
      baseFolder + "DR5_train_ModLinear_double.scenario",
      baseFolder + "DR5_train_ModLinear_float.scenario"};

  std::vector<std::string> parameters = {
      "chess_4d_500000_Linear_double_StreamingOCLMultiPlatform_tuned.cfg",
      "chess_4d_500000_Linear_float_StreamingOCLMultiPlatform_tuned.cfg",
      "chess_4d_500000_ModLinear_double_StreamingModOCLMaskMultiPlatform_tuned.cfg",
      "chess_4d_500000_ModLinear_float_StreamingModOCLMaskMultiPlatform_tuned.cfg",
      "friedman1_10d_500000_Linear_double_StreamingOCLMultiPlatform_tuned.cfg",
      "friedman1_10d_500000_Linear_float_StreamingOCLMultiPlatform_tuned.cfg",
      "friedman1_10d_500000_ModLinear_double_StreamingModOCLMaskMultiPlatform_tuned.cfg",
      "friedman1_10d_500000_ModLinear_float_StreamingModOCLMaskMultiPlatform_tuned.cfg",
      "DR5_train_Linear_double_StreamingOCLMultiPlatform_tuned.cfg",
      "DR5_train_Linear_float_StreamingOCLMultiPlatform_tuned.cfg",
      "DR5_train_ModLinear_double_StreamingModOCLMaskMultiPlatform_tuned.cfg",
      "DR5_train_ModLinear_float_StreamingModOCLMaskMultiPlatform_tuned.cfg"};

  std::vector<bool> isDouble = {true, false, true, false, true, false,
                                true, false, true, false, true, false};

  std::ofstream outFile(std::string(argv[1]) + "_consistency.log");
  for (size_t i = 0; i < scenarios.size(); i++) {
    std::string &scenarioFileName = scenarios[i];
    std::string &parameterFile = parameters[i];

    if (floatOnly && isDouble[i]) {
      std::cout << "skipping..." << std::endl;
      continue;
    }

    std::cout << "scenario: " << scenarioFileName << std::endl;
    std::cout << argv[0] << std::endl;
    std::cout << argv[1] << std::endl;

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

    for (size_t repeat = 0; repeat < 20; repeat++) {
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
