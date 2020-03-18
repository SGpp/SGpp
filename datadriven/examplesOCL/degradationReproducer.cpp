// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/datadriven/application/MetaLearner.hpp>
#include <sgpp/datadriven/operation/hash/DatadrivenOperationCommon.hpp>
#include <sgpp/base/opencl/OCLOperationConfiguration.hpp>
#include <sgpp/datadriven/application/LearnerScenario.hpp>
#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/exception/application_exception.hpp>

#include <fstream>
#include <iostream>
#include <string>


using sgpp::base::DataVector;
using sgpp::base::DataMatrix;

std::string baseFolder = "datadriven/performanceTests/scenarios/";

void verifyLearned(sgpp::datadriven::TestsetConfiguration &testsetConfiguration,
                   DataVector &alpha) {
  std::cout << "opening file" << std::endl;
  DataVector alphaReference = DataVector::fromFile(testsetConfiguration.alphaReferenceFileName);
  std::cout << "opening done" << std::endl;
  if (alphaReference.getSize() != alpha.getSize()) {
    std::cout << "error: size of reference vector doesn't match" << std::endl;
  }

  std::ofstream outFile("differences.log");

  double mse = 0.0;
  double largestDifference = 0.0;
  for (size_t i = 0; i < alpha.getSize(); i++) {
    double difference = fabs(alphaReference[i] - alpha[i]);

    outFile << "alpha[" << i << "] = " << alpha[i] << ", alphaReference[" << i
            << "] = " << alphaReference[i] << " difference: " << difference << std::endl;
    //    outFile << alpha[i] << std::endl;

    if (difference > largestDifference) {
      largestDifference = difference;
    }

    mse += difference * difference;
  }
  mse /= static_cast<double>(alpha.getSize());

  outFile.close();

  if (mse > testsetConfiguration.expectedMSE ||
      largestDifference > testsetConfiguration.expectedLargestDifference) {
    std::string message("error: violated the expected error, mse: " + std::to_string(mse) +
                        " (excepted: " + std::to_string(testsetConfiguration.expectedMSE) +
                        ") largestDifference: " + std::to_string(largestDifference) +
                        " (excepted: " +
                        std::to_string(testsetConfiguration.expectedLargestDifference) + ")");
    //    throw sgpp::base::application_exception(message.c_str());
    std::cout << message.c_str() << std::endl;
  } else {
    std::cout << "mse: " << mse << " ok!" << std::endl;
    std::cout << "largestDifference: " << largestDifference << " ok!" << std::endl;
  }
}

int main(int argc, char **argv) {
  //  std::string scenarioFileName =
  //      baseFolder + "friedman2_4d_300000_StreamingModOCLMaskMultiPlatform_double.scenario";
  //  std::string scenarioFileName =
  //      baseFolder + "friedman2_4d_300000_StreamingModOCLFastMultiPlatform_double.scenario";
  std::string scenarioFileName =
      baseFolder + "friedman2_4d_300000_StreamingOCLMultiPlatform_double.scenario";
  //  std::string scenarioFileName =
  //      baseFolder + "friedman2_4d_300000_StreamingModOCLOpt_double.scenario";

  //  std::string scenarioFileName =
  //      baseFolder + "friedman2_4d_300000_StreamingModOCLMaskMultiPlatform_float.scenario";
  // std::string scenarioFileName =
  //     baseFolder + "friedman2_4d_300000_StreamingModOCLFastMultiPlatform_float.scenario";
  //  std::string scenarioFileName =
  //      baseFolder + "friedman2_4d_300000_StreamingOCLMultiPlatform_float.scenario";
  //  std::string scenarioFileName =
  //      baseFolder + "friedman2_4d_300000_StreamingModOCLOpt_float.scenario";

  std::string parameterFile("hawaii_double_tuned.cfg");

  sgpp::base::OCLOperationConfiguration parameters(parameterFile);
  sgpp::datadriven::OperationMultipleEvalConfiguration configuration(
      sgpp::datadriven::OperationMultipleEvalType::STREAMING,
      sgpp::datadriven::OperationMultipleEvalSubType::OCLMP, parameters);

  sgpp::datadriven::LearnerScenario scenario(scenarioFileName);

  bool verbose = true;
  sgpp::datadriven::MetaLearner learner(
      scenario.getGridConfig(), scenario.getSolverConfigurationRefine(),
      scenario.getSolverConfigurationFinal(), scenario.getAdaptivityConfiguration(),
      scenario.getLambda(), verbose);

  std::string datasetFile = scenario.getDatasetFileName();
  for (size_t i = 0; i < 100; i++) {
    try {
      learner.learn(configuration, datasetFile);
      //  learner.learnAndCompare(configuration, datasetFile, 4);

      //      DataVector &alpha = learner.getLearnedAlpha();
      sgpp::datadriven::TestsetConfiguration testsetConfiguration =
          scenario.getTestsetConfiguration();
      //      std::cout << "before verify" << std::endl;
      //      verifyLearned(testsetConfiguration, alpha);
      //      std::cout << "after verify" << std::endl;
    } catch (sgpp::base::operation_exception &e) {
      std::cout << "exception caught: " << e.what() << std::endl;
    }
  }
  return EXIT_SUCCESS;
}
