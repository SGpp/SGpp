// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#if USE_OCL == 1

#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/predicate.hpp>

#include <limits>
#include <algorithm>
#include <vector>
#include <string>
#include <sstream>

#include "sgpp/globaldef.hpp"
#include "sgpp/base/opencl/OCLOperationConfiguration.hpp"
#include "sgpp/datadriven/application/StaticParameterTuner.hpp"
#include "sgpp/datadriven/application/MetaLearner.hpp"
#include "sgpp/base/exception/application_exception.hpp"
#include "sgpp/datadriven/DatadrivenOpFactory.hpp"
#include "sgpp/datadriven/tools/Dataset.hpp"
#include "sgpp/datadriven/tools/ARFFTools.hpp"

namespace SGPP {
namespace datadriven {

StaticParameterTuner::StaticParameterTuner(SGPP::base::OCLOperationConfiguration &fixedParameters,
                                           bool collectStatistics, bool verbose)
    : verbose(verbose),
      collectStatistics(collectStatistics),
      statisticsFolder("."),
      fixedParameters(fixedParameters) {}

void StaticParameterTuner::addParameter(const std::string &name,
                                        const std::vector<std::string> &valueRange) {
  this->tunableParameters.push_back(TunableParameter(name, valueRange));
}

SGPP::base::OCLOperationConfiguration StaticParameterTuner::tuneEverything(
    SGPP::datadriven::LearnerScenario &scenario, const std::string &kernelName) {
  if (scenario.getInternalPrecision() == InternalPrecision::Float) {
    this->fixedParameters.replaceIDAttr("INTERNAL_PRECISION", "float");
  } else {
    this->fixedParameters.replaceIDAttr("INTERNAL_PRECISION", "double");
  }

  std::vector<std::string> platformsCopy = this->fixedParameters["PLATFORMS"].keys();
  for (const std::string &platformName : platformsCopy) {
    json::Node &platformNode = this->fixedParameters["PLATFORMS"][platformName];
    std::vector<std::string> devicesCopy = platformNode["DEVICES"].keys();

    // temporarily remove all other platforms
    std::vector<std::string> otherPlatformsName;
    std::vector<std::unique_ptr<json::Node>> otherPlatforms;
    for (const std::string &otherPlatformName : platformsCopy) {
      if (otherPlatformName.compare(platformName) == 0) {
        continue;
      }
      otherPlatformsName.push_back(otherPlatformName);
      auto otherPlatformNode = this->fixedParameters["PLATFORMS"][otherPlatformName].erase();
      otherPlatforms.push_back(std::move(otherPlatformNode));
    }

    for (const std::string &deviceName : devicesCopy) {
      json::Node &deviceNode = platformNode["DEVICES"][deviceName];

      if (deviceNode.contains("COUNT")) {
        if (deviceNode["COUNT"].getUInt() == 0) {
          continue;
        }
      }

      std::cout << "tuning for device: " << deviceName << std::endl;

      // temporarily remove all other devices
      std::vector<std::string> otherDevicesName;
      std::vector<std::unique_ptr<json::Node>> otherDevices;
      for (const std::string &otherDeviceName : devicesCopy) {
        if (otherDeviceName.compare(deviceName) == 0) {
          continue;
        }
        otherDevicesName.push_back(otherDeviceName);
        auto otherDeviceNode = platformNode["DEVICES"][otherDeviceName].erase();
        otherDevices.push_back(std::move(otherDeviceNode));
      }

      // limit the number of devices used to 1 for tuning
      bool addedDeviceLimit = false;
      size_t oldLimitValue = 0;
      if (!deviceNode.contains("COUNT")) {
        deviceNode.addIDAttr("COUNT", 1ul);
        addedDeviceLimit = true;
      } else {
        oldLimitValue = deviceNode["COUNT"].getUInt();
      }

      // add an attribute for the kernel if none exists
      if (!deviceNode["KERNELS"].contains(kernelName)) {
        deviceNode["KERNELS"].addDictAttr(kernelName);
      }

      if (!deviceNode["KERNELS"][kernelName].contains("REUSE_SOURCE")) {
        deviceNode["KERNELS"][kernelName].addIDAttr("REUSE_SOURCE", false);
      }

      if (!deviceNode["KERNELS"][kernelName].contains("WRITE_SOURCE")) {
        deviceNode["KERNELS"][kernelName].addIDAttr("WRITE_SOURCE", false);
      }

      //            if (useDoublePrecision) {
      //                deviceNode["KERNELS"][kernelName].replaceTextAttr("INTERNAL_PRECISION",
      //                "double");
      //            } else {
      //                deviceNode["KERNELS"][kernelName].replaceTextAttr("INTERNAL_PRECISION",
      //                "float");
      //            }

      this->tuneParameters(scenario, platformName, deviceName, kernelName);

      if (collectStatistics) {
        std::string safePlatformName = platformName;
        std::replace(safePlatformName.begin(), safePlatformName.end(), ' ', '_');
        std::string safeDeviceName = deviceName;
        std::replace(safeDeviceName.begin(), safeDeviceName.end(), ' ', '_');
        std::string statisticsFileName =
            "statistics_" + safePlatformName + "_" + safeDeviceName + "_" + kernelName + ".csv";
        this->writeStatisticsToFile(statisticsFileName, platformName, deviceName, kernelName);
      }

      if (addedDeviceLimit) {
        deviceNode["COUNT"].erase();
      } else {
        deviceNode["COUNT"].setUInt(oldLimitValue);
      }

      // add the removed devices again for the next iteration
      for (size_t i = 0; i < otherDevices.size(); i++) {
        platformNode["DEVICES"].addAttribute(otherDevicesName[i], std::move(otherDevices[i]));
      }
    }

    // add the removed platforms again for the next iteration
    for (size_t i = 0; i < otherPlatforms.size(); i++) {
      this->fixedParameters["PLATFORMS"].addAttribute(otherPlatformsName[i],
                                                      std::move(otherPlatforms[i]));
    }
  }

  std::cout << "after tuning: " << std::endl;
  std::cout << this->fixedParameters << std::endl;

  return this->fixedParameters;
}

void StaticParameterTuner::setStatisticsFolder(const std::string &statisticsFolder) {
  this->statisticsFolder = statisticsFolder;
}

void StaticParameterTuner::tuneParameters(SGPP::datadriven::LearnerScenario &scenario,
                                          const std::string &platformName,
                                          const std::string &deviceName,
                                          const std::string &kernelName) {
  if (collectStatistics) {
    this->statistics.clear();
  }

  //    SGPP::base::OCLOperationConfiguration &currentParameters =
  //    fixedParameters;

  for (std::string &platformKey : fixedParameters["PLATFORMS"].keys()) {
    if (platformName.compare(platformKey) != 0) {
      throw;
    }
    for (std::string &deviceKey : fixedParameters["PLATFORMS"][platformKey]["DEVICES"].keys()) {
      if (deviceName.compare(deviceKey) != 0) {
        throw;
      }
      if (!fixedParameters["PLATFORMS"][platformName]["DEVICES"][deviceName]["KERNELS"].contains(
              kernelName)) {
        throw;
      }
    }
  }

  json::Node &kernelNode =
      fixedParameters["PLATFORMS"][platformName]["DEVICES"][deviceName]["KERNELS"][kernelName];

  // create initial parameter combination
  std::vector<size_t> valueIndices(tunableParameters.size());
  for (size_t i = 0; i < valueIndices.size(); i++) {
    valueIndices[i] = 0;
    TunableParameter &parameter = tunableParameters[i];
    kernelNode.replaceIDAttr(parameter.getName(), parameter.getValues()[0]);
  }

  std::cout << "-----------------------------------" << std::endl;
  for (std::string &key : kernelNode.keys()) {
    std::cout << "key: " << key << " value: " << kernelNode[key].get() << std::endl;
  }

  // evaluate initial parameter combination
  double shortestDuration;
  double highestGFlops;
  evaluateSetup(scenario, fixedParameters, kernelName, shortestDuration, highestGFlops);
  if (collectStatistics) {
    this->statistics.emplace_back(fixedParameters, shortestDuration, highestGFlops);
  }

  std::unique_ptr<json::Node> bestParameters(kernelNode.clone());

  size_t parameterIndex = 0;
  while (parameterIndex < tunableParameters.size()) {
    // enumerate next parameter combination
    // can the current index be increased?
    TunableParameter &parameter = tunableParameters[parameterIndex];
    if (valueIndices[parameterIndex] + 1 < parameter.getValues().size()) {
      valueIndices[parameterIndex] += 1;
      kernelNode[parameter.getName()].set(parameter.getValues()[valueIndices[parameterIndex]]);
      // reset lower indices
      for (size_t i = 0; i < parameterIndex; i++) {
        valueIndices[i] = 0;
        TunableParameter &parameterForReset = tunableParameters[i];
        kernelNode[parameterForReset.getName()].set(parameterForReset.getValues()[0]);
      }
      parameterIndex = 0;

      std::cout << "-----------------------------------" << std::endl;
      for (std::string &key : kernelNode.keys()) {
        std::cout << "key: " << key << " value: " << kernelNode[key].get() << std::endl;
      }

      // evaluate current parameter combination

      double duration;
      double GFlops;
      evaluateSetup(scenario, fixedParameters, kernelName, duration, GFlops);

      if (duration < shortestDuration) {
        std::cout << "new best combination! old: " << shortestDuration << " new: " << duration
                  << std::endl;
        shortestDuration = duration;
        highestGFlops = GFlops;
        bestParameters = std::unique_ptr<json::Node>(kernelNode.clone());
        std::cout << *bestParameters << std::endl;
      }
      if (collectStatistics) {
        this->statistics.emplace_back(fixedParameters, duration, GFlops);
      }
    } else {
      parameterIndex += 1;
    }
  }

  std::cout << "overall best parameters:" << std::endl;
  std::cout << *bestParameters << std::endl;

  kernelNode = *bestParameters;

  std::cout << "written parameters:" << std::endl;
  std::cout << this->fixedParameters["PLATFORMS"][platformName]["DEVICES"][deviceName]["KERNELS"]
                                    [kernelName] << std::endl;
}

double StaticParameterTuner::evaluateSetup(SGPP::datadriven::LearnerScenario &scenario,
                                           SGPP::base::OCLOperationConfiguration &currentParameters,
                                           const std::string &kernelName, double &duration,
                                           double &GFlops) {
  SGPP::datadriven::MetaLearner learner(
      scenario.getGridConfig(), scenario.getSolverConfigurationRefine(),
      scenario.getSolverConfigurationFinal(), scenario.getAdaptivityConfiguration(),
      scenario.getLambda(), verbose);

  SGPP::datadriven::OperationMultipleEvalType operationType;
  SGPP::datadriven::OperationMultipleEvalSubType operationSubType;

  if (kernelName.compare("StreamingOCLMultiPlatform") == 0) {
    operationType = SGPP::datadriven::OperationMultipleEvalType::STREAMING;
    operationSubType = SGPP::datadriven::OperationMultipleEvalSubType::OCLMP;
  } else if (kernelName.compare("StreamingModOCLFastMultiPlatform") == 0) {
    operationType = SGPP::datadriven::OperationMultipleEvalType::STREAMING;
    operationSubType = SGPP::datadriven::OperationMultipleEvalSubType::OCLFASTMP;
  } else if (kernelName.compare("StreamingModOCLMaskMultiPlatform") == 0) {
    operationType = SGPP::datadriven::OperationMultipleEvalType::STREAMING;
    operationSubType = SGPP::datadriven::OperationMultipleEvalSubType::OCLMASKMP;
  } else {
    throw SGPP::base::application_exception(
        "error: configured kernel is not known to static parameter tuner");
  }

  SGPP::datadriven::OperationMultipleEvalConfiguration configuration(
      operationType, operationSubType, currentParameters);

  std::string fileName = scenario.getDatasetFileName();

  duration = std::numeric_limits<double>::max();
  try {
    std::cout << "evaluating parameter combination" << std::endl;

    learner.learn(configuration, fileName);

    LearnerTiming timing = learner.getLearnerTiming();

    duration = timing.timeComplete_;

    GFlops = timing.GFlop_ / timing.timeComplete_;

    TestsetConfiguration testsetConfiguration = scenario.getTestsetConfiguration();

    if (testsetConfiguration.hasTestDataset) {
      this->verifyLearned(testsetConfiguration, learner.getLearnedAlpha());
    }
  } catch (SGPP::base::operation_exception &exception) {
    if (verbose) {
      std::cout << "invalid combination detected" << std::endl;
    }
  }

  if (verbose) {
    std::cout << "duration: " << duration << std::endl;
  }
  return duration;
}

void StaticParameterTuner::writeStatisticsToFile(const std::string &statisticsFileName,
                                                 const std::string &platformName,
                                                 const std::string &deviceName,
                                                 const std::string &kernelName) {
  if (!collectStatistics) {
    throw;
  }

  std::ofstream file(statisticsFolder + "/" + statisticsFileName);

  bool first = true;
  for (TunableParameter &columnParameter : this->tunableParameters) {
    if (!first) {
      file << ", ";
    } else {
      first = false;
    }
    file << columnParameter.getName();
  }
  if (this->tunableParameters.size() > 0) {
    file << ", ";
  }
  file << "duration, GFlops" << std::endl;

  for (auto &parameterDurationGFlopsTuple : this->statistics) {
    SGPP::base::OCLOperationConfiguration &parameter = std::get<0>(parameterDurationGFlopsTuple);
    json::Node &kernelNode =
        parameter["PLATFORMS"][platformName]["DEVICES"][deviceName]["KERNELS"][kernelName];
    double duration = std::get<1>(parameterDurationGFlopsTuple);
    double GFlops = std::get<2>(parameterDurationGFlopsTuple);

    first = true;
    for (TunableParameter &columnParameter : this->tunableParameters) {
      if (!first) {
        file << ", ";
      } else {
        first = false;
      }
      file << kernelNode[columnParameter.getName()].get();
    }
    if (this->tunableParameters.size() > 0) {
      file << ", ";
    }
    file << duration << ", " << GFlops << std::endl;
  }
}

void StaticParameterTuner::verifyLearned(TestsetConfiguration &testsetConfiguration,
                                         base::DataVector &alpha) {
  base::DataVector alphaReference =
      base::DataVector::fromFile(testsetConfiguration.alphaReferenceFileName);

  if (alphaReference.getSize() != alpha.getSize()) {
    throw base::application_exception("error: size of reference vector doesn't match");
  }

  double mse = 0.0;
  double largestDifference = 0.0;
  for (size_t i = 0; i < alpha.getSize(); i++) {
    double difference = fabs(alphaReference[i] - alpha[i]);

    if (difference > largestDifference) {
      largestDifference = difference;
    }

    mse += difference * difference;
  }
  mse /= static_cast<double>(alpha.getSize());

  if (mse > testsetConfiguration.expectedMSE ||
      largestDifference > testsetConfiguration.expectedLargestDifference) {
    std::string message("error: violated the expected error, mse: " + std::to_string(mse) +
                        " (excepted: " + std::to_string(testsetConfiguration.expectedMSE) +
                        ") largestDifference: " + std::to_string(largestDifference) +
                        " (excepted: " +
                        std::to_string(testsetConfiguration.expectedLargestDifference) + ")");
    throw base::application_exception(message.c_str());
  }
}

}  // namespace datadriven
}  // namespace SGPP

#endif
