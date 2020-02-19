// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#if USE_OCL == 1

#include <sgpp/base/exception/application_exception.hpp>
#include <sgpp/base/opencl/OCLOperationConfiguration.hpp>
#include <sgpp/datadriven/DatadrivenOpFactory.hpp>
#include <sgpp/datadriven/application/MetaLearner.hpp>
#include <sgpp/datadriven/application/StaticParameterTuner.hpp>
#include <sgpp/datadriven/tools/ARFFTools.hpp>
#include <sgpp/datadriven/tools/Dataset.hpp>
#include <sgpp/globaldef.hpp>

#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/predicate.hpp>

#include <algorithm>
#include <limits>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

namespace sgpp {
namespace datadriven {

StaticParameterTuner::StaticParameterTuner(sgpp::base::OCLOperationConfiguration &fixedParameters,
                                           bool verbose)
    : verbose(verbose),
      collectStatistics(false),
      statisticsFolder("."),
      scenarioName(""),
      fixedParameters(fixedParameters),
      configuredExperiments(0),
      currentExperiment(0) {}

void StaticParameterTuner::addParameter(const std::string &name,
                                        const std::vector<std::string> &valueRange) {
  this->tunableParameters.push_back(TunableParameter(name, valueRange));
}

sgpp::base::OCLOperationConfiguration StaticParameterTuner::tuneEverything(
    sgpp::datadriven::LearnerScenario &scenario, const std::string &kernelName) {
  uint64_t configuredDevices = 0;
  for (const std::string &platformName : this->fixedParameters["PLATFORMS"].keys()) {
    json::Node &platformNode = this->fixedParameters["PLATFORMS"][platformName];
    for (const std::string &deviceName : platformNode["DEVICES"].keys()) {
      json::Node &deviceNode = platformNode["DEVICES"][deviceName];
      if (!deviceNode.contains("COUNT") || deviceNode["COUNT"].getUInt() > 0) {
        configuredDevices += 1;
      }
    }
  }
  uint64_t currentDevice = 1;

  // use the precision configuration from the scenario and overwrite the precision from the
  // submitted configuration
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

      if (verbose) {
        std::cout << "tuning for device: " << deviceName << " (" << currentDevice << "/"
                  << configuredDevices << ")" << std::endl;
      }
      currentDevice += 1;

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
      bool addedCountLimit = false;
      bool hasOldCountLimit = false;
      size_t oldCountLimit = 0;
      if (!deviceNode.contains("SELECT")) {
        if (!deviceNode.contains("COUNT")) {
          deviceNode.addIDAttr("COUNT", UINT64_C(1));
          addedCountLimit = true;
        } else {
          oldCountLimit = deviceNode["COUNT"].getUInt();
          hasOldCountLimit = true;
        }
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

      // if there is no explicit specification, set schedule size to a very large value
      // this makes tuning easier, as kernels are running longer
      //      bool addedScheduleSize = false;
      if (!deviceNode["KERNELS"][kernelName].contains("KERNEL_SCHEDULE_SIZE")) {
        //        addedScheduleSize = true;
        // TODO(pfandedd): improve, for now multiples of 1024 should run with any kernel
        deviceNode["KERNELS"][kernelName].addIDAttr("KERNEL_SCHEDULE_SIZE", UINT64_C(1024000));
      }

      // special case for intel:
      // if "OPTIMIZATION_FLAGS" contains "-cl-strict-aliasing", the kernel cannot be build on
      // intel,
      // even though this is a standard flag
      std::unique_ptr<TunableParameter> optimizationFlagsCopy;
      if (platformName.compare("Intel(R) OpenCL") == 0) {
        for (TunableParameter &parameter : tunableParameters) {
          if (parameter.getName().compare("OPTIMIZATION_FLAGS") == 0) {
            // save the old values
            optimizationFlagsCopy = std::make_unique<TunableParameter>(parameter);
            // change the new values
            for (std::string &parameterValue : parameter.getValues()) {
              size_t found = parameterValue.find("-cl-strict-aliasing");
              if (found != std::string::npos) {
                parameterValue.replace(found, 19, "");
              }
            }
          }
        }
      }

      this->tuneParameters(scenario, platformName, deviceName, kernelName);

      // replace the "OPTIMIZATION_FLAGS" values with the original ones for the intel platform
      if (optimizationFlagsCopy.operator bool()) {
        for (TunableParameter &parameter : tunableParameters) {
          if (parameter.getName().compare("OPTIMIZATION_FLAGS") == 0) {
            // reset to the old values
            for (size_t i = 0; i < parameter.getValues().size(); i++) {
              parameter.getValues()[i] = optimizationFlagsCopy->getValues()[i];
            }
          }
        }
      }

      if (collectStatistics) {
        std::string safePlatformName = platformName;
        std::replace(safePlatformName.begin(), safePlatformName.end(), ' ', '_');
        std::string safeDeviceName = deviceName;
        std::replace(safeDeviceName.begin(), safeDeviceName.end(), ' ', '_');

        std::string statisticsFileName =
            "statistics_" + safePlatformName + "_" + safeDeviceName + "_" + scenarioName + ".csv";
        this->writeStatisticsToFile(statisticsFileName, platformName, deviceName, kernelName);
      }

      // keep the schedule size for now, makes other experiments easier
      //      if (addedScheduleSize) {
      //        deviceNode["KERNELS"][kernelName]["KERNEL_SCHEDULE_SIZE"].erase();
      //      }

      if (addedCountLimit) {
        deviceNode["COUNT"].erase();
      } else {
        if (hasOldCountLimit) {
          deviceNode["COUNT"].setUInt(oldCountLimit);
        }
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

  if (verbose) {
    std::cout << "after tuning: " << std::endl;
    std::cout << this->fixedParameters << std::endl;
  }

  return this->fixedParameters;
}

void StaticParameterTuner::enableStatistics(const std::string &statisticsFolder,
                                            const std::string &scenarioName) {
  this->collectStatistics = true;
  this->statisticsFolder = statisticsFolder;
  this->scenarioName = scenarioName;
}

void StaticParameterTuner::tuneParameters(sgpp::datadriven::LearnerScenario &scenario,
                                          const std::string &platformName,
                                          const std::string &deviceName,
                                          const std::string &kernelName) {
  if (collectStatistics) {
    this->statistics.clear();
  }

  for (size_t i = 0; i < tunableParameters.size(); i++) {
    TunableParameter &tunableParameter = tunableParameters[i];
    std::cout << tunableParameter.getName() << " values: " << tunableParameter.getValues().size()
              << std::endl;
    if (i == 0) {
      configuredExperiments = tunableParameter.getValues().size();
    }
    configuredExperiments *= tunableParameter.getValues().size();
  }
  currentExperiment = 1;

  //    sgpp::base::OCLOperationConfiguration &currentParameters =
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
    kernelNode.replaceTextAttr(parameter.getName(), parameter.getValues()[0]);
  }

  if (verbose) {
    std::cout << "-----------------------------------" << std::endl;
    for (std::string &key : kernelNode.keys()) {
      std::cout << "key: " << key << " value: " << kernelNode[key].get() << std::endl;
    }
    std::cout << "-----------------------------------" << std::endl;
    std::cout << "evaluating combination (" << currentExperiment << "/" << configuredExperiments
              << ") for current device" << std::endl;
  }
  currentExperiment += 1;

  // evaluate initial parameter combination
  double shortestDuration;
  double shortestDurationOperations;
  double shortestDurationKernels;
  double highestGFlops;
  evaluateSetup(scenario, fixedParameters, kernelName, shortestDuration, shortestDurationOperations,
                shortestDurationKernels, highestGFlops);
  if (collectStatistics) {
    this->statistics.emplace_back(fixedParameters, shortestDuration, shortestDurationOperations,
                                  shortestDurationKernels, highestGFlops);
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

      if (verbose) {
        std::cout << "evaluating combination (" << currentExperiment << "/" << configuredExperiments
                  << ") for current device" << std::endl;
      }
      currentExperiment += 1;

      // evaluate current parameter combination
      double duration;
      double durationOperations;
      double durationKernels;
      double GFlops;
      evaluateSetup(scenario, fixedParameters, kernelName, duration, durationOperations,
                    durationKernels, GFlops);

      if (duration < shortestDuration) {
        std::cout << "new best combination! old: " << shortestDuration << " new: " << duration
                  << std::endl;
        shortestDuration = duration;
        highestGFlops = GFlops;
        bestParameters = std::unique_ptr<json::Node>(kernelNode.clone());
        //        std::cout << *bestParameters << std::endl;
      }
      if (collectStatistics) {
        this->statistics.emplace_back(fixedParameters, duration, durationOperations,
                                      durationKernels, GFlops);
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

double StaticParameterTuner::evaluateSetup(sgpp::datadriven::LearnerScenario &scenario,
                                           sgpp::base::OCLOperationConfiguration &currentParameters,
                                           const std::string &kernelName, double &duration,
                                           double &durationOperation, double &durationKernel,
                                           double &GFlops) {
  sgpp::datadriven::MetaLearner learner(
      scenario.getGridConfig(), scenario.getSolverConfigurationRefine(),
      scenario.getSolverConfigurationFinal(), scenario.getAdaptivityConfiguration(),
      scenario.getLambda(), verbose);

  sgpp::datadriven::OperationMultipleEvalType operationType;
  sgpp::datadriven::OperationMultipleEvalSubType operationSubType;

  if (kernelName.compare("StreamingOCLMultiPlatform") == 0) {
    operationType = sgpp::datadriven::OperationMultipleEvalType::STREAMING;
    operationSubType = sgpp::datadriven::OperationMultipleEvalSubType::OCLMP;
  } else if (kernelName.compare("StreamingModOCLFastMultiPlatform") == 0) {
    operationType = sgpp::datadriven::OperationMultipleEvalType::STREAMING;
    operationSubType = sgpp::datadriven::OperationMultipleEvalSubType::OCLFASTMP;
  } else if (kernelName.compare("StreamingModOCLMaskMultiPlatform") == 0) {
    operationType = sgpp::datadriven::OperationMultipleEvalType::STREAMING;
    operationSubType = sgpp::datadriven::OperationMultipleEvalSubType::OCLMASKMP;
  } else if (kernelName.compare("StreamingModOCLOpt") == 0) {
    operationType = sgpp::datadriven::OperationMultipleEvalType::STREAMING;
    operationSubType = sgpp::datadriven::OperationMultipleEvalSubType::OCLOPT;

  } else {
    throw sgpp::base::application_exception(
        "error: configured kernel is not known to static parameter tuner");
  }

  sgpp::datadriven::OperationMultipleEvalConfiguration configuration(
      operationType, operationSubType, currentParameters);

  std::string fileName = scenario.getDatasetFileName();

  duration = std::numeric_limits<double>::max();
  try {
    std::cout << "evaluating parameter combination" << std::endl;

    learner.learn(configuration, fileName);

    LearnerTiming timing = learner.getLearnerTiming();

    duration = timing.timeComplete_;
    durationOperation = timing.timeMultComplete_ + timing.timeMultTransComplete_;
    durationKernel = timing.timeMultCompute_ + timing.timeMultTransCompute_;

    GFlops = timing.GFlop_ / timing.timeComplete_;

    TestsetConfiguration testsetConfiguration = scenario.getTestsetConfiguration();

    if (testsetConfiguration.hasTestDataset) {
      this->verifyLearned(testsetConfiguration, learner.getLearnedAlpha());
    }
  } catch (sgpp::base::operation_exception &exception) {
    if (verbose) {
      std::cout << "invalid combination detected:" << exception.what() << std::endl;
    }
  }

  if (verbose) {
    std::cout << "duration: " << duration << std::endl;
  }
  return duration;
}

bool StaticParameterTunerStatisticsTupleComparer(
    std::tuple<sgpp::base::OCLOperationConfiguration, double, double, double, double> left,
    std::tuple<sgpp::base::OCLOperationConfiguration, double, double, double, double> right) {
  return std::get<4>(left) < std::get<4>(right);
}

void StaticParameterTuner::writeStatisticsToFile(const std::string &statisticsFileName,
                                                 const std::string &platformName,
                                                 const std::string &deviceName,
                                                 const std::string &kernelName) {
  if (!collectStatistics) {
    throw;
  }

  std::sort(this->statistics.begin(), this->statistics.end(),
            StaticParameterTunerStatisticsTupleComparer);

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
  file << "duration, durationOperations, durationKernels, GFlops" << std::endl;

  for (auto &parameterDurationGFlopsTuple : this->statistics) {
    sgpp::base::OCLOperationConfiguration &parameter = std::get<0>(parameterDurationGFlopsTuple);
    json::Node &kernelNode =
        parameter["PLATFORMS"][platformName]["DEVICES"][deviceName]["KERNELS"][kernelName];
    double duration = std::get<1>(parameterDurationGFlopsTuple);
    double durationOperations = std::get<2>(parameterDurationGFlopsTuple);
    double durationKernels = std::get<3>(parameterDurationGFlopsTuple);
    double GFlops = std::get<4>(parameterDurationGFlopsTuple);

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
    file << duration << ", " << durationOperations << ", " << durationKernels << ", " << GFlops
         << std::endl;
  }
}

void StaticParameterTuner::verifyLearned(TestsetConfiguration &testsetConfiguration,
                                         base::DataVector &alpha) {
  base::DataVector alphaReference;
  try {
    alphaReference = base::DataVector::fromFile(testsetConfiguration.alphaReferenceFileName);
  } catch (std::exception &e) {
    throw base::application_exception("error: coult not open alpha verification file");
  }

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
    std::stringstream errorStream;
    errorStream << std::scientific;
    errorStream << "error: violated the expected error, mse: " << mse;
    errorStream << " (excepted: " << testsetConfiguration.expectedMSE;
    errorStream << ") largestDifference: ";
    errorStream << largestDifference;
    errorStream << " (excepted: " << testsetConfiguration.expectedLargestDifference << ")";
    //    std::string message("error: violated the expected error, mse: " + std::to_string(mse) +
    //                        " (excepted: " + std::to_string(testsetConfiguration.expectedMSE) +
    //                        ") largestDifference: " + std::to_string(largestDifference) +
    //                        " (excepted: " +
    //                        std::to_string(testsetConfiguration.expectedLargestDifference) + ")");
    std::string message(errorStream.str());
    throw base::application_exception(message.c_str());
  } else {
    if (verbose) {
      std::cout << "verification passed (mse: " << mse
                << ", largestDifference: " << largestDifference << ")" << std::endl;
    }
  }
}

}  // namespace datadriven
}  // namespace sgpp

#endif
