// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once
#if USE_OCL == 1

#include <utility>
#include <string>
#include <vector>
#include <tuple>

#include <sgpp/globaldef.hpp>
#include <sgpp/base/opencl/OCLOperationConfiguration.hpp>
#include <sgpp/datadriven/application/TunableParameter.hpp>
#include <sgpp/datadriven/application/LearnerScenario.hpp>
#include <sgpp/datadriven/operation/hash/DatadrivenOperationCommon.hpp>

namespace sgpp {
namespace datadriven {

class StaticParameterTuner {
 private:
  bool verbose;

  bool collectStatistics;
  std::string statisticsFolder;
  std::string scenarioName;

  std::vector<std::tuple<sgpp::base::OCLOperationConfiguration, double, double, double, double>>
      statistics;

  sgpp::base::OCLOperationConfiguration fixedParameters;

  std::vector<TunableParameter> tunableParameters;
  uint64_t configuredExperiments;
  uint64_t currentExperiment;

  double evaluateSetup(sgpp::datadriven::LearnerScenario &scenario,
                       sgpp::base::OCLOperationConfiguration &currentParameters,
                       const std::string &kernelName, double &duration, double &durationOperation,
                       double &durationKernel, double &GFlops);

  void writeStatisticsToFile(const std::string &statisticsFileName, const std::string &platformName,
                             const std::string &deviceName, const std::string &kernelName);

  void verifyLearned(TestsetConfiguration &testsetConfiguration, base::DataVector &alpha);

 public:
  StaticParameterTuner(sgpp::base::OCLOperationConfiguration &fixedParameters,
                       bool verbose = false);

  void enableStatistics(const std::string &statisticsFolder, const std::string &scenarioName);

  void addParameter(const std::string &name, const std::vector<std::string> &valueRange);

  sgpp::base::OCLOperationConfiguration tuneEverything(sgpp::datadriven::LearnerScenario &scenario,
                                                       const std::string &kernelName);

  void tuneParameters(sgpp::datadriven::LearnerScenario &scenario, const std::string &platformName,
                      const std::string &deviceName, const std::string &kernelName);
};
}  // namespace datadriven
}  // namespace sgpp

#endif
