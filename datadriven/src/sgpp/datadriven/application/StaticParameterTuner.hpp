// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#if USE_OCL == 1

#include <utility>
#include <string>
#include <vector>

#include "sgpp/globaldef.hpp"
#include "sgpp/base/opencl/OCLOperationConfiguration.hpp"
#include "sgpp/datadriven/application/TunableParameter.hpp"
#include "sgpp/datadriven/application/LearnerScenario.hpp"
#include "sgpp/datadriven/operation/hash/simple/DatadrivenOperationCommon.hpp"

namespace sgpp {
namespace datadriven {

class StaticParameterTuner {
 private:
  bool collectStatistics;
  bool verbose;
  std::vector<std::pair<sgpp::base::OCLOperationConfiguration, double>> statistics;

  sgpp::base::OCLOperationConfiguration fixedParameters;

  std::vector<TunableParameter> tunableParameters;

  double evaluateSetup(sgpp::datadriven::LearnerScenario &scenario,
                       sgpp::base::OCLOperationConfiguration &currentParameters,
                       const std::string &kernelName);

  void writeStatisticsToFile(const std::string &statisticsFileName, const std::string &platformName,
                             const std::string &deviceName, const std::string &kernelName);

  void verifyLearned(TestsetConfiguration &testsetConfiguration, base::DataVector &alpha);

 public:
  StaticParameterTuner(sgpp::base::OCLOperationConfiguration &fixedParameters,
                       bool collectStatistics = false, bool verbose = false);

  void addParameter(const std::string &name, const std::vector<std::string> &valueRange);

  sgpp::base::OCLOperationConfiguration tuneEverything(sgpp::datadriven::LearnerScenario &scenario,
                                                       const std::string &kernelName);

  void tuneParameters(sgpp::datadriven::LearnerScenario &scenario, const std::string &platformName,
                      const std::string &deviceName, const std::string &kernelName);
};
}  // namespace datadriven
}  // namespace sgpp

#endif
