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

namespace SGPP {
namespace datadriven {

class StaticParameterTuner {
 private:
  bool collectStatistics;
  bool verbose;
  std::vector<std::pair<SGPP::base::OCLOperationConfiguration, double>> statistics;

  SGPP::base::OCLOperationConfiguration fixedParameters;

  std::vector<TunableParameter> tunableParameters;

  double evaluateSetup(SGPP::datadriven::LearnerScenario &scenario,
                       SGPP::base::OCLOperationConfiguration &currentParameters,
                       const std::string &kernelName);

  void writeStatisticsToFile(const std::string &statisticsFileName, const std::string &platformName,
                             const std::string &deviceName, const std::string &kernelName);

  void verifyLearned(TestsetConfiguration &testsetConfiguration, base::DataVector &alpha);

 public:
  StaticParameterTuner(SGPP::base::OCLOperationConfiguration &fixedParameters,
                       bool collectStatistics = false, bool verbose = false);

  void addParameter(const std::string &name, const std::vector<std::string> &valueRange);

  SGPP::base::OCLOperationConfiguration tuneEverything(SGPP::datadriven::LearnerScenario &scenario,
                                                       const std::string &kernelName);

  void tuneParameters(SGPP::datadriven::LearnerScenario &scenario, const std::string &platformName,
                      const std::string &deviceName, const std::string &kernelName);
};
}  // namespace datadriven
}  // namespace SGPP

#endif
