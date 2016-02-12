// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#if USE_OCL == 1

#include <sgpp/globaldef.hpp>

#include <sgpp/base/opencl/OCLOperationConfiguration.hpp>
#include <sgpp/datadriven/application/TunableParameter.hpp>
#include <sgpp/datadriven/application/LearnerScenario.hpp>

namespace SGPP {
namespace datadriven {

class StaticParameterTuner {
 private:
  bool collectStatistics;
  bool verbose;
  std::vector<std::pair<SGPP::base::OCLOperationConfiguration, double>>
      statistics;

  SGPP::base::OCLOperationConfiguration fixedParameters;

  std::vector<TunableParameter> tunableParameters;

  //    json::Node &kernelNode;

  double evaluateSetup(SGPP::datadriven::LearnerScenario &scenario,
                       SGPP::base::OCLOperationConfiguration &currentParameters,
                       const std::string &kernelName);

  void writeStatisticsToFile(const std::string &statisticsFileName,
                             const std::string &platformName,
                             const std::string &deviceName,
                             const std::string &kernelName);

 public:
  StaticParameterTuner(SGPP::base::OCLOperationConfiguration &fixedParameters,
                       bool collectStatistics = false, bool verbose = false);

  ////    //write a file with detailed stats for the optimization
  //    StaticParameterTuner(const std::string &tunerFileName,
  //    SGPP::base::OCLOperationConfiguration &fixedParameters,
  //            bool collectStatistics = false, bool verbose = false);

  //    void addFixedParameter(const std::string &name, const std::string
  //    &value, const ParameterType type);

  void addParameter(const std::string &name,
                    const std::vector<std::string> &valueRange);

  SGPP::base::OCLOperationConfiguration tuneEverything(
      SGPP::datadriven::LearnerScenario &scenario,
      const std::string &kernelName);

  void tuneParameters(SGPP::datadriven::LearnerScenario &scenario,
                      const std::string &platformName,
                      const std::string &deviceName,
                      const std::string &kernelName);

  //    void writeToFile(const std::string &fileName);
  //
  //    void readFromFile(const std::string &fileName);
};
}
}

#endif
