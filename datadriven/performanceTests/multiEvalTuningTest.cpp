// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#if USE_OCL == 1

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <boost/filesystem.hpp>

#include <zlib.h>

#include <random>
#include <fstream>
#include <iostream>
#include <chrono>
#include <map>
#include <string>
#include <vector>

#include "testsCommon.hpp"

#include "sgpp/datadriven/application/LearnerScenario.hpp"
#include "sgpp/datadriven/application/StaticParameterTuner.hpp"

BOOST_AUTO_TEST_SUITE(AutoTuningPaper)

BOOST_AUTO_TEST_CASE(Friedman2_4d_Linear_Float) {
  // internal precision is specified by the scenario, the parameter configuration is overwritten
  std::string scenarioFileName = "friedman2_4d_300000_StreamingOCLMultiPlatform_float.scenario";
  std::string parameterConfigurationFileName = "platformFloat.cfg";
  std::string kernelName = "StreamingOCLMultiPlatform";
  bool collectStatistics = true;

  size_t dotPosition = scenarioFileName.find('.');
  std::string scenarioFileNamePrefix = scenarioFileName.substr(0, dotPosition);
  std::string outputFileName = scenarioFileNamePrefix + "_tuned.cfg";

  sgpp::datadriven::LearnerScenario scenario(scenarioFileName);
  sgpp::base::OCLOperationConfiguration parameter(parameterConfigurationFileName);
  sgpp::datadriven::StaticParameterTuner staticParameterTuner(parameter, true);

  std::string statisticsFolderName = "statistics";

  if (collectStatistics) {
    staticParameterTuner.enableStatistics(statisticsFolderName, scenarioFileNamePrefix);
    try {
      if (boost::filesystem::create_directory(statisticsFolderName)) {
        BOOST_MESSAGE("created output directory: " << statisticsFolderName);
      }
    } catch (boost::filesystem::filesystem_error &e) {
      BOOST_FAIL("could not create statistics output folder: " << statisticsFolderName << ": "
                                                               << e.what());
    }
  }

  staticParameterTuner.addParameter("KERNEL_USE_LOCAL_MEMORY", {"true", "false"});
  staticParameterTuner.addParameter("KERNEL_DATA_BLOCK_SIZE", {"1", "2", "4"});
  staticParameterTuner.addParameter("KERNEL_TRANS_GRID_BLOCK_SIZE", {"1", "2", "4"});
  staticParameterTuner.addParameter("KERNEL_STORE_DATA", {"register", "array"});
  staticParameterTuner.addParameter("KERNEL_MAX_DIM_UNROLL", {"10", "1"});
  staticParameterTuner.addParameter("LOCAL_SIZE", {"128", "256"});
  staticParameterTuner.addParameter("VERBOSE", {"true"});

  sgpp::base::OCLOperationConfiguration bestParameters =
      staticParameterTuner.tuneEverything(scenario, kernelName);

  bestParameters.serialize(outputFileName);
}

BOOST_AUTO_TEST_CASE(Friedman2_4d_Linear_Double) {
  // internal precision is specified by the scenario, the parameter configuration is overwritten
  std::string scenarioFileName = "friedman2_4d_300000_StreamingOCLMultiPlatform_double.scenario";
  std::string parameterConfigurationFileName = "platformDouble.cfg";
  std::string kernelName = "StreamingOCLMultiPlatform";
  bool collectStatistics = true;

  size_t dotPosition = scenarioFileName.find('.');
  std::string scenarioFileNamePrefix = scenarioFileName.substr(0, dotPosition);
  std::string outputFileName = scenarioFileNamePrefix + "_tuned.cfg";

  sgpp::datadriven::LearnerScenario scenario(scenarioFileName);
  sgpp::base::OCLOperationConfiguration parameter(parameterConfigurationFileName);
  sgpp::datadriven::StaticParameterTuner staticParameterTuner(parameter, true);

  std::string statisticsFolderName = "statistics";

  if (collectStatistics) {
    staticParameterTuner.enableStatistics(statisticsFolderName, scenarioFileNamePrefix);
    try {
      if (boost::filesystem::create_directory(statisticsFolderName)) {
        BOOST_MESSAGE("created output directory: " << statisticsFolderName);
      }
    } catch (boost::filesystem::filesystem_error &e) {
      BOOST_FAIL("could not create statistics output folder: " << statisticsFolderName << ": "
                                                               << e.what());
    }
  }

  staticParameterTuner.addParameter("KERNEL_USE_LOCAL_MEMORY", {"true", "false"});
  staticParameterTuner.addParameter("KERNEL_DATA_BLOCK_SIZE", {"1", "2", "4"});
  staticParameterTuner.addParameter("KERNEL_TRANS_GRID_BLOCK_SIZE", {"1", "2", "4"});
  staticParameterTuner.addParameter("KERNEL_STORE_DATA", {"register", "array"});
  staticParameterTuner.addParameter("KERNEL_MAX_DIM_UNROLL", {"10", "1"});
  staticParameterTuner.addParameter("LOCAL_SIZE", {"128", "256"});
  staticParameterTuner.addParameter("VERBOSE", {"true"});

  sgpp::base::OCLOperationConfiguration bestParameters =
      staticParameterTuner.tuneEverything(scenario, kernelName);

  bestParameters.serialize(outputFileName);
}

BOOST_AUTO_TEST_CASE(Friedman2_4d_ModLinearMask_Float) {
  // internal precision is specified by the scenario, the parameter configuration is overwritten
  std::string scenarioFileName =
      "friedman2_4d_300000_StreamingModOCLMaskMultiPlatform_float.scenario";
  std::string parameterConfigurationFileName = "platformFloat.cfg";
  std::string kernelName = "StreamingModOCLMaskMultiPlatform";
  bool collectStatistics = true;

  size_t dotPosition = scenarioFileName.find('.');
  std::string scenarioFileNamePrefix = scenarioFileName.substr(0, dotPosition);
  std::string outputFileName = scenarioFileNamePrefix + "_tuned.cfg";

  sgpp::datadriven::LearnerScenario scenario(scenarioFileName);
  sgpp::base::OCLOperationConfiguration parameter(parameterConfigurationFileName);
  sgpp::datadriven::StaticParameterTuner staticParameterTuner(parameter, true);

  std::string statisticsFolderName = "statistics";

  if (collectStatistics) {
    staticParameterTuner.enableStatistics(statisticsFolderName, scenarioFileNamePrefix);
    try {
      if (boost::filesystem::create_directory(statisticsFolderName)) {
        BOOST_MESSAGE("created output directory: " << statisticsFolderName);
      }
    } catch (boost::filesystem::filesystem_error &e) {
      BOOST_FAIL("could not create statistics output folder: " << statisticsFolderName << ": "
                                                               << e.what());
    }
  }

  staticParameterTuner.addParameter("KERNEL_USE_LOCAL_MEMORY", {"false", "true"});
  staticParameterTuner.addParameter("KERNEL_DATA_BLOCK_SIZE", {"1", "2", "4"});
  staticParameterTuner.addParameter("KERNEL_TRANS_GRID_BLOCK_SIZE", {"1", "2", "4"});
  staticParameterTuner.addParameter("KERNEL_STORE_DATA", {"register", "array"});
  staticParameterTuner.addParameter("KERNEL_MAX_DIM_UNROLL", {"10", "1"});
  staticParameterTuner.addParameter("LOCAL_SIZE", {"128", "256"});
  staticParameterTuner.addParameter("VERBOSE", {"true"});

  sgpp::base::OCLOperationConfiguration bestParameters =
      staticParameterTuner.tuneEverything(scenario, kernelName);

  bestParameters.serialize(outputFileName);
}

BOOST_AUTO_TEST_CASE(Friedman2_4d_ModLinearMask_Double) {
  // internal precision is specified by the scenario, the parameter configuration is overwritten
  std::string scenarioFileName =
      "friedman2_4d_300000_StreamingModOCLMaskMultiPlatform_double.scenario";
  std::string parameterConfigurationFileName = "platformDouble.cfg";
  std::string kernelName = "StreamingModOCLMaskMultiPlatform";
  bool collectStatistics = true;

  size_t dotPosition = scenarioFileName.find('.');
  std::string scenarioFileNamePrefix = scenarioFileName.substr(0, dotPosition);
  std::string outputFileName = scenarioFileNamePrefix + "_tuned.cfg";

  sgpp::datadriven::LearnerScenario scenario(scenarioFileName);
  sgpp::base::OCLOperationConfiguration parameter(parameterConfigurationFileName);
  sgpp::datadriven::StaticParameterTuner staticParameterTuner(parameter, true);

  std::string statisticsFolderName = "statistics";

  if (collectStatistics) {
    staticParameterTuner.enableStatistics(statisticsFolderName, scenarioFileNamePrefix);
    try {
      if (boost::filesystem::create_directory(statisticsFolderName)) {
        BOOST_MESSAGE("created output directory: " << statisticsFolderName);
      }
    } catch (boost::filesystem::filesystem_error &e) {
      BOOST_FAIL("could not create statistics output folder: " << statisticsFolderName << ": "
                                                               << e.what());
    }
  }

  staticParameterTuner.addParameter("KERNEL_USE_LOCAL_MEMORY", {"false", "true"});
  staticParameterTuner.addParameter("KERNEL_DATA_BLOCK_SIZE", {"1", "2", "4"});
  staticParameterTuner.addParameter("KERNEL_TRANS_GRID_BLOCK_SIZE", {"1", "2", "4"});
  staticParameterTuner.addParameter("KERNEL_STORE_DATA", {"register", "array"});
  staticParameterTuner.addParameter("KERNEL_MAX_DIM_UNROLL", {"10", "1"});
  staticParameterTuner.addParameter("LOCAL_SIZE", {"128", "256"});
  staticParameterTuner.addParameter("VERBOSE", {"true"});

  sgpp::base::OCLOperationConfiguration bestParameters =
      staticParameterTuner.tuneEverything(scenario, kernelName);

  bestParameters.serialize(outputFileName);
}

BOOST_AUTO_TEST_CASE(Friedman2_4d_ModLinearFast_Float) {
  // internal precision is specified by the scenario, the parameter configuration is overwritten
  std::string scenarioFileName =
      "friedman2_4d_300000_StreamingModOCLFastMultiPlatform_float.scenario";
  std::string parameterConfigurationFileName = "platformFloat.cfg";
  std::string kernelName = "StreamingModOCLFastMultiPlatform";
  bool collectStatistics = true;

  size_t dotPosition = scenarioFileName.find('.');
  std::string scenarioFileNamePrefix = scenarioFileName.substr(0, dotPosition);
  std::string outputFileName = scenarioFileNamePrefix + "_tuned.cfg";

  sgpp::datadriven::LearnerScenario scenario(scenarioFileName);
  sgpp::base::OCLOperationConfiguration parameter(parameterConfigurationFileName);
  sgpp::datadriven::StaticParameterTuner staticParameterTuner(parameter, true);

  std::string statisticsFolderName = "statistics";

  if (collectStatistics) {
    staticParameterTuner.enableStatistics(statisticsFolderName, scenarioFileNamePrefix);
    try {
      if (boost::filesystem::create_directory(statisticsFolderName)) {
        BOOST_MESSAGE("created output directory: " << statisticsFolderName);
      }
    } catch (boost::filesystem::filesystem_error &e) {
      BOOST_FAIL("could not create statistics output folder: " << statisticsFolderName << ": "
                                                               << e.what());
    }
  }

  staticParameterTuner.addParameter("KERNEL_USE_LOCAL_MEMORY", {"false", "true"});
  staticParameterTuner.addParameter("KERNEL_DATA_BLOCK_SIZE", {"1", "2", "4"});
  staticParameterTuner.addParameter("KERNEL_TRANS_DATA_BLOCK_SIZE", {"1", "2", "4"});
  staticParameterTuner.addParameter("KERNEL_TRANS_GRID_BLOCK_SIZE", {"1", "2", "4"});
  staticParameterTuner.addParameter("KERNEL_STORE_DATA", {"register", "array"});
  staticParameterTuner.addParameter("KERNEL_MAX_DIM_UNROLL", {"10", "1"});
  staticParameterTuner.addParameter("LOCAL_SIZE", {"128", "256"});
  staticParameterTuner.addParameter("VERBOSE", {"true"});

  sgpp::base::OCLOperationConfiguration bestParameters =
      staticParameterTuner.tuneEverything(scenario, kernelName);

  bestParameters.serialize(outputFileName);
}

BOOST_AUTO_TEST_CASE(Friedman2_4d_ModLinearFast_Double) {
  // internal precision is specified by the scenario, the parameter configuration is overwritten
  std::string scenarioFileName =
      "friedman2_4d_300000_StreamingModOCLFastMultiPlatform_double.scenario";
  std::string parameterConfigurationFileName = "platformDouble.cfg";
  std::string kernelName = "StreamingModOCLFastMultiPlatform";
  bool collectStatistics = true;

  size_t dotPosition = scenarioFileName.find('.');
  std::string scenarioFileNamePrefix = scenarioFileName.substr(0, dotPosition);
  std::string outputFileName = scenarioFileNamePrefix + "_tuned.cfg";

  sgpp::datadriven::LearnerScenario scenario(scenarioFileName);
  sgpp::base::OCLOperationConfiguration parameter(parameterConfigurationFileName);
  sgpp::datadriven::StaticParameterTuner staticParameterTuner(parameter, true);

  std::string statisticsFolderName = "statistics";

  if (collectStatistics) {
    staticParameterTuner.enableStatistics(statisticsFolderName, scenarioFileNamePrefix);
    try {
      if (boost::filesystem::create_directory(statisticsFolderName)) {
        BOOST_MESSAGE("created output directory: " << statisticsFolderName);
      }
    } catch (boost::filesystem::filesystem_error &e) {
      BOOST_FAIL("could not create statistics output folder: " << statisticsFolderName << ": "
                                                               << e.what());
    }
  }

  staticParameterTuner.addParameter("KERNEL_USE_LOCAL_MEMORY", {"false", "true"});
  staticParameterTuner.addParameter("KERNEL_DATA_BLOCK_SIZE", {"1", "2", "4"});
  staticParameterTuner.addParameter("KERNEL_TRANS_DATA_BLOCK_SIZE", {"1", "2", "4"});
  staticParameterTuner.addParameter("KERNEL_TRANS_GRID_BLOCK_SIZE", {"1", "2", "4"});
  staticParameterTuner.addParameter("KERNEL_STORE_DATA", {"register", "array"});
  staticParameterTuner.addParameter("KERNEL_MAX_DIM_UNROLL", {"10", "1"});
  staticParameterTuner.addParameter("LOCAL_SIZE", {"128", "256"});
  staticParameterTuner.addParameter("VERBOSE", {"true"});

  sgpp::base::OCLOperationConfiguration bestParameters =
      staticParameterTuner.tuneEverything(scenario, kernelName);

  bestParameters.serialize(outputFileName);
}

BOOST_AUTO_TEST_SUITE_END()

#endif
