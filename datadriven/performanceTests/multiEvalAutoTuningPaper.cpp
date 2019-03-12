// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifdef ZLIB
#if USE_OCL == 1

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <boost/filesystem.hpp>

#include <zlib.h>

#include <chrono>
#include <fstream>
#include <iostream>
#include <map>
#include <random>
#include <string>
#include <vector>

#include "testsCommon.hpp"

#include "sgpp/datadriven/application/LearnerScenario.hpp"
#include "sgpp/datadriven/application/MetaLearner.hpp"
#include "sgpp/datadriven/application/StaticParameterTuner.hpp"

std::string scenarioBaseDir = "datadriven/performanceTests/scenarios/";

BOOST_AUTO_TEST_SUITE(AutoTuningPaper)

BOOST_AUTO_TEST_CASE(Chess_4d_Linear_Float) {
  // internal precision is specified by the scenario, the parameter configuration is overwritten
  std::string scenarioFileName = "chess_4d_500000_Linear_float.scenario";
  std::string parameterConfigurationFileName = "platformFloat.cfg";
  std::string kernelName = "StreamingOCLMultiPlatform";
  bool collectStatistics = true;

  size_t dotPosition = scenarioFileName.find('.');
  std::string scenarioFileNamePrefix = scenarioFileName.substr(0, dotPosition);
  scenarioFileNamePrefix = scenarioFileNamePrefix.append("_" + kernelName);
  std::string outputFileName = scenarioFileNamePrefix + "_tuned.cfg";

  sgpp::datadriven::LearnerScenario scenario(scenarioBaseDir + scenarioFileName);
  sgpp::base::OCLOperationConfiguration parameter(parameterConfigurationFileName);
  sgpp::datadriven::StaticParameterTuner staticParameterTuner(parameter, true);

  std::string statisticsFolderName = "statistics";

  if (collectStatistics) {
    staticParameterTuner.enableStatistics(statisticsFolderName, scenarioFileNamePrefix);
    try {
      if (boost::filesystem::create_directory(statisticsFolderName)) {
        BOOST_TEST_MESSAGE("created output directory: " << statisticsFolderName);
      }
    } catch (boost::filesystem::filesystem_error &e) {
      BOOST_FAIL("could not create statistics output folder: " << statisticsFolderName << ": "
                                                               << e.what());
    }
  }

  staticParameterTuner.addParameter("KERNEL_USE_LOCAL_MEMORY", {"true", "false"});
  staticParameterTuner.addParameter("KERNEL_DATA_BLOCK_SIZE", {"1", "2", "4"});
  staticParameterTuner.addParameter("KERNEL_TRANS_GRID_BLOCK_SIZE", {"1", "2", "4"});
  staticParameterTuner.addParameter("KERNEL_STORE_DATA", {"array"});
  staticParameterTuner.addParameter("KERNEL_MAX_DIM_UNROLL", {"10", "1"});
  staticParameterTuner.addParameter("LOCAL_SIZE", {"128", "256"});
  staticParameterTuner.addParameter("VERBOSE", {"true"});
  staticParameterTuner.addParameter("OPTIMIZATION_FLAGS",
                                    {"-cl-strict-aliasing -cl-fast-relaxed-math"});

  sgpp::base::OCLOperationConfiguration bestParameters =
      staticParameterTuner.tuneEverything(scenario, kernelName);

  bestParameters.serialize(outputFileName);
}

BOOST_AUTO_TEST_CASE(Chess_4d_Linear_Double) {
  // internal precision is specified by the scenario, the parameter configuration is overwritten
  std::string scenarioFileName = "chess_4d_500000_Linear_double.scenario";
  std::string parameterConfigurationFileName = "platformDouble.cfg";
  std::string kernelName = "StreamingOCLMultiPlatform";
  bool collectStatistics = true;

  size_t dotPosition = scenarioFileName.find('.');
  std::string scenarioFileNamePrefix = scenarioFileName.substr(0, dotPosition);
  scenarioFileNamePrefix = scenarioFileNamePrefix.append("_" + kernelName);
  std::string outputFileName = scenarioFileNamePrefix + "_tuned.cfg";

  sgpp::datadriven::LearnerScenario scenario(scenarioBaseDir + scenarioFileName);
  sgpp::base::OCLOperationConfiguration parameter(parameterConfigurationFileName);
  sgpp::datadriven::StaticParameterTuner staticParameterTuner(parameter, true);

  std::string statisticsFolderName = "statistics";

  if (collectStatistics) {
    staticParameterTuner.enableStatistics(statisticsFolderName, scenarioFileNamePrefix);
    try {
      if (boost::filesystem::create_directory(statisticsFolderName)) {
        BOOST_TEST_MESSAGE("created output directory: " << statisticsFolderName);
      }
    } catch (boost::filesystem::filesystem_error &e) {
      BOOST_FAIL("could not create statistics output folder: " << statisticsFolderName << ": "
                                                               << e.what());
    }
  }

  staticParameterTuner.addParameter("KERNEL_USE_LOCAL_MEMORY", {"true", "false"});
  staticParameterTuner.addParameter("KERNEL_DATA_BLOCK_SIZE", {"1", "2", "4"});
  staticParameterTuner.addParameter("KERNEL_TRANS_GRID_BLOCK_SIZE", {"1", "2", "4"});
  staticParameterTuner.addParameter("KERNEL_STORE_DATA", {"array"});
  staticParameterTuner.addParameter("KERNEL_MAX_DIM_UNROLL", {"10", "1"});
  staticParameterTuner.addParameter("LOCAL_SIZE", {"128", "256"});
  staticParameterTuner.addParameter("VERBOSE", {"true"});
  staticParameterTuner.addParameter("OPTIMIZATION_FLAGS",
                                    {"-cl-strict-aliasing -cl-fast-relaxed-math"});

  sgpp::base::OCLOperationConfiguration bestParameters =
      staticParameterTuner.tuneEverything(scenario, kernelName);

  bestParameters.serialize(outputFileName);
}

BOOST_AUTO_TEST_CASE(Chess_4d_ModLinear_Float) {
  // internal precision is specified by the scenario, the parameter configuration is overwritten
  std::string scenarioFileName = "chess_4d_500000_ModLinear_float.scenario";
  std::string parameterConfigurationFileName = "platformFloat.cfg";
  std::string kernelName = "StreamingModOCLMaskMultiPlatform";
  bool collectStatistics = true;

  size_t dotPosition = scenarioFileName.find('.');
  std::string scenarioFileNamePrefix = scenarioFileName.substr(0, dotPosition);
  scenarioFileNamePrefix = scenarioFileNamePrefix.append("_" + kernelName);
  std::string outputFileName = scenarioFileNamePrefix + "_tuned.cfg";

  sgpp::datadriven::LearnerScenario scenario(scenarioBaseDir + scenarioFileName);
  sgpp::base::OCLOperationConfiguration parameter(parameterConfigurationFileName);
  sgpp::datadriven::StaticParameterTuner staticParameterTuner(parameter, true);

  std::string statisticsFolderName = "statistics";

  if (collectStatistics) {
    staticParameterTuner.enableStatistics(statisticsFolderName, scenarioFileNamePrefix);
    try {
      if (boost::filesystem::create_directory(statisticsFolderName)) {
        BOOST_TEST_MESSAGE("created output directory: " << statisticsFolderName);
      }
    } catch (boost::filesystem::filesystem_error &e) {
      BOOST_FAIL("could not create statistics output folder: " << statisticsFolderName << ": "
                                                               << e.what());
    }
  }

  staticParameterTuner.addParameter("KERNEL_USE_LOCAL_MEMORY", {"false", "true"});
  staticParameterTuner.addParameter("KERNEL_DATA_BLOCK_SIZE", {"1", "2", "4"});
  staticParameterTuner.addParameter("KERNEL_TRANS_GRID_BLOCK_SIZE", {"1", "2", "4"});
  staticParameterTuner.addParameter("KERNEL_STORE_DATA", {"array"});
  staticParameterTuner.addParameter("KERNEL_MAX_DIM_UNROLL", {"10", "1"});
  staticParameterTuner.addParameter("LOCAL_SIZE", {"128", "256"});
  staticParameterTuner.addParameter("VERBOSE", {"true"});
  staticParameterTuner.addParameter("OPTIMIZATION_FLAGS",
                                    {"-cl-strict-aliasing -cl-fast-relaxed-math"});

  sgpp::base::OCLOperationConfiguration bestParameters =
      staticParameterTuner.tuneEverything(scenario, kernelName);

  bestParameters.serialize(outputFileName);
}

BOOST_AUTO_TEST_CASE(Chess_4d_ModLinear_Double) {
  // internal precision is specified by the scenario, the parameter configuration is overwritten
  std::string scenarioFileName = "chess_4d_500000_ModLinear_double.scenario";
  std::string parameterConfigurationFileName = "platformDouble.cfg";
  std::string kernelName = "StreamingModOCLMaskMultiPlatform";
  bool collectStatistics = true;

  size_t dotPosition = scenarioFileName.find('.');
  std::string scenarioFileNamePrefix = scenarioFileName.substr(0, dotPosition);
  scenarioFileNamePrefix = scenarioFileNamePrefix.append("_" + kernelName);
  std::string outputFileName = scenarioFileNamePrefix + "_tuned.cfg";

  sgpp::datadriven::LearnerScenario scenario(scenarioBaseDir + scenarioFileName);
  sgpp::base::OCLOperationConfiguration parameter(parameterConfigurationFileName);
  sgpp::datadriven::StaticParameterTuner staticParameterTuner(parameter, true);

  std::string statisticsFolderName = "statistics";

  if (collectStatistics) {
    staticParameterTuner.enableStatistics(statisticsFolderName, scenarioFileNamePrefix);
    try {
      if (boost::filesystem::create_directory(statisticsFolderName)) {
        BOOST_TEST_MESSAGE("created output directory: " << statisticsFolderName);
      }
    } catch (boost::filesystem::filesystem_error &e) {
      BOOST_FAIL("could not create statistics output folder: " << statisticsFolderName << ": "
                                                               << e.what());
    }
  }

  staticParameterTuner.addParameter("KERNEL_USE_LOCAL_MEMORY", {"false", "true"});
  staticParameterTuner.addParameter("KERNEL_DATA_BLOCK_SIZE", {"1", "2", "4"});
  staticParameterTuner.addParameter("KERNEL_TRANS_GRID_BLOCK_SIZE", {"1", "2", "4"});
  staticParameterTuner.addParameter("KERNEL_STORE_DATA", {"array"});
  staticParameterTuner.addParameter("KERNEL_MAX_DIM_UNROLL", {"10", "1"});
  staticParameterTuner.addParameter("LOCAL_SIZE", {"128", "256"});
  staticParameterTuner.addParameter("VERBOSE", {"true"});
  staticParameterTuner.addParameter("OPTIMIZATION_FLAGS",
                                    {"-cl-strict-aliasing -cl-fast-relaxed-math"});

  sgpp::base::OCLOperationConfiguration bestParameters =
      staticParameterTuner.tuneEverything(scenario, kernelName);

  bestParameters.serialize(outputFileName);
}

BOOST_AUTO_TEST_CASE(Friedman2_4d_Linear_Float) {
  // internal precision is specified by the scenario, the parameter configuration is overwritten
  std::string scenarioFileName = "friedman2_4d_500000_Linear_float.scenario";
  std::string parameterConfigurationFileName = "platformFloat.cfg";
  std::string kernelName = "StreamingOCLMultiPlatform";
  bool collectStatistics = true;

  size_t dotPosition = scenarioFileName.find('.');
  std::string scenarioFileNamePrefix = scenarioFileName.substr(0, dotPosition);
  scenarioFileNamePrefix = scenarioFileNamePrefix.append("_" + kernelName);
  std::string outputFileName = scenarioFileNamePrefix + "_tuned.cfg";

  sgpp::datadriven::LearnerScenario scenario(scenarioBaseDir + scenarioFileName);
  sgpp::base::OCLOperationConfiguration parameter(parameterConfigurationFileName);
  sgpp::datadriven::StaticParameterTuner staticParameterTuner(parameter, true);

  std::string statisticsFolderName = "statistics";

  if (collectStatistics) {
    staticParameterTuner.enableStatistics(statisticsFolderName, scenarioFileNamePrefix);
    try {
      if (boost::filesystem::create_directory(statisticsFolderName)) {
        BOOST_TEST_MESSAGE("created output directory: " << statisticsFolderName);
      }
    } catch (boost::filesystem::filesystem_error &e) {
      BOOST_FAIL("could not create statistics output folder: " << statisticsFolderName << ": "
                                                               << e.what());
    }
  }

  staticParameterTuner.addParameter("KERNEL_USE_LOCAL_MEMORY", {"true", "false"});
  staticParameterTuner.addParameter("KERNEL_DATA_BLOCK_SIZE", {"1", "2", "4"});
  staticParameterTuner.addParameter("KERNEL_TRANS_GRID_BLOCK_SIZE", {"1", "2", "4"});
  staticParameterTuner.addParameter("KERNEL_STORE_DATA", {"array"});
  staticParameterTuner.addParameter("KERNEL_MAX_DIM_UNROLL", {"10", "1"});
  staticParameterTuner.addParameter("LOCAL_SIZE", {"128", "256"});
  staticParameterTuner.addParameter("VERBOSE", {"true"});
  staticParameterTuner.addParameter("OPTIMIZATION_FLAGS",
                                    {"-cl-strict-aliasing -cl-fast-relaxed-math"});

  sgpp::base::OCLOperationConfiguration bestParameters =
      staticParameterTuner.tuneEverything(scenario, kernelName);

  bestParameters.serialize(outputFileName);
}

BOOST_AUTO_TEST_CASE(Friedman2_4d_Linear_Double) {
  // internal precision is specified by the scenario, the parameter configuration is overwritten
  std::string scenarioFileName = "friedman2_4d_500000_Linear_double.scenario";
  std::string parameterConfigurationFileName = "platformDouble.cfg";
  std::string kernelName = "StreamingOCLMultiPlatform";
  bool collectStatistics = true;

  size_t dotPosition = scenarioFileName.find('.');
  std::string scenarioFileNamePrefix = scenarioFileName.substr(0, dotPosition);
  scenarioFileNamePrefix = scenarioFileNamePrefix.append("_" + kernelName);
  std::string outputFileName = scenarioFileNamePrefix + "_tuned.cfg";

  sgpp::datadriven::LearnerScenario scenario(scenarioBaseDir + scenarioFileName);
  sgpp::base::OCLOperationConfiguration parameter(parameterConfigurationFileName);
  sgpp::datadriven::StaticParameterTuner staticParameterTuner(parameter, true);

  std::string statisticsFolderName = "statistics";

  if (collectStatistics) {
    staticParameterTuner.enableStatistics(statisticsFolderName, scenarioFileNamePrefix);
    try {
      if (boost::filesystem::create_directory(statisticsFolderName)) {
        BOOST_TEST_MESSAGE("created output directory: " << statisticsFolderName);
      }
    } catch (boost::filesystem::filesystem_error &e) {
      BOOST_FAIL("could not create statistics output folder: " << statisticsFolderName << ": "
                                                               << e.what());
    }
  }

  staticParameterTuner.addParameter("KERNEL_USE_LOCAL_MEMORY", {"true", "false"});
  staticParameterTuner.addParameter("KERNEL_DATA_BLOCK_SIZE", {"1", "2", "4"});
  staticParameterTuner.addParameter("KERNEL_TRANS_GRID_BLOCK_SIZE", {"1", "2", "4"});
  staticParameterTuner.addParameter("KERNEL_STORE_DATA", {"array"});
  staticParameterTuner.addParameter("KERNEL_MAX_DIM_UNROLL", {"10", "1"});
  staticParameterTuner.addParameter("LOCAL_SIZE", {"128", "256"});
  staticParameterTuner.addParameter("VERBOSE", {"true"});
  staticParameterTuner.addParameter("OPTIMIZATION_FLAGS",
                                    {"-cl-strict-aliasing -cl-fast-relaxed-math"});

  sgpp::base::OCLOperationConfiguration bestParameters =
      staticParameterTuner.tuneEverything(scenario, kernelName);

  bestParameters.serialize(outputFileName);
}

BOOST_AUTO_TEST_CASE(Friedman2_4d_ModLinear_Float) {
  // internal precision is specified by the scenario, the parameter configuration is overwritten
  std::string scenarioFileName = "friedman2_4d_500000_ModLinear_float.scenario";
  std::string parameterConfigurationFileName = "platformFloat.cfg";
  std::string kernelName = "StreamingModOCLMaskMultiPlatform";
  bool collectStatistics = true;

  size_t dotPosition = scenarioFileName.find('.');
  std::string scenarioFileNamePrefix = scenarioFileName.substr(0, dotPosition);
  scenarioFileNamePrefix = scenarioFileNamePrefix.append("_" + kernelName);
  std::string outputFileName = scenarioFileNamePrefix + "_tuned.cfg";

  sgpp::datadriven::LearnerScenario scenario(scenarioBaseDir + scenarioFileName);
  sgpp::base::OCLOperationConfiguration parameter(parameterConfigurationFileName);
  sgpp::datadriven::StaticParameterTuner staticParameterTuner(parameter, true);

  std::string statisticsFolderName = "statistics";

  if (collectStatistics) {
    staticParameterTuner.enableStatistics(statisticsFolderName, scenarioFileNamePrefix);
    try {
      if (boost::filesystem::create_directory(statisticsFolderName)) {
        BOOST_TEST_MESSAGE("created output directory: " << statisticsFolderName);
      }
    } catch (boost::filesystem::filesystem_error &e) {
      BOOST_FAIL("could not create statistics output folder: " << statisticsFolderName << ": "
                                                               << e.what());
    }
  }

  staticParameterTuner.addParameter("KERNEL_USE_LOCAL_MEMORY", {"false", "true"});
  staticParameterTuner.addParameter("KERNEL_DATA_BLOCK_SIZE", {"1", "2", "4"});
  staticParameterTuner.addParameter("KERNEL_TRANS_GRID_BLOCK_SIZE", {"1", "2", "4"});
  staticParameterTuner.addParameter("KERNEL_STORE_DATA", {"array"});
  staticParameterTuner.addParameter("KERNEL_MAX_DIM_UNROLL", {"10", "1"});
  staticParameterTuner.addParameter("LOCAL_SIZE", {"128", "256"});
  staticParameterTuner.addParameter("VERBOSE", {"true"});
  staticParameterTuner.addParameter("OPTIMIZATION_FLAGS",
                                    {"-cl-strict-aliasing -cl-fast-relaxed-math"});

  sgpp::base::OCLOperationConfiguration bestParameters =
      staticParameterTuner.tuneEverything(scenario, kernelName);

  bestParameters.serialize(outputFileName);
}

BOOST_AUTO_TEST_CASE(Friedman2_4d_ModLinear_Double) {
  // internal precision is specified by the scenario, the parameter configuration is overwritten
  std::string scenarioFileName = "friedman2_4d_500000_ModLinear_double.scenario";
  std::string parameterConfigurationFileName = "platformDouble.cfg";
  std::string kernelName = "StreamingModOCLMaskMultiPlatform";
  bool collectStatistics = true;

  size_t dotPosition = scenarioFileName.find('.');
  std::string scenarioFileNamePrefix = scenarioFileName.substr(0, dotPosition);
  scenarioFileNamePrefix = scenarioFileNamePrefix.append("_" + kernelName);
  std::string outputFileName = scenarioFileNamePrefix + "_tuned.cfg";

  sgpp::datadriven::LearnerScenario scenario(scenarioBaseDir + scenarioFileName);
  sgpp::base::OCLOperationConfiguration parameter(parameterConfigurationFileName);
  sgpp::datadriven::StaticParameterTuner staticParameterTuner(parameter, true);

  std::string statisticsFolderName = "statistics";

  if (collectStatistics) {
    staticParameterTuner.enableStatistics(statisticsFolderName, scenarioFileNamePrefix);
    try {
      if (boost::filesystem::create_directory(statisticsFolderName)) {
        BOOST_TEST_MESSAGE("created output directory: " << statisticsFolderName);
      }
    } catch (boost::filesystem::filesystem_error &e) {
      BOOST_FAIL("could not create statistics output folder: " << statisticsFolderName << ": "
                                                               << e.what());
    }
  }

  staticParameterTuner.addParameter("KERNEL_USE_LOCAL_MEMORY", {"false", "true"});
  staticParameterTuner.addParameter("KERNEL_DATA_BLOCK_SIZE", {"1", "2", "4"});
  staticParameterTuner.addParameter("KERNEL_TRANS_GRID_BLOCK_SIZE", {"1", "2", "4"});
  staticParameterTuner.addParameter("KERNEL_STORE_DATA", {"array"});
  staticParameterTuner.addParameter("KERNEL_MAX_DIM_UNROLL", {"10", "1"});
  staticParameterTuner.addParameter("LOCAL_SIZE", {"128", "256"});
  staticParameterTuner.addParameter("VERBOSE", {"true"});
  staticParameterTuner.addParameter("OPTIMIZATION_FLAGS",
                                    {"-cl-strict-aliasing -cl-fast-relaxed-math"});

  sgpp::base::OCLOperationConfiguration bestParameters =
      staticParameterTuner.tuneEverything(scenario, kernelName);

  bestParameters.serialize(outputFileName);
}

/*
BOOST_AUTO_TEST_CASE(Friedman2_4d_ModLinearFast_Float) {
  // internal precision is specified by the scenario, the parameter configuration is overwritten
  std::string scenarioFileName =
      "friedman2_4d_500000_StreamingModOCLFastMultiPlatform_float.scenario";
  std::string parameterConfigurationFileName = "platformFloat.cfg";
  std::string kernelName = "StreamingModOCLFastMultiPlatform";
  bool collectStatistics = true;

  size_t dotPosition = scenarioFileName.find('.');
  std::string scenarioFileNamePrefix = scenarioFileName.substr(0, dotPosition);
  std::string outputFileName = scenarioFileNamePrefix + "_tuned.cfg";

  sgpp::datadriven::LearnerScenario scenario(scenarioBaseDir + scenarioFileName);
  sgpp::base::OCLOperationConfiguration parameter(parameterConfigurationFileName);
  sgpp::datadriven::StaticParameterTuner staticParameterTuner(parameter, true);

  std::string statisticsFolderName = "statistics";

  if (collectStatistics) {
    staticParameterTuner.enableStatistics(statisticsFolderName, scenarioFileNamePrefix);
    try {
      if (boost::filesystem::create_directory(statisticsFolderName)) {
        BOOST_TEST_MESSAGE("created output directory: " << statisticsFolderName);
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
      "friedman2_4d_500000_StreamingModOCLFastMultiPlatform_double.scenario";
  std::string parameterConfigurationFileName = "platformDouble.cfg";
  std::string kernelName = "StreamingModOCLFastMultiPlatform";
  bool collectStatistics = true;

  size_t dotPosition = scenarioFileName.find('.');
  std::string scenarioFileNamePrefix = scenarioFileName.substr(0, dotPosition);
  std::string outputFileName = scenarioFileNamePrefix + "_tuned.cfg";

  sgpp::datadriven::LearnerScenario scenario(scenarioBaseDir + scenarioFileName);
  sgpp::base::OCLOperationConfiguration parameter(parameterConfigurationFileName);
  sgpp::datadriven::StaticParameterTuner staticParameterTuner(parameter, true);

  std::string statisticsFolderName = "statistics";

  if (collectStatistics) {
    staticParameterTuner.enableStatistics(statisticsFolderName, scenarioFileNamePrefix);
    try {
      if (boost::filesystem::create_directory(statisticsFolderName)) {
        BOOST_TEST_MESSAGE("created output directory: " << statisticsFolderName);
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

BOOST_AUTO_TEST_CASE(Friedman2_4d_ModLinearOpt_Double) {
  // internal precision is specified by the scenario, the parameter configuration is overwritten
  std::string scenarioFileName = "friedman2_4d_500000_StreamingModOCLOpt_double.scenario";
  std::string parameterConfigurationFileName = "platformDouble.cfg";
  std::string kernelName = "StreamingModOCLOpt";
  bool collectStatistics = true;

  size_t dotPosition = scenarioFileName.find('.');
  std::string scenarioFileNamePrefix = scenarioFileName.substr(0, dotPosition);
  std::string outputFileName = scenarioFileNamePrefix + "_tuned.cfg";

  sgpp::datadriven::LearnerScenario scenario(scenarioBaseDir + scenarioFileName);
  sgpp::base::OCLOperationConfiguration parameter(parameterConfigurationFileName);
  sgpp::datadriven::StaticParameterTuner staticParameterTuner(parameter, true);

  std::string statisticsFolderName = "statistics";

  if (collectStatistics) {
    staticParameterTuner.enableStatistics(statisticsFolderName, scenarioFileNamePrefix);
    try {
      if (boost::filesystem::create_directory(statisticsFolderName)) {
        BOOST_TEST_MESSAGE("created output directory: " << statisticsFolderName);
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

BOOST_AUTO_TEST_CASE(Friedman2_4d_ModLinearOpt_Float) {
  // internal precision is specified by the scenario, the parameter configuration is overwritten
  std::string scenarioFileName = "friedman2_4d_500000_StreamingModOCLOpt_float.scenario";
  std::string parameterConfigurationFileName = "platformFloat.cfg";
  std::string kernelName = "StreamingModOCLOpt";
  bool collectStatistics = true;

  size_t dotPosition = scenarioFileName.find('.');
  std::string scenarioFileNamePrefix = scenarioFileName.substr(0, dotPosition);
  std::string outputFileName = scenarioFileNamePrefix + "_tuned.cfg";

  sgpp::datadriven::LearnerScenario scenario(scenarioBaseDir + scenarioFileName);
  sgpp::base::OCLOperationConfiguration parameter(parameterConfigurationFileName);
  sgpp::datadriven::StaticParameterTuner staticParameterTuner(parameter, true);

  std::string statisticsFolderName = "statistics";

  if (collectStatistics) {
    staticParameterTuner.enableStatistics(statisticsFolderName, scenarioFileNamePrefix);
    try {
      if (boost::filesystem::create_directory(statisticsFolderName)) {
        BOOST_TEST_MESSAGE("created output directory: " << statisticsFolderName);
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
*/

BOOST_AUTO_TEST_CASE(Friedman1_10d_Linear_Float) {
  // internal precision is specified by the scenario, the parameter configuration is overwritten
  std::string scenarioFileName = "friedman1_10d_500000_Linear_float.scenario";
  std::string parameterConfigurationFileName = "platformFloat.cfg";
  std::string kernelName = "StreamingOCLMultiPlatform";
  bool collectStatistics = true;

  size_t dotPosition = scenarioFileName.find('.');
  std::string scenarioFileNamePrefix = scenarioFileName.substr(0, dotPosition);
  scenarioFileNamePrefix = scenarioFileNamePrefix.append("_" + kernelName);
  std::string outputFileName = scenarioFileNamePrefix + "_tuned.cfg";

  sgpp::datadriven::LearnerScenario scenario(scenarioBaseDir + scenarioFileName);
  sgpp::base::OCLOperationConfiguration parameter(parameterConfigurationFileName);
  sgpp::datadriven::StaticParameterTuner staticParameterTuner(parameter, true);

  std::string statisticsFolderName = "statistics";

  if (collectStatistics) {
    staticParameterTuner.enableStatistics(statisticsFolderName, scenarioFileNamePrefix);
    try {
      if (boost::filesystem::create_directory(statisticsFolderName)) {
        BOOST_TEST_MESSAGE("created output directory: " << statisticsFolderName);
      }
    } catch (boost::filesystem::filesystem_error &e) {
      BOOST_FAIL("could not create statistics output folder: " << statisticsFolderName << ": "
                                                               << e.what());
    }
  }

  staticParameterTuner.addParameter("KERNEL_USE_LOCAL_MEMORY", {"true", "false"});
  staticParameterTuner.addParameter("KERNEL_DATA_BLOCK_SIZE", {"1", "2", "4"});
  staticParameterTuner.addParameter("KERNEL_TRANS_GRID_BLOCK_SIZE", {"1", "2", "4"});
  staticParameterTuner.addParameter("KERNEL_STORE_DATA", {"array"});
  staticParameterTuner.addParameter("KERNEL_MAX_DIM_UNROLL", {"10", "4", "1"});
  staticParameterTuner.addParameter("LOCAL_SIZE", {"128", "256"});
  staticParameterTuner.addParameter("VERBOSE", {"true"});
  staticParameterTuner.addParameter("KERNEL_PREFETCH_SIZE", {"32", "64"});
  staticParameterTuner.addParameter("KERNEL_TRANS_PREFETCH_SIZE", {"32", "64"});
  staticParameterTuner.addParameter("OPTIMIZATION_FLAGS",
                                    {"-cl-strict-aliasing -cl-fast-relaxed-math"});

  sgpp::base::OCLOperationConfiguration bestParameters =
      staticParameterTuner.tuneEverything(scenario, kernelName);

  bestParameters.serialize(outputFileName);
}

BOOST_AUTO_TEST_CASE(Friedman1_10d_Linear_Double) {
  // internal precision is specified by the scenario, the parameter configuration is overwritten
  std::string scenarioFileName = "friedman1_10d_500000_Linear_double.scenario";
  std::string parameterConfigurationFileName = "platformDouble.cfg";
  std::string kernelName = "StreamingOCLMultiPlatform";
  bool collectStatistics = true;

  size_t dotPosition = scenarioFileName.find('.');
  std::string scenarioFileNamePrefix = scenarioFileName.substr(0, dotPosition);
  scenarioFileNamePrefix = scenarioFileNamePrefix.append("_" + kernelName);
  std::string outputFileName = scenarioFileNamePrefix + "_tuned.cfg";

  sgpp::datadriven::LearnerScenario scenario(scenarioBaseDir + scenarioFileName);
  sgpp::base::OCLOperationConfiguration parameter(parameterConfigurationFileName);
  sgpp::datadriven::StaticParameterTuner staticParameterTuner(parameter, true);

  std::string statisticsFolderName = "statistics";

  if (collectStatistics) {
    staticParameterTuner.enableStatistics(statisticsFolderName, scenarioFileNamePrefix);
    try {
      if (boost::filesystem::create_directory(statisticsFolderName)) {
        BOOST_TEST_MESSAGE("created output directory: " << statisticsFolderName);
      }
    } catch (boost::filesystem::filesystem_error &e) {
      BOOST_FAIL("could not create statistics output folder: " << statisticsFolderName << ": "
                                                               << e.what());
    }
  }

  staticParameterTuner.addParameter("KERNEL_USE_LOCAL_MEMORY", {"true", "false"});
  staticParameterTuner.addParameter("KERNEL_DATA_BLOCK_SIZE", {"1", "2", "4"});
  staticParameterTuner.addParameter("KERNEL_TRANS_GRID_BLOCK_SIZE", {"1", "2", "4"});
  staticParameterTuner.addParameter("KERNEL_STORE_DATA", {"array"});
  staticParameterTuner.addParameter("KERNEL_MAX_DIM_UNROLL", {"10", "4", "1"});
  staticParameterTuner.addParameter("LOCAL_SIZE", {"128", "256"});
  staticParameterTuner.addParameter("VERBOSE", {"true"});
  staticParameterTuner.addParameter("KERNEL_PREFETCH_SIZE", {"32", "64"});
  staticParameterTuner.addParameter("KERNEL_TRANS_PREFETCH_SIZE", {"32", "64"});
  staticParameterTuner.addParameter("OPTIMIZATION_FLAGS",
                                    {"-cl-strict-aliasing -cl-fast-relaxed-math"});

  sgpp::base::OCLOperationConfiguration bestParameters =
      staticParameterTuner.tuneEverything(scenario, kernelName);

  bestParameters.serialize(outputFileName);
}

BOOST_AUTO_TEST_CASE(Friedman1_10d_ModLinear_Float) {
  // internal precision is specified by the scenario, the parameter configuration is overwritten
  std::string scenarioFileName = "friedman1_10d_500000_ModLinear_float.scenario";
  std::string parameterConfigurationFileName = "platformFloat.cfg";
  std::string kernelName = "StreamingModOCLMaskMultiPlatform";
  bool collectStatistics = true;

  size_t dotPosition = scenarioFileName.find('.');
  std::string scenarioFileNamePrefix = scenarioFileName.substr(0, dotPosition);
  scenarioFileNamePrefix = scenarioFileNamePrefix.append("_" + kernelName);
  std::string outputFileName = scenarioFileNamePrefix + "_tuned.cfg";

  sgpp::datadriven::LearnerScenario scenario(scenarioBaseDir + scenarioFileName);
  sgpp::base::OCLOperationConfiguration parameter(parameterConfigurationFileName);
  sgpp::datadriven::StaticParameterTuner staticParameterTuner(parameter, true);

  std::string statisticsFolderName = "statistics";

  if (collectStatistics) {
    staticParameterTuner.enableStatistics(statisticsFolderName, scenarioFileNamePrefix);
    try {
      if (boost::filesystem::create_directory(statisticsFolderName)) {
        BOOST_TEST_MESSAGE("created output directory: " << statisticsFolderName);
      }
    } catch (boost::filesystem::filesystem_error &e) {
      BOOST_FAIL("could not create statistics output folder: " << statisticsFolderName << ": "
                                                               << e.what());
    }
  }

  staticParameterTuner.addParameter("KERNEL_USE_LOCAL_MEMORY", {"false", "true"});
  staticParameterTuner.addParameter("KERNEL_DATA_BLOCK_SIZE", {"1", "2", "4"});
  staticParameterTuner.addParameter("KERNEL_TRANS_GRID_BLOCK_SIZE", {"1", "2", "4"});
  staticParameterTuner.addParameter("KERNEL_STORE_DATA", {"array"});
  staticParameterTuner.addParameter("KERNEL_MAX_DIM_UNROLL", {"10", "4", "1"});
  staticParameterTuner.addParameter("LOCAL_SIZE", {"128", "256"});
  staticParameterTuner.addParameter("VERBOSE", {"true"});
  staticParameterTuner.addParameter("KERNEL_PREFETCH_SIZE", {"32", "64"});
  staticParameterTuner.addParameter("KERNEL_TRANS_PREFETCH_SIZE", {"32", "64"});
  staticParameterTuner.addParameter("OPTIMIZATION_FLAGS",
                                    {"-cl-strict-aliasing -cl-fast-relaxed-math"});

  sgpp::base::OCLOperationConfiguration bestParameters =
      staticParameterTuner.tuneEverything(scenario, kernelName);

  bestParameters.serialize(outputFileName);
}

BOOST_AUTO_TEST_CASE(Friedman1_10d_ModLinear_Double) {
  // internal precision is specified by the scenario, the parameter configuration is overwritten
  std::string scenarioFileName = "friedman1_10d_500000_ModLinear_double.scenario";
  std::string parameterConfigurationFileName = "platformDouble.cfg";
  std::string kernelName = "StreamingModOCLMaskMultiPlatform";
  bool collectStatistics = true;

  size_t dotPosition = scenarioFileName.find('.');
  std::string scenarioFileNamePrefix = scenarioFileName.substr(0, dotPosition);
  scenarioFileNamePrefix = scenarioFileNamePrefix.append("_" + kernelName);
  std::string outputFileName = scenarioFileNamePrefix + "_tuned.cfg";

  sgpp::datadriven::LearnerScenario scenario(scenarioBaseDir + scenarioFileName);
  sgpp::base::OCLOperationConfiguration parameter(parameterConfigurationFileName);
  sgpp::datadriven::StaticParameterTuner staticParameterTuner(parameter, true);

  std::string statisticsFolderName = "statistics";

  if (collectStatistics) {
    staticParameterTuner.enableStatistics(statisticsFolderName, scenarioFileNamePrefix);
    try {
      if (boost::filesystem::create_directory(statisticsFolderName)) {
        BOOST_TEST_MESSAGE("created output directory: " << statisticsFolderName);
      }
    } catch (boost::filesystem::filesystem_error &e) {
      BOOST_FAIL("could not create statistics output folder: " << statisticsFolderName << ": "
                                                               << e.what());
    }
  }

  staticParameterTuner.addParameter("KERNEL_USE_LOCAL_MEMORY", {"false", "true"});
  staticParameterTuner.addParameter("KERNEL_DATA_BLOCK_SIZE", {"1", "2", "4"});
  staticParameterTuner.addParameter("KERNEL_TRANS_GRID_BLOCK_SIZE", {"1", "2", "4"});
  staticParameterTuner.addParameter("KERNEL_STORE_DATA", {"array"});
  staticParameterTuner.addParameter("KERNEL_MAX_DIM_UNROLL", {"10", "4", "1"});
  staticParameterTuner.addParameter("LOCAL_SIZE", {"128", "256"});
  staticParameterTuner.addParameter("VERBOSE", {"true"});
  staticParameterTuner.addParameter("KERNEL_PREFETCH_SIZE", {"32", "64"});
  staticParameterTuner.addParameter("KERNEL_TRANS_PREFETCH_SIZE", {"32", "64"});
  staticParameterTuner.addParameter("OPTIMIZATION_FLAGS",
                                    {"-cl-strict-aliasing -cl-fast-relaxed-math"});

  sgpp::base::OCLOperationConfiguration bestParameters =
      staticParameterTuner.tuneEverything(scenario, kernelName);

  bestParameters.serialize(outputFileName);
}

BOOST_AUTO_TEST_CASE(DR5_Linear_Float) {
  // internal precision is specified by the scenario, the parameter configuration is overwritten
  std::string scenarioFileName = "DR5_train_Linear_float.scenario";
  std::string parameterConfigurationFileName = "platformFloat.cfg";
  std::string kernelName = "StreamingOCLMultiPlatform";
  bool collectStatistics = true;

  size_t dotPosition = scenarioFileName.find('.');
  std::string scenarioFileNamePrefix = scenarioFileName.substr(0, dotPosition);
  scenarioFileNamePrefix = scenarioFileNamePrefix.append("_" + kernelName);
  std::string outputFileName = scenarioFileNamePrefix + "_tuned.cfg";

  sgpp::datadriven::LearnerScenario scenario(scenarioBaseDir + scenarioFileName);
  sgpp::base::OCLOperationConfiguration parameter(parameterConfigurationFileName);
  sgpp::datadriven::StaticParameterTuner staticParameterTuner(parameter, true);

  std::string statisticsFolderName = "statistics";

  if (collectStatistics) {
    staticParameterTuner.enableStatistics(statisticsFolderName, scenarioFileNamePrefix);
    try {
      if (boost::filesystem::create_directory(statisticsFolderName)) {
        BOOST_TEST_MESSAGE("created output directory: " << statisticsFolderName);
      }
    } catch (boost::filesystem::filesystem_error &e) {
      BOOST_FAIL("could not create statistics output folder: " << statisticsFolderName << ": "
                                                               << e.what());
    }
  }

  staticParameterTuner.addParameter("KERNEL_USE_LOCAL_MEMORY", {"true", "false"});
  staticParameterTuner.addParameter("KERNEL_DATA_BLOCK_SIZE", {"1", "2", "4"});
  staticParameterTuner.addParameter("KERNEL_TRANS_GRID_BLOCK_SIZE", {"1", "2", "4"});
  staticParameterTuner.addParameter("KERNEL_STORE_DATA", {"array"});
  staticParameterTuner.addParameter("KERNEL_MAX_DIM_UNROLL", {"10", "4", "1"});
  staticParameterTuner.addParameter("LOCAL_SIZE", {"128", "256"});
  staticParameterTuner.addParameter("VERBOSE", {"true"});
  staticParameterTuner.addParameter("KERNEL_PREFETCH_SIZE", {"64"});
  staticParameterTuner.addParameter("KERNEL_TRANS_PREFETCH_SIZE", {"64"});
  staticParameterTuner.addParameter("OPTIMIZATION_FLAGS",
                                    {"-cl-strict-aliasing -cl-fast-relaxed-math"});

  sgpp::base::OCLOperationConfiguration bestParameters =
      staticParameterTuner.tuneEverything(scenario, kernelName);

  bestParameters.serialize(outputFileName);
}

BOOST_AUTO_TEST_CASE(DR5_Linear_Double) {
  // internal precision is specified by the scenario, the parameter configuration is overwritten
  std::string scenarioFileName = "DR5_train_Linear_double.scenario";
  std::string parameterConfigurationFileName = "platformDouble.cfg";
  std::string kernelName = "StreamingOCLMultiPlatform";
  bool collectStatistics = true;

  size_t dotPosition = scenarioFileName.find('.');
  std::string scenarioFileNamePrefix = scenarioFileName.substr(0, dotPosition);
  scenarioFileNamePrefix = scenarioFileNamePrefix.append("_" + kernelName);
  std::string outputFileName = scenarioFileNamePrefix + "_tuned.cfg";

  sgpp::datadriven::LearnerScenario scenario(scenarioBaseDir + scenarioFileName);
  sgpp::base::OCLOperationConfiguration parameter(parameterConfigurationFileName);
  sgpp::datadriven::StaticParameterTuner staticParameterTuner(parameter, true);

  std::string statisticsFolderName = "statistics";

  if (collectStatistics) {
    staticParameterTuner.enableStatistics(statisticsFolderName, scenarioFileNamePrefix);
    try {
      if (boost::filesystem::create_directory(statisticsFolderName)) {
        BOOST_TEST_MESSAGE("created output directory: " << statisticsFolderName);
      }
    } catch (boost::filesystem::filesystem_error &e) {
      BOOST_FAIL("could not create statistics output folder: " << statisticsFolderName << ": "
                                                               << e.what());
    }
  }

  staticParameterTuner.addParameter("KERNEL_USE_LOCAL_MEMORY", {"true", "false"});
  staticParameterTuner.addParameter("KERNEL_DATA_BLOCK_SIZE", {"1", "2", "4"});
  staticParameterTuner.addParameter("KERNEL_TRANS_GRID_BLOCK_SIZE", {"1", "2", "4"});
  staticParameterTuner.addParameter("KERNEL_STORE_DATA", {"array"});
  staticParameterTuner.addParameter("KERNEL_MAX_DIM_UNROLL", {"10", "4", "1"});
  staticParameterTuner.addParameter("LOCAL_SIZE", {"128", "256"});
  staticParameterTuner.addParameter("VERBOSE", {"true"});
  staticParameterTuner.addParameter("KERNEL_PREFETCH_SIZE", {"64"});
  staticParameterTuner.addParameter("KERNEL_TRANS_PREFETCH_SIZE", {"64"});
  staticParameterTuner.addParameter("OPTIMIZATION_FLAGS",
                                    {"-cl-strict-aliasing -cl-fast-relaxed-math"});

  sgpp::base::OCLOperationConfiguration bestParameters =
      staticParameterTuner.tuneEverything(scenario, kernelName);

  bestParameters.serialize(outputFileName);
}

BOOST_AUTO_TEST_CASE(DR5_ModLinear_Float) {
  // internal precision is specified by the scenario, the parameter configuration is overwritten
  std::string scenarioFileName = "DR5_train_ModLinear_float.scenario";
  std::string parameterConfigurationFileName = "platformFloat.cfg";
  std::string kernelName = "StreamingModOCLMaskMultiPlatform";
  bool collectStatistics = true;

  size_t dotPosition = scenarioFileName.find('.');
  std::string scenarioFileNamePrefix = scenarioFileName.substr(0, dotPosition);
  scenarioFileNamePrefix = scenarioFileNamePrefix.append("_" + kernelName);
  std::string outputFileName = scenarioFileNamePrefix + "_tuned.cfg";

  sgpp::datadriven::LearnerScenario scenario(scenarioBaseDir + scenarioFileName);
  sgpp::base::OCLOperationConfiguration parameter(parameterConfigurationFileName);
  sgpp::datadriven::StaticParameterTuner staticParameterTuner(parameter, true);

  std::string statisticsFolderName = "statistics";

  if (collectStatistics) {
    staticParameterTuner.enableStatistics(statisticsFolderName, scenarioFileNamePrefix);
    try {
      if (boost::filesystem::create_directory(statisticsFolderName)) {
        BOOST_TEST_MESSAGE("created output directory: " << statisticsFolderName);
      }
    } catch (boost::filesystem::filesystem_error &e) {
      BOOST_FAIL("could not create statistics output folder: " << statisticsFolderName << ": "
                                                               << e.what());
    }
  }

  staticParameterTuner.addParameter("KERNEL_USE_LOCAL_MEMORY", {"false", "true"});
  staticParameterTuner.addParameter("KERNEL_DATA_BLOCK_SIZE", {"1", "2", "4"});
  staticParameterTuner.addParameter("KERNEL_TRANS_GRID_BLOCK_SIZE", {"1", "2", "4"});
  staticParameterTuner.addParameter("KERNEL_STORE_DATA", {"array"});
  staticParameterTuner.addParameter("KERNEL_MAX_DIM_UNROLL", {"10", "4", "1"});
  staticParameterTuner.addParameter("LOCAL_SIZE", {"128", "256"});
  staticParameterTuner.addParameter("VERBOSE", {"true"});
  staticParameterTuner.addParameter("KERNEL_PREFETCH_SIZE", {"64"});
  staticParameterTuner.addParameter("KERNEL_TRANS_PREFETCH_SIZE", {"64"});
  staticParameterTuner.addParameter("OPTIMIZATION_FLAGS",
                                    {"-cl-strict-aliasing -cl-fast-relaxed-math"});

  sgpp::base::OCLOperationConfiguration bestParameters =
      staticParameterTuner.tuneEverything(scenario, kernelName);

  bestParameters.serialize(outputFileName);
}

BOOST_AUTO_TEST_CASE(DR5_ModLinear_Double) {
  // internal precision is specified by the scenario, the parameter configuration is overwritten
  std::string scenarioFileName = "DR5_train_ModLinear_double.scenario";
  std::string parameterConfigurationFileName = "platformDouble.cfg";
  std::string kernelName = "StreamingModOCLMaskMultiPlatform";
  bool collectStatistics = true;

  size_t dotPosition = scenarioFileName.find('.');
  std::string scenarioFileNamePrefix = scenarioFileName.substr(0, dotPosition);
  scenarioFileNamePrefix = scenarioFileNamePrefix.append("_" + kernelName);
  std::string outputFileName = scenarioFileNamePrefix + "_tuned.cfg";

  sgpp::datadriven::LearnerScenario scenario(scenarioBaseDir + scenarioFileName);
  sgpp::base::OCLOperationConfiguration parameter(parameterConfigurationFileName);
  sgpp::datadriven::StaticParameterTuner staticParameterTuner(parameter, true);

  std::string statisticsFolderName = "statistics";

  if (collectStatistics) {
    staticParameterTuner.enableStatistics(statisticsFolderName, scenarioFileNamePrefix);
    try {
      if (boost::filesystem::create_directory(statisticsFolderName)) {
        BOOST_TEST_MESSAGE("created output directory: " << statisticsFolderName);
      }
    } catch (boost::filesystem::filesystem_error &e) {
      BOOST_FAIL("could not create statistics output folder: " << statisticsFolderName << ": "
                                                               << e.what());
    }
  }

  staticParameterTuner.addParameter("KERNEL_USE_LOCAL_MEMORY", {"false", "true"});
  staticParameterTuner.addParameter("KERNEL_DATA_BLOCK_SIZE", {"1", "2", "4"});
  staticParameterTuner.addParameter("KERNEL_TRANS_GRID_BLOCK_SIZE", {"1", "2", "4"});
  staticParameterTuner.addParameter("KERNEL_STORE_DATA", {"array"});
  staticParameterTuner.addParameter("KERNEL_MAX_DIM_UNROLL", {"10", "4", "1"});
  staticParameterTuner.addParameter("LOCAL_SIZE", {"128", "256"});
  staticParameterTuner.addParameter("VERBOSE", {"true"});
  staticParameterTuner.addParameter("KERNEL_PREFETCH_SIZE", {"64"});
  staticParameterTuner.addParameter("KERNEL_TRANS_PREFETCH_SIZE", {"64"});
  staticParameterTuner.addParameter("OPTIMIZATION_FLAGS",
                                    {"-cl-strict-aliasing -cl-fast-relaxed-math"});

  sgpp::base::OCLOperationConfiguration bestParameters =
      staticParameterTuner.tuneEverything(scenario, kernelName);

  bestParameters.serialize(outputFileName);
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(AutoTuningPaperIntrinsicsComparison)

BOOST_AUTO_TEST_CASE(StreamingIntrinsicsComparison) {
  std::string baseFolder = "datadriven/performanceTests/scenarios/";
  std::vector<std::string> scenarios = {
      baseFolder + "chess_4d_500000_Linear_double.scenario",
      baseFolder + "chess_4d_500000_Linear_float.scenario",
      baseFolder + "chess_4d_500000_ModLinear_double.scenario",
      baseFolder + "chess_4d_500000_ModLinear_float.scenario",
      baseFolder + "friedman2_4d_500000_Linear_double.scenario",
      baseFolder + "friedman2_4d_500000_Linear_float.scenario",
      baseFolder + "friedman2_4d_500000_ModLinear_double.scenario",
      baseFolder + "friedman2_4d_500000_ModLinear_float.scenario",
      baseFolder + "friedman1_10d_500000_Linear_double.scenario",
      baseFolder + "friedman1_10d_500000_Linear_float.scenario",
      baseFolder + "friedman1_10d_500000_ModLinear_double.scenario",
      baseFolder + "friedman1_10d_500000_ModLinear_float.scenario",
      baseFolder + "DR5_train_Linear_double.scenario",
      baseFolder + "DR5_train_Linear_float.scenario",
      baseFolder + "DR5_train_ModLinear_double.scenario",
      baseFolder + "DR5_train_ModLinear_float.scenario"};

  std::ofstream outFile("statistics/intrinsics.log");
  for (size_t i = 0; i < scenarios.size(); i++) {
    std::string &scenarioFileName = scenarios[i];

    std::cout << "scenario: " << scenarioFileName << std::endl;

    outFile << scenarioFileName << std::endl;

    sgpp::datadriven::LearnerScenario scenario(scenarioFileName);

    sgpp::datadriven::OperationMultipleEvalConfiguration configuration(
        sgpp::datadriven::OperationMultipleEvalType::STREAMING,
        sgpp::datadriven::OperationMultipleEvalSubType::DEFAULT);

    bool verbose = true;
    sgpp::datadriven::MetaLearner learner(
        scenario.getGridConfig(), scenario.getSolverConfigurationRefine(),
        scenario.getSolverConfigurationFinal(), scenario.getAdaptivityConfiguration(),
        scenario.getLambda(), verbose);

    std::string datasetFile = scenario.getDatasetFileName();
    try {
      learner.learn(configuration, datasetFile);
      //  learner.learnAndCompare(configuration, datasetFile, 4);

      sgpp::datadriven::LearnerTiming timing = learner.getLearnerTiming();
      std::cout << "time complete: " << timing.timeComplete_ << std::endl;
      outFile << timing.timeComplete_ << std::endl;
    } catch (sgpp::base::operation_exception &e) {
      std::cout << "exception caught: " << e.what() << std::endl;
    }
  }
  outFile.close();
}

BOOST_AUTO_TEST_SUITE_END()

#endif
#endif
