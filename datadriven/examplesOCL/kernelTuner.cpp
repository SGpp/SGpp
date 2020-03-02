// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/opencl/OCLOperationConfiguration.hpp>
#include <sgpp/datadriven/application/LearnerScenario.hpp>
#include <sgpp/datadriven/application/MetaLearner.hpp>
#include <sgpp/datadriven/application/StaticParameterTuner.hpp>
#include <sgpp/datadriven/operation/hash/DatadrivenOperationCommon.hpp>

#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>

#include <iostream>
#include <string>


int main(int argc, char **argv) {
  std::string scenarioFileName;
  std::string parameterConfigurationFileName;  // <- includes kernel type, subtype
  std::string kernelName;
  std::string outputFileName;
  bool collectStatistics = false;
  //    bool useDoublePrecision = false;
  //    std::vector<std::string> devices;

  boost::program_options::options_description description("Allowed options");
  description.add_options()("help", "display help")(
      "scenario", boost::program_options::value<std::string>(&scenarioFileName),
      "the scenario file to be used (serialized LearnerScenario)")(
      "parameterConfiguration",
      boost::program_options::value<std::string>(&parameterConfigurationFileName),
      "the parameter configuration file to be used (serialized "
      "StaticParameterTuner)")("kernel", boost::program_options::value<std::string>(&kernelName),
                               "name of the kernel to be tuned for")(
      "outputFileName", boost::program_options::value<std::string>(&outputFileName),
      "output file for optimized parameters")(
      "collectStatistics", boost::program_options::value<bool>(&collectStatistics),
      "collect statistics for each device and kernel optimized and write them to csv files");

  boost::program_options::variables_map variables_map;

  boost::program_options::parsed_options options = parse_command_line(argc, argv, description);
  boost::program_options::store(options, variables_map);
  boost::program_options::notify(variables_map);

  if (variables_map.count("help")) {
    std::cout << description << std::endl;
    return 0;
  }

  // check whether all files exist
  std::ifstream scenarioFile(scenarioFileName);
  if (!scenarioFile.good()) {
    std::cout << "scenario file not found" << std::endl;
    return 1;
  }
  scenarioFile.close();
  std::ifstream parameterConfigurationFile(parameterConfigurationFileName);
  if (!parameterConfigurationFile.good()) {
    std::cout << "parameter file not found" << std::endl;
    return 1;
  }
  parameterConfigurationFile.close();

  sgpp::datadriven::LearnerScenario scenario(scenarioFileName);

  sgpp::base::OCLOperationConfiguration parameter(parameterConfigurationFileName);

  sgpp::datadriven::StaticParameterTuner staticParameterTuner(parameter, true);

  size_t dotPos = scenarioFileName.find('.');
  std::string scenarioFileNamePrefix = scenarioFileName.substr(0, dotPos);
  std::string statisticsFolderName = "statistics";

  if (collectStatistics) {
    staticParameterTuner.enableStatistics(statisticsFolderName, scenarioFileNamePrefix);

    if (boost::filesystem::exists(boost::filesystem::path(statisticsFolderName))) {
      std::cerr << "error: statistics output directory already exists: " << statisticsFolderName
                << std::endl;
      return 1;
    }

    try {
      if (boost::filesystem::create_directory(statisticsFolderName)) {
        std::cout << "created output directory: " << statisticsFolderName << std::endl;
      }
    } catch (boost::filesystem::filesystem_error &e) {
      std::cerr << "could not create statistics output folder: " << statisticsFolderName << ": "
                << e.what() << std::endl;
      return 1;
    }
  }

  if (kernelName.compare("StreamingOCLMultiPlatform") == 0) {
    staticParameterTuner.addParameter("KERNEL_USE_LOCAL_MEMORY", {"true", "false"});
    staticParameterTuner.addParameter("KERNEL_DATA_BLOCK_SIZE", {"1", "2", "4"});
    staticParameterTuner.addParameter("KERNEL_TRANS_GRID_BLOCK_SIZE", {"1", "2", "4"});
    staticParameterTuner.addParameter("KERNEL_STORE_DATA", {"register", "array"});
    staticParameterTuner.addParameter("KERNEL_MAX_DIM_UNROLL", {"10", "1"});
    staticParameterTuner.addParameter("LOCAL_SIZE", {"128", "256"});
  } else if (kernelName.compare("StreamingModOCLFastMultiPlatform") == 0) {
    staticParameterTuner.addParameter("KERNEL_USE_LOCAL_MEMORY", {"false", "true"});
    staticParameterTuner.addParameter("KERNEL_DATA_BLOCK_SIZE", {"1", "2", "4"});
    staticParameterTuner.addParameter("KERNEL_TRANS_DATA_BLOCK_SIZE", {"1", "2", "4"});
    staticParameterTuner.addParameter("KERNEL_TRANS_GRID_BLOCK_SIZE", {"1", "2", "4"});
    staticParameterTuner.addParameter("KERNEL_STORE_DATA", {"register", "array"});
    staticParameterTuner.addParameter("KERNEL_MAX_DIM_UNROLL", {"10", "1"});
    staticParameterTuner.addParameter("LOCAL_SIZE", {"128", "256"});
  } else if (kernelName.compare("StreamingModOCLMaskMultiPlatform") == 0) {
    staticParameterTuner.addParameter("KERNEL_USE_LOCAL_MEMORY", {"false", "true"});
    staticParameterTuner.addParameter("KERNEL_DATA_BLOCK_SIZE", {"1", "2", "4", "8"});
    staticParameterTuner.addParameter("KERNEL_TRANS_GRID_BLOCK_SIZE", {"1", "2", "4"});
    staticParameterTuner.addParameter("KERNEL_STORE_DATA", {"register", "array"});
    staticParameterTuner.addParameter("KERNEL_MAX_DIM_UNROLL", {"10", "1"});
    staticParameterTuner.addParameter("LOCAL_SIZE", {"128", "256"});
  } else {
    std::cout << "error: kernel name \"" << kernelName << "\" not recognized" << std::endl;
    return 1;
  }

  //  staticParameterTuner.setupReferenceValues(referenceDataset, referenceValues, expectedMSE,
  //                                            expectedLargestDifference);

  sgpp::base::OCLOperationConfiguration bestParameters =
      staticParameterTuner.tuneEverything(scenario, kernelName);

  bestParameters.serialize(outputFileName);

  std::cout << "-------------- all done! --------------" << std::endl;
  return 0;
}
