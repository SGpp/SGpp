// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <iostream>
#include <string>

#include "sgpp/datadriven/application/MetaLearner.hpp"
#include "sgpp/datadriven/operation/hash/simple/DatadrivenOperationCommon.hpp"
#include "sgpp/base/opencl/OCLOperationConfiguration.hpp"
#include "sgpp/datadriven/application/StaticParameterTuner.hpp"
#include "sgpp/datadriven/application/TunableParameter.hpp"
#include "sgpp/datadriven/application/LearnerScenario.hpp"

int main(int argc, char** argv) {
  // Specify scenario for performance optimization

  int maxLevel = 10;
  //    std::string fileName = "friedman_4d_small.arff";
  std::string fileName = "friedman_4d.arff";
  SGPP::base::RegularGridConfiguration gridConfig;
  SGPP::solver::SLESolverConfiguration SLESolverConfigRefine;
  SGPP::solver::SLESolverConfiguration SLESolverConfigFinal;
  SGPP::base::AdpativityConfiguration adaptConfig;

  // setup grid
  gridConfig.dim_ = 0;  // dim is inferred from the data
  gridConfig.level_ = maxLevel;
  gridConfig.type_ = SGPP::base::GridType::Linear;

  // Set Adaptivity
  adaptConfig.maxLevelType_ = false;
  adaptConfig.noPoints_ = 80;
  adaptConfig.numRefinements_ = 0;
  adaptConfig.percent_ = 200.0;
  adaptConfig.threshold_ = 0.0;

  // Set solver during refinement
  SLESolverConfigRefine.eps_ = 0;
  SLESolverConfigRefine.maxIterations_ = 1;
  SLESolverConfigRefine.threshold_ = -1.0;
  SLESolverConfigRefine.type_ = SGPP::solver::SLESolverType::CG;

  // Set solver for final step
  SLESolverConfigFinal.eps_ = 0;
  SLESolverConfigFinal.maxIterations_ = 1;
  SLESolverConfigFinal.threshold_ = -1.0;
  SLESolverConfigFinal.type_ = SGPP::solver::SLESolverType::CG;

  double lambda = 0.000001;

  SGPP::datadriven::LearnerScenario scenario(fileName, lambda, gridConfig, SLESolverConfigRefine,
                                             SLESolverConfigFinal, adaptConfig);

  std::string kernelName("StreamingOCLMultiPlatform");

  SGPP::base::OCLOperationConfiguration parameter("skeleton.cfg");
  //    json::Node &kernelNode =
  //    parameter["PLATFORMS"][platformName]["DEVICES"][deviceName]["KERNELS]"][kernelName];

  SGPP::datadriven::StaticParameterTuner staticParameterTuner(parameter, true, true);

  //    staticParameterTuner.addFixedParameter("OCL_MANAGER_VERBOSE", "false",
  //    SGPP::datadriven::ParameterType::ID);
  //    staticParameterTuner.addFixedParameter("VERBOSE", "false",
  //    SGPP::datadriven::ParameterType::ID);
  //    staticParameterTuner.addFixedParameter("PLATFORM", "NVIDIA CUDA",
  //    SGPP::datadriven::ParameterType::TEXT);
  //    staticParameterTuner.addFixedParameter("INTERNAL_PRECISION", "float",
  //    SGPP::datadriven::ParameterType::TEXT);

  staticParameterTuner.addParameter("KERNEL_USE_LOCAL_MEMORY", {"false"});        // , "true"
  staticParameterTuner.addParameter("KERNEL_DATA_BLOCK_SIZE", {"1", "4"});        // , "2", "4", "8"
  staticParameterTuner.addParameter("KERNEL_TRANS_GRID_BLOCK_SIZE", {"1", "4"});  // , "2", "4", "8"
  staticParameterTuner.addParameter("KERNEL_STORE_DATA", {"register"});           // "array",
  staticParameterTuner.addParameter("KERNEL_MAX_DIM_UNROLL", {"4"});              // "1", "8", "16"

  //    staticParameterTuner.writeToFile("testOut.tuner");
  //
  //    SGPP::datadriven::StaticParameterTuner
  //    readStaticParameterTuner("testOut.tuner");

  SGPP::base::OCLOperationConfiguration bestParameters =
      staticParameterTuner.tuneEverything(scenario, kernelName);

  //    json::Node &kernelNode =
  //    parameter["PLATFORMS"][platformName]["DEVICES"][deviceName]["KERNELS"][kernelName];
  //    std::vector<std::string> keys = kernelNode.keys();
  //    std::cout << "--------------------------------------" << std::endl;
  //    std::cout << "best parameters:" << std::endl;
  //    for (std::string key : keys) {
  //        if (key.compare(0, 7, "KERNEL_", 0, 7) == 0) {
  //            std::cout << "key: " << key << " value: " <<
  //            kernelNode[key].get() << std::endl;
  //        }
  //    }

  bestParameters.serialize("bestParameters.cfg");

  //    staticParameterTuner.writeStatisticsToFile("statistics.csv",
  //    platformName, deviceName, kernelName);

  std::cout << "-------------- all done! --------------" << std::endl;

  return 0;
}
