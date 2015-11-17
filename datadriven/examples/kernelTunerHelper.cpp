/*
 * kernelTuner.cpp
 *
 *  Created on: Oct 27, 2015
 *      Author: pfandedd
 */
#include <iostream>

#if USE_OCL == 1
#include "sgpp/datadriven/application/MetaLearner.hpp"
#include "sgpp/datadriven/operation/hash/simple/DatadrivenOperationCommon.hpp"
#include <sgpp/base/opencl/OCLOperationConfiguration.hpp>
#include <sgpp/datadriven/application/StaticParameterTuner.hpp>
#include <sgpp/datadriven/application/TunableParameter.hpp>
#include <sgpp/datadriven/application/LearnerScenario.hpp>

int main(int argc, char **argv) {

//TODO: export parameter names and value ranges
//TODO: export constraints (or constraint checkers -> function pointers)

    // Specify scenario for performance optimization

    int maxLevel = 1;
//    std::string fileName = "friedman_4d_small.arff";
    std::string fileName = "friedman_4d.arff";
    sg::base::RegularGridConfiguration gridConfig;
    sg::solver::SLESolverConfiguration SLESolverConfigRefine;
    sg::solver::SLESolverConfiguration SLESolverConfigFinal;
    sg::base::AdpativityConfiguration adaptConfig;

    // setup grid
    gridConfig.dim_ = 0; //dim is inferred from the data
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
    SLESolverConfigRefine.maxIterations_ = 5;
    SLESolverConfigRefine.threshold_ = -1.0;
    SLESolverConfigRefine.type_ = SGPP::solver::SLESolverType::CG;

    // Set solver for final step
    SLESolverConfigFinal.eps_ = 0;
    SLESolverConfigFinal.maxIterations_ = 20;
    SLESolverConfigFinal.threshold_ = -1.0;
    SLESolverConfigFinal.type_ = SGPP::solver::SLESolverType::CG;

    double lambda = 0.000001;

    SGPP::datadriven::LearnerScenario scenario(fileName, lambda, gridConfig, SLESolverConfigRefine, SLESolverConfigFinal, adaptConfig);

    SGPP::datadriven::StaticParameterTuner staticParameterTuner(true, true);

    staticParameterTuner.addFixedParameter("OCL_MANAGER_VERBOSE", "false", SGPP::datadriven::ParameterType::BOOL);
    staticParameterTuner.addFixedParameter("VERBOSE", "false", SGPP::datadriven::ParameterType::BOOL);
    staticParameterTuner.addFixedParameter("PLATFORM", "NVIDIA CUDA", SGPP::datadriven::ParameterType::TEXT);
    staticParameterTuner.addFixedParameter("SELECT_SPECIFIC_DEVICE", "DISABLED", SGPP::datadriven::ParameterType::TEXT);
    staticParameterTuner.addFixedParameter("MAX_DEVICES", "1", SGPP::datadriven::ParameterType::TEXT);
    staticParameterTuner.addFixedParameter("INTERNAL_PRECISION", "float", SGPP::datadriven::ParameterType::DOUBLE);

    staticParameterTuner.addParameter("KERNEL_USE_LOCAL_MEMORY", {"false", "true"}, SGPP::datadriven::ParameterType::BOOL);
    staticParameterTuner.addParameter("KERNEL_DATA_BLOCKING_SIZE", {"1", "2", "4", "8"}, SGPP::datadriven::ParameterType::UINT);
    staticParameterTuner.addParameter("KERNEL_TRANS_GRID_BLOCKING_SIZE", {"1", "2", "4", "8"}, SGPP::datadriven::ParameterType::UINT);
    staticParameterTuner.addParameter("KERNEL_STORE_DATA", {"register"}, SGPP::datadriven::ParameterType::TEXT); //"array",
    staticParameterTuner.addParameter("KERNEL_MAX_DIM_UNROLL", {"4"}, SGPP::datadriven::ParameterType::UINT); //"1", "8", "16"

//    staticParameterTuner.writeToFile("testOut.tuner");
//
//    SGPP::datadriven::StaticParameterTuner readStaticParameterTuner("testOut.tuner");

    SGPP::base::OCLOperationConfiguration bestParameters = staticParameterTuner.tuneParameters(scenario);

    std::vector<std::string> keys = bestParameters.keys();
    std::cout << "--------------------------------------" << std::endl;
    std::cout << "best parameters:" << std::endl;
    for (std::string key : keys) {
        if (key.compare(0, 7, "KERNEL_", 0, 7) == 0) {
            std::cout << "key: " << key << " value: " << bestParameters[key].get() << std::endl;
        }
    }

    bestParameters.serialize("bestParameters.cfg");

    staticParameterTuner.writeStatisticsToFile("statistics.csv");

    return 0;
}
#else
int main(int argc, char **argv) {
  std::cout << "no OpenCL support" << std::endl;
}
#endif

