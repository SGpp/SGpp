/*
 * kernelTuner.cpp
 *
 *  Created on: Oct 27, 2015
 *      Author: pfandedd
 */

#include "sgpp/datadriven/application/MetaLearner.hpp"
#include "sgpp/datadriven/operation/hash/simple/DatadrivenOperationCommon.hpp"
#include <sgpp/datadriven/opencl/OCLConfigurationParameters.hpp>

int main(int argc, char **argv) {

//TODO: export parameter names and value ranges
//TODO: export constraints (or constraint checkers -> function pointers)

    // Specify scenario for performance optimization

    int maxLevel = 9;
    std::string fileName = "DR5_train.arff";
    sg::base::RegularGridConfiguration gridConfig;
    sg::solver::SLESolverConfiguration SLESolverConfigRefine;
    sg::solver::SLESolverConfiguration SLESolverConfigFinal;
    sg::base::AdpativityConfiguration adaptConfig;

    // setup grid
    gridConfig.dim_ = 0; //dim is inferred from the data
    gridConfig.level_ = maxLevel;
    gridConfig.type_ = sg::base::Linear;

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
    SLESolverConfigRefine.type_ = sg::solver::CG;

    // Set solver for final step
    SLESolverConfigFinal.eps_ = 0;
    SLESolverConfigFinal.maxIterations_ = 5;
    SLESolverConfigFinal.threshold_ = -1.0;
    SLESolverConfigFinal.type_ = sg::solver::CG;

    double lambda = 0.000001;

    bool verbose = true;
    SGPP::datadriven::MetaLearner learner(gridConfig, SLESolverConfigRefine, SLESolverConfigFinal, adaptConfig, lambda,
            verbose);

    SGPP::base::OCLConfigurationParameters parameters;

    parameters.set("OCL_MANAGER_VERBOSE", "false");
    parameters.set("VERBOSE", "false");
    parameters.set("PLATFORM", "first");
    parameters.set("SELECT_SPECIFIC_DEVICE", "DISABLED");
    parameters.set("MAX_DEVICES", "1");
    parameters.set("INTERNAL_PRECISION", "float");

    parameters.set("KERNEL_USE_LOCAL_MEMORY", "false");
    parameters.set("KERNEL_DATA_BLOCKING_SIZE", "4");
    parameters.set("KERNEL_TRANS_GRID_BLOCKING_SIZE", "4");
    parameters.set("KERNEL_STORE_DATA", "array");
    parameters.set("KERNEL_MAX_DIM_UNROLL", "1");

    SGPP::datadriven::OperationMultipleEvalConfiguration configuration(
    SGPP::datadriven::OperationMultipleEvalType::STREAMING,
    SGPP::datadriven::OperationMultipleEvalSubType::OCLMP, parameters);

    learner.learn(configuration, fileName);

    return 0;
}


