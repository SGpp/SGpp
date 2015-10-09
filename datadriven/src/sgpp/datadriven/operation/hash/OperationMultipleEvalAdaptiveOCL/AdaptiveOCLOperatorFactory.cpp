/*
 * OCLOperatorFactory.hpp
 *
 *  Created on: Mar 25, 2015
 *      Author: pfandedd
 */

#include <sgpp/base/operation/hash/OperationMultipleEval.hpp>
#include <sgpp/base/exception/factory_exception.hpp>
#include <sgpp/globaldef.hpp>
#include "../../../opencl/OCLConfigurationParameters.hpp"
#include "OperationMultiEvalAdaptiveOCL.hpp"

namespace SGPP {
namespace datadriven {

base::OperationMultipleEval *createAdaptiveOCLConfigured(base::Grid& grid, base::DataMatrix& dataset,
    SGPP::datadriven::OperationMultipleEvalConfiguration &configuration) {


    std::shared_ptr<base::OCLConfigurationParameters> parameters(new base::OCLConfigurationParameters());

    if (configuration.getParameters().operator bool()) {
        parameters = std::static_pointer_cast<base::OCLConfigurationParameters>(configuration.getParameters()->clone());
    } else {
        parameters->set("KERNEL_USE_LOCAL_MEMORY", "false");
        parameters->set("LOCAL_SIZE", "128");
        parameters->set("KERNEL_DATA_BLOCKING_SIZE", "1");
        parameters->set("LINEAR_LOAD_BALANCING_VERBOSE", "false");
        parameters->set("KERNEL_TRANS_DATA_BLOCK_SIZE", "1");
        parameters->set("ADAPTIVE_STREAMING_HARD_LIMIT", "10"); //absolute value
        parameters->set("ADAPTIVE_STREAMING_DENSITY", "5"); //In percent

        parameters->readFromFile("AdaptiveOCL.cfg");
    }

    if (parameters->getAsBoolean("VERBOSE")) {
        std::cout << "are optimizations on: " << parameters->getAsBoolean("ENABLE_OPTIMIZATIONS") << std::endl;
        std::cout << "is local memory on: " << parameters->getAsBoolean("KERNEL_USE_LOCAL_MEMORY") << std::endl;
        std::cout << "local size: " << parameters->getAsUnsigned("LOCAL_SIZE") << std::endl;
        std::cout << "internal precision: " << parameters->get("INTERNAL_PRECISION") << std::endl;
        std::cout << "platform is: " << parameters->get("PLATFORM") << std::endl;
        std::cout << "device type is: " << parameters->get("DEVICE_TYPE") << std::endl;
        std::cout << "hard limit: " << parameters->get("ADAPTIVE_STREAMING_HARD_LIMIT") << std::endl;
        std::cout << "soft limit: " << parameters->get("ADAPTIVE_STREAMING_DENSITY") << std::endl;
    }

    if (parameters->get("INTERNAL_PRECISION") == "float") {
        return new datadriven::OperationMultiEvalAdaptiveOCL<float>(grid, dataset, parameters);
    } else if (parameters->get("INTERNAL_PRECISION") == "double") {
        return new datadriven::OperationMultiEvalAdaptiveOCL<double>(grid, dataset, parameters);
    } else {
        throw base::factory_exception(
                "Error creating operation\"OperationMultiEvalAdaptiveOCL\": invalid value for parameter \"INTERNAL_PRECISION\"");
    }

    //TODO: make parameters changable through api
    /*std::shared_ptr<base::OCLConfigurationParameters> parameters;
    parameters->set("KERNEL_USE_LOCAL_MEMORY", "true");
    parameters->set("KERNEL_MAX_DIM_UNROLL", "10");
    parameters->set("LINEAR_LOAD_BALANCING_VERBOSE", "false");

    parameters->readFromFile("AdaptiveOCL.cfg");

//  std::cout << "are optimizations on: " << parameters->getAsBoolean("STREAMING_OCL_ENABLE_OPTIMIZATIONS") << std::endl;
//  std::cout << "is local memory on: " << parameters->getAsBoolean("STREAMING_OCL_USE_LOCAL_MEMORY") << std::endl;
//  std::cout << "local size: " << parameters->getAsUnsigned("STREAMING_OCL_LOCAL_SIZE") << std::endl;
//  std::cout << "max dim unroll: " << parameters->getAsUnsigned("STREAMING_OCL_MAX_DIM_UNROLL") << std::endl;
//  std::cout << "internal precision: " << parameters->get("STREAMING_OCL_INTERNAL_PRECISION") << std::endl;
//  std::cout << "platform is: " << parameters->get("STREAMING_OCL_PLATFORM") << std::endl;
//  std::cout << "device type is: " << parameters->get("STREAMING_OCL_DEVICE_TYPE") << std::endl;

    if (parameters->get("INTERNAL_PRECISION") == "float") {
        return new datadriven::OperationMultiEvalAdaptiveOCL<float>(grid, dataset, parameters);
    } else if (parameters->get("INTERNAL_PRECISION") == "double") {
        return new datadriven::OperationMultiEvalAdaptiveOCL<double>(grid, dataset, parameters);
    } else {
        throw base::factory_exception(
                "Error creating operation\"OperationMultiEvalAdaptiveOCL\": invalid value for parameter \"INTERNAL_PRECISION\"");
    }*/
}

}
}
