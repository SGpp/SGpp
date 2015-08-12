/*
 * OCLOperatorFactory.hpp
 *
 *  Created on: Mar 25, 2015
 *      Author: pfandedd
 */

#include <sgpp/base/operation/hash/OperationMultipleEval.hpp>
#include <sgpp/base/exception/factory_exception.hpp>
#include <sgpp/base/opencl/OCLConfigurationParameters.hpp>
#include <sgpp/globaldef.hpp>
#include "OperationMultiEvalAdaptiveOCL.hpp"

namespace SGPP {
namespace datadriven {

base::OperationMultipleEval *createAdaptiveOCLConfigured(base::Grid& grid, base::DataMatrix& dataset) {

    //TODO: make parameters changable through api
    base::OCLConfigurationParameters parameters;
    parameters.set("KERNEL_USE_LOCAL_MEMORY", "true");
    parameters.set("KERNEL_MAX_DIM_UNROLL", "10");
    parameters.set("LINEAR_LOAD_BALANCING_VERBOSE", "false");

    parameters.readFromFile("AdaptiveOCL.cfg");

//  std::cout << "are optimizations on: " << parameters.getAsBoolean("STREAMING_OCL_ENABLE_OPTIMIZATIONS") << std::endl;
//  std::cout << "is local memory on: " << parameters.getAsBoolean("STREAMING_OCL_USE_LOCAL_MEMORY") << std::endl;
//  std::cout << "local size: " << parameters.getAsUnsigned("STREAMING_OCL_LOCAL_SIZE") << std::endl;
//  std::cout << "max dim unroll: " << parameters.getAsUnsigned("STREAMING_OCL_MAX_DIM_UNROLL") << std::endl;
//  std::cout << "internal precision: " << parameters.get("STREAMING_OCL_INTERNAL_PRECISION") << std::endl;
//  std::cout << "platform is: " << parameters.get("STREAMING_OCL_PLATFORM") << std::endl;
//  std::cout << "device type is: " << parameters.get("STREAMING_OCL_DEVICE_TYPE") << std::endl;

    if (parameters.get("INTERNAL_PRECISION") == "float") {
        return new datadriven::OperationMultiEvalAdaptiveOCL<float>(grid, dataset, parameters);
    } else if (parameters.get("INTERNAL_PRECISION") == "double") {
        return new datadriven::OperationMultiEvalAdaptiveOCL<double>(grid, dataset, parameters);
    } else {
        throw base::factory_exception(
                "Error creating operation\"OperationMultiEvalAdaptiveOCL\": invalid value for parameter \"INTERNAL_PRECISION\"");
    }
}

}
}
