/*
 * OCLOperatorFactory.hpp
 *
 *  Created on: Mar 25, 2015
 *      Author: pfandedd
 */

#include <sgpp/base/operation/hash/OperationMultipleEval.hpp>
#include <sgpp/base/exception/factory_exception.hpp>
#include <sgpp/base/opencl/OCLConfigurationParameters.hpp>
#include "OperationMultiEvalStreamingOCLMultiPlatform.hpp"

#include <sgpp/globaldef.hpp>

namespace SGPP {
namespace datadriven {

base::OperationMultipleEval* createStreamingOCLMultiPlatformConfigured(base::Grid& grid, base::DataMatrix& dataset,
        base::OCLConfigurationParameters *parameters) {

    if (parameters == nullptr) {
        std::map<std::string, std::string> defaultParameter;
        defaultParameter["KERNEL_USE_LOCAL_MEMORY"] = "false";
        defaultParameter["KERNEL_STORE_DATA"] = "array";
        defaultParameter["KERNEL_MAX_DIM_UNROLL"] = "10";
        defaultParameter["LINEAR_LOAD_BALANCING_VERBOSE"] = "false";
        defaultParameter["KERNEL_DATA_BLOCKING_SIZE"] = "1";
        defaultParameter["KERNEL_TRANS_GRID_BLOCKING_SIZE"] = "1";
        defaultParameter["LINEAR_LOAD_BALANCING_VERBOSE"] = "false";

        parameters = new base::OCLConfigurationParameters("StreamingOCLMultiPlatform.cfg", defaultParameter);
    }

    if (parameters->getAsBoolean("KERNEL_VERBOSE")) {
        std::cout << "are optimizations on: " << parameters->getAsBoolean("ENABLE_OPTIMIZATIONS") << std::endl;
        std::cout << "is local memory on: " << parameters->getAsBoolean("KERNEL_USE_LOCAL_MEMORY") << std::endl;
        std::cout << "local size: " << parameters->getAsUnsigned("LOCAL_SIZE") << std::endl;
        std::cout << "internal precision: " << (*parameters)["INTERNAL_PRECISION"] << std::endl;
        std::cout << "platform is: " << (*parameters)["PLATFORM"] << std::endl;
        std::cout << "device type is: " << (*parameters)["DEVICE_TYPE"] << std::endl;
        std::cout << "KERNEL_USE_LOCAL_MEMORY: " << (*parameters)["KERNEL_USE_LOCAL_MEMORY"] << std::endl;
    }

    if ((*parameters)["INTERNAL_PRECISION"] == "float") {
        return new datadriven::OperationMultiEvalStreamingOCLMultiPlatform<float>(grid, dataset, *parameters);
    } else if ((*parameters)["INTERNAL_PRECISION"] == "double") {
        return new datadriven::OperationMultiEvalStreamingOCLMultiPlatform<double>(grid, dataset, *parameters);
    } else {
        throw base::factory_exception(
                "Error creating operation\"OperationMultiEvalStreamingOCLMultiPlatform\": invalid value for parameter \"INTERNAL_PRECISION\"");
    }
}

}
}
