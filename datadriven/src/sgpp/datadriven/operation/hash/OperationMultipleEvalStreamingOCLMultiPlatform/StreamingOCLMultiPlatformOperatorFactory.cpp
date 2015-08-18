/*
 * OCLOperatorFactory.hpp
 *
 *  Created on: Mar 25, 2015
 *      Author: pfandedd
 */

#include <sgpp/base/operation/hash/OperationMultipleEval.hpp>
#include <sgpp/base/exception/factory_exception.hpp>
#include "OperationMultiEvalStreamingOCLMultiPlatform.hpp"

#include <sgpp/globaldef.hpp>
#include "../../../opencl/OCLConfigurationParameters.hpp"
#include <sgpp/datadriven/operation/hash/simple/DatadrivenOperationCommon.hpp>

namespace SGPP {
namespace datadriven {

base::OperationMultipleEval* createStreamingOCLMultiPlatformConfigured(base::Grid& grid, base::DataMatrix& dataset,
SGPP::datadriven::OperationMultipleEvalConfiguration &configuration) {

    std::shared_ptr<base::OCLConfigurationParameters> parameters = std::static_pointer_cast<
            base::OCLConfigurationParameters>(configuration.getParameters()->clone());

    if (parameters->doUseDefaults()) {
        parameters->set("KERNEL_USE_LOCAL_MEMORY", "false");
        parameters->set("KERNEL_STORE_DATA", "array");
        parameters->set("KERNEL_MAX_DIM_UNROLL", "10");
        parameters->set("LINEAR_LOAD_BALANCING_VERBOSE", "false");
        parameters->set("KERNEL_DATA_BLOCKING_SIZE", "1");
        parameters->set("KERNEL_TRANS_GRID_BLOCKING_SIZE", "1");
        parameters->set("LINEAR_LOAD_BALANCING_VERBOSE", "false");

        parameters->readFromFile("StreamingOCLMultiPlatform.cfg");
    }

    if (parameters->getAsBoolean("VERBOSE")) {
        std::cout << "are optimizations on: " << parameters->getAsBoolean("ENABLE_OPTIMIZATIONS") << std::endl;
        std::cout << "is local memory on: " << parameters->getAsBoolean("KERNEL_USE_LOCAL_MEMORY") << std::endl;
        std::cout << "local size: " << parameters->getAsUnsigned("LOCAL_SIZE") << std::endl;
        std::cout << "internal precision: " << parameters->get("INTERNAL_PRECISION") << std::endl;
        std::cout << "platform is: " << parameters->get("PLATFORM") << std::endl;
        std::cout << "device type is: " << parameters->get("DEVICE_TYPE") << std::endl;
        std::cout << "KERNEL_USE_LOCAL_MEMORY: " << parameters->get("KERNEL_USE_LOCAL_MEMORY") << std::endl;
    }

    if (parameters->get("INTERNAL_PRECISION") == "float") {
        return new datadriven::OperationMultiEvalStreamingOCLMultiPlatform<float>(grid, dataset, parameters);
    } else if (parameters->get("INTERNAL_PRECISION") == "double") {
        return new datadriven::OperationMultiEvalStreamingOCLMultiPlatform<double>(grid, dataset, parameters);
    } else {
        throw base::factory_exception(
                "Error creating operation\"OperationMultiEvalStreamingOCLMultiPlatform\": invalid value for parameter \"INTERNAL_PRECISION\"");
    }
}

}
}
