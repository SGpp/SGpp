/*
 * OCLOperatorFactory.hpp
 *
 *  Created on: Mar 25, 2015
 *      Author: pfandedd
 */

#include "StreamingBSplineOCLOperatorFactory.hpp"

#include <sgpp/base/operation/hash/OperationMultipleEval.hpp>
#include <sgpp/base/exception/factory_exception.hpp>
#include <sgpp/globaldef.hpp>
#include <sgpp/base/opencl/OCLConfigurationParameters.hpp>

namespace SGPP {
namespace datadriven {

base::OperationMultipleEval* createStreamingBSplineOCLConfigured(base::Grid& grid, base::DataMatrix& dataset,
SGPP::datadriven::OperationMultipleEvalConfiguration &configuration) {

    std::shared_ptr<base::OCLConfigurationParameters> parameters;

    if (configuration.getParameters().operator bool()) {
        parameters = std::dynamic_pointer_cast<base::OCLConfigurationParameters>(configuration.getParameters()->clone());
    } else {
        parameters = std::make_shared<base::OCLConfigurationParameters>();
        parameters->set("KERNEL_USE_LOCAL_MEMORY", "true");
        parameters->set("KERNEL_DATA_BLOCKING_SIZE", "1");
        parameters->set("LINEAR_LOAD_BALANCING_VERBOSE", "false");
        parameters->set("KERNEL_TRANS_DATA_BLOCK_SIZE", "1");

        parameters->readFromFile("StreamingBSplineOCL.cfg");
    }

    if (parameters->getAsBoolean("VERBOSE")) {
        std::cout << "are optimizations on: " << parameters->getAsBoolean("ENABLE_OPTIMIZATIONS") << std::endl;
        std::cout << "is local memory on: " << parameters->getAsBoolean("KERNEL_USE_LOCAL_MEMORY") << std::endl;
        std::cout << "local size: " << parameters->getAsUnsigned("LOCAL_SIZE") << std::endl;
        std::cout << "internal precision: " << parameters->get("INTERNAL_PRECISION") << std::endl;
        std::cout << "platform is: " << parameters->get("PLATFORM") << std::endl;
        std::cout << "device type is: " << parameters->get("DEVICE_TYPE") << std::endl;
    }

    if (parameters->get("INTERNAL_PRECISION") == "float") {
        return new datadriven::OperationMultiEvalStreamingBSplineOCL<float>(grid, dataset, parameters);
    } else if (parameters->get("INTERNAL_PRECISION") == "double") {
        return new datadriven::OperationMultiEvalStreamingBSplineOCL<double>(grid, dataset, parameters);
    } else {
        throw base::factory_exception(
                "Error creating operation\"OperationMultiEvalStreamingBSplineOCL\": invalid value for parameter \"INTERNAL_PRECISION\"");
    }
}

}
}
