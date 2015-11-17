/*
 * OCLOperatorFactory.hpp
 *
 *  Created on: Mar 25, 2015
 *      Author: pfandedd
 */

#include <sgpp/base/operation/hash/OperationMultipleEval.hpp>
#include <sgpp/base/exception/factory_exception.hpp>
#include <sgpp/globaldef.hpp>
#include "StreamingModOCLFastMultiPlatformOperatorFactory.hpp"
#include <sgpp/base/opencl/OCLOperationConfiguration.hpp>

namespace SGPP {
namespace datadriven {

base::OperationMultipleEval* createStreamingModOCLFastMultiPlatformConfigured(base::Grid& grid,
        base::DataMatrix& dataset,
        SGPP::datadriven::OperationMultipleEvalConfiguration &configuration) {

    std::shared_ptr<base::OCLOperationConfiguration> parameters;

    if (configuration.getParameters().operator bool()) {
        base::OCLOperationConfiguration *cloned = dynamic_cast<base::OCLOperationConfiguration *>(configuration.getParameters()->clone());
        parameters = std::shared_ptr<base::OCLOperationConfiguration>(cloned);
    } else {
        parameters = std::make_shared<base::OCLOperationConfiguration>("StreamingModOCLFastMultiPlatform.cfg");

        if ((*parameters).contains("KERNEL_USE_LOCAL_MEMORY") == false) {
            (*parameters).addIDAttr("KERNEL_USE_LOCAL_MEMORY", true);
        }

        if ((*parameters).contains("KERNEL_DATA_BLOCKING_SIZE") == false) {
            (*parameters).addIDAttr("KERNEL_DATA_BLOCKING_SIZE", 1ul);
        }

        if ((*parameters).contains("LINEAR_LOAD_BALANCING_VERBOSE") == false) {
            (*parameters).addIDAttr("LINEAR_LOAD_BALANCING_VERBOSE", false);
        }

        if ((*parameters).contains("KERNEL_TRANS_DATA_BLOCK_SIZE") == false) {
            (*parameters).addIDAttr("KERNEL_TRANS_DATA_BLOCK_SIZE", 1ul);
        }

        if ((*parameters).contains("KERNEL_TRANS_UNROLL_1D") == false) {
            (*parameters).addIDAttr("KERNEL_TRANS_UNROLL_1D", true);
        }

        if ((*parameters).contains("KERNEL_STORE_DATA") == false) {
            (*parameters).addTextAttr("KERNEL_STORE_DATA", "array");
        }
    }

    if ((*parameters)["VERBOSE"].getBool()) {
        std::cout << "are optimizations on: " << (*parameters)["ENABLE_OPTIMIZATIONS"].getBool() << std::endl;
        std::cout << "is local memory on: " << (*parameters)["KERNEL_USE_LOCAL_MEMORY"].getBool() << std::endl;
        std::cout << "local size: " << (*parameters)["LOCAL_SIZE"].getUInt() << std::endl;
        std::cout << "internal precision: " << (*parameters)["INTERNAL_PRECISION"].get() << std::endl;
        std::cout << "platform is: " << (*parameters)["PLATFORM"].get() << std::endl;
        std::cout << "device type is: " << (*parameters)["DEVICE_TYPE"].get() << std::endl;
        std::cout << "KERNEL_USE_LOCAL_MEMORY: " << (*parameters)["KERNEL_USE_LOCAL_MEMORY"].get() << std::endl;
        std::cout << "KERNEL_DATA_BLOCKING_SIZE: " << (*parameters)["KERNEL_DATA_BLOCKING_SIZE"].get() << std::endl;
        std::cout << "KERNEL_TRANS_DATA_BLOCK_SIZE: " << (*parameters)["KERNEL_TRANS_DATA_BLOCK_SIZE"].get()
                << std::endl;
        std::cout << "KERNEL_TRANS_GRID_BLOCK_SIZE: " << (*parameters)["KERNEL_TRANS_GRID_BLOCK_SIZE"].get()
                << std::endl;
        std::cout << "KERNEL_MAX_DIM_UNROLL: " << (*parameters)["KERNEL_MAX_DIM_UNROLL"].get() << std::endl;
        std::cout << "KERNEL_STORE_DATA: " << (*parameters)["KERNEL_STORE_DATA"].get() << std::endl;
    }

    if ((*parameters)["INTERNAL_PRECISION"].get() == "float") {
        return new datadriven::OperationMultiEvalStreamingModOCLFastMultiPlatform<float>(grid, dataset, parameters);
    } else if ((*parameters)["INTERNAL_PRECISION"].get() == "double") {
        return new datadriven::OperationMultiEvalStreamingModOCLFastMultiPlatform<double>(grid, dataset, parameters);
    } else {
        throw base::factory_exception(
                "Error creating operation\"OperationMultiEvalStreamingModOCL\": invalid value for parameter \"INTERNAL_PRECISION\"");
    }
}

}
}
