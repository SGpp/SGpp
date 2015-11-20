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
#include <sgpp/datadriven/operation/hash/simple/DatadrivenOperationCommon.hpp>
#include <sgpp/base/opencl/OCLOperationConfiguration.hpp>
#include <sgpp/datadriven/operation/hash/OperationMultipleEvalStreamingOCLMultiPlatform/StreamingOCLMultiPlatformConfiguration.hpp>

namespace SGPP {
namespace datadriven {

void augmentDefaultParameters(SGPP::base::OCLOperationConfiguration &parameters) {

    for (std::string &platformName : parameters["PLATFORMS"].keys()) {
        json::Node &platformNode = parameters["PLATFORMS"][platformName];
        for (std::string &deviceName : platformNode["DEVICES"].keys()) {
            json::Node &deviceNode = platformNode["DEVICES"][deviceName];

            const std::string &kernelName = SGPP::datadriven::StreamingOCLMultiPlatformConfiguration::getKernelName();

            json::Node &kernelNode =
                    deviceNode["KERNELS"].contains(kernelName) ?
                            deviceNode["KERNELS"][kernelName] : deviceNode["KERNELS"].addDictAttr(kernelName);
//            std::cout << "in kernel augment" << std::endl;
//            std::cout << "-----------------------------------" << std::endl;
//            for (std::string &key: kernelNode.keys()) {
//                std::cout << "key: " << key << " value: " << kernelNode[key].get() << std::endl;
//            }

            if (kernelNode.contains("LOCAL_SIZE") == false) {
                kernelNode.addIDAttr("LOCAL_SIZE", 128ul);
            }

            if (kernelNode.contains("KERNEL_USE_LOCAL_MEMORY") == false) {
                kernelNode.addIDAttr("KERNEL_USE_LOCAL_MEMORY", false);
            }

            if (kernelNode.contains("KERNEL_STORE_DATA") == false) {
                kernelNode.addTextAttr("KERNEL_STORE_DATA", "array");
            }

            if (kernelNode.contains("KERNEL_MAX_DIM_UNROLL") == false) {
                kernelNode.addIDAttr("KERNEL_MAX_DIM_UNROLL", 10ul);
            }

//            if (kernelNode.contains("LINEAR_LOAD_BALANCING_VERBOSE") == false) {
//                kernelNode.addIDAttr("LINEAR_LOAD_BALANCING_VERBOSE", false);
//            }

            if (kernelNode.contains("KERNEL_DATA_BLOCKING_SIZE") == false) {
                kernelNode.addIDAttr("KERNEL_DATA_BLOCKING_SIZE", 1ul);
            }

            if (kernelNode.contains("KERNEL_TRANS_GRID_BLOCKING_SIZE") == false) {
                kernelNode.addIDAttr("KERNEL_TRANS_GRID_BLOCKING_SIZE", 1ul);
            }



//            if (kernelNode.contains("LINEAR_LOAD_BALANCING_VERBOSE") == false) {
//                kernelNode.addIDAttr("LINEAR_LOAD_BALANCING_VERBOSE", false);
//            }
        }
    }
}

base::OperationMultipleEval* createStreamingOCLMultiPlatformConfigured(base::Grid& grid, base::DataMatrix& dataset,
SGPP::datadriven::OperationMultipleEvalConfiguration &configuration) {

    std::shared_ptr<base::OCLManagerMultiPlatform> manager;

    std::shared_ptr<base::OCLOperationConfiguration> parameters;
    if (configuration.getParameters().operator bool()) {
        base::OCLOperationConfiguration *cloned =
                dynamic_cast<base::OCLOperationConfiguration *>(configuration.getParameters()->clone());
        parameters = std::shared_ptr<base::OCLOperationConfiguration>(cloned);
        manager = std::make_shared<base::OCLManagerMultiPlatform>(parameters);
    } else {
//        parameters = std::make_shared<base::OCLOperationConfiguration>("StreamingOCLMultiPlatform.cfg");
        manager = std::make_shared<base::OCLManagerMultiPlatform>();
        parameters = manager->getConfiguration();
    }

    augmentDefaultParameters(*parameters);

    std::string &firstPlatformName = (*parameters)["PLATFORMS"].keys()[0];
    std::string &firstDeviceName = (*parameters)["PLATFORMS"][firstPlatformName]["DEVICES"].keys()[0];
    json::Node &deviceNode = (*parameters)["PLATFORMS"][firstPlatformName]["DEVICES"][firstDeviceName];
    json::Node &firstDeviceConfig = deviceNode["KERNELS"][StreamingOCLMultiPlatformConfiguration::getKernelName()];

    if (firstDeviceConfig["INTERNAL_PRECISION"].get().compare("float") == 0) {
        return new datadriven::OperationMultiEvalStreamingOCLMultiPlatform<float>(grid, dataset, manager, parameters, firstDeviceConfig);
    } else if (firstDeviceConfig["INTERNAL_PRECISION"].get().compare("double") == 0) {
        return new datadriven::OperationMultiEvalStreamingOCLMultiPlatform<double>(grid, dataset, manager, parameters, firstDeviceConfig);
    } else {
        throw base::factory_exception(
                "Error creating operation\"OperationMultiEvalStreamingOCLMultiPlatform\": invalid value for parameter \"INTERNAL_PRECISION\"");
    }
}

}
}

//    if ((*parameters)["VERBOSE"].getBool()) {
//        std::cout << "are optimizations on: " << (*parameters)["ENABLE_OPTIMIZATIONS"].getBool() << std::endl;
//        std::cout << "is local memory on: " << (*parameters)["KERNEL_USE_LOCAL_MEMORY"].getBool() << std::endl;
//        std::cout << "local size: " << (*parameters)["LOCAL_SIZE"].getUInt() << std::endl;
//        std::cout << "internal precision: " << (*parameters)["INTERNAL_PRECISION"].get() << std::endl;
//        std::cout << "platform is: " << (*parameters)["PLATFORM"].get() << std::endl;
//        std::cout << "device type is: " << (*parameters)["DEVICE_TYPE"].get() << std::endl;
//        std::cout << "KERNEL_USE_LOCAL_MEMORY: " << (*parameters)["KERNEL_USE_LOCAL_MEMORY"].get() << std::endl;
//        std::cout << "KERNEL_DATA_BLOCKING_SIZE: " << (*parameters)["KERNEL_DATA_BLOCKING_SIZE"].get() << std::endl;
//        std::cout << "KERNEL_TRANS_GRID_BLOCKING_SIZE: " << (*parameters)["KERNEL_TRANS_GRID_BLOCKING_SIZE"].get()
//                << std::endl;
//        std::cout << "LOAD_BALANCING_VERBOSE: " << (*parameters)["LOAD_BALANCING_VERBOSE"].get() << std::endl;
//        std::cout << "KERNEL_STORE_DATA: " << (*parameters)["KERNEL_STORE_DATA"].get() << std::endl;
//    }
