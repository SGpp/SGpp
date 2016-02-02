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
#include "Configuration.hpp"

namespace SGPP {
namespace datadriven {

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
        manager = std::make_shared<base::OCLManagerMultiPlatform>();
        parameters = manager->getConfiguration();
    }

    //TODO: filter devices that are disabled (COUNT=0)

    StreamingOCLMultiPlatform::Configuration::augmentDefaultParameters(*parameters);

//    std::string &firstPlatformName = (*parameters)["PLATFORMS"].keys()[0];
//    std::string &firstDeviceName = (*parameters)["PLATFORMS"][firstPlatformName]["DEVICES"].keys()[0];
//    json::Node &deviceNode = (*parameters)["PLATFORMS"][firstPlatformName]["DEVICES"][firstDeviceName];
//    json::Node &firstDeviceConfig = deviceNode["KERNELS"][StreamingOCLMultiPlatform::Configuration::getKernelName()];

    if ((*parameters)["INTERNAL_PRECISION"].get().compare("float") == 0) {
        return new datadriven::StreamingOCLMultiPlatform::OperationMultiEvalStreamingOCLMultiPlatform<float>(grid,
                dataset, manager, parameters, (*parameters));
    } else if ((*parameters)["INTERNAL_PRECISION"].get().compare("double") == 0) {
        return new datadriven::StreamingOCLMultiPlatform::OperationMultiEvalStreamingOCLMultiPlatform<double>(grid,
                dataset, manager, parameters, (*parameters));
    } else {
        throw base::factory_exception(
                "Error creating operation\"OperationMultiEvalStreamingOCLMultiPlatform\": invalid value for parameter \"INTERNAL_PRECISION\"");
    }
}

}
}
