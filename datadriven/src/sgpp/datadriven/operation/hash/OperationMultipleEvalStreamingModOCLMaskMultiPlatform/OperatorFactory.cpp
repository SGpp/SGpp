/*
 * OCLOperatorFactory.hpp
 *
 *  Created on: Mar 25, 2015
 *      Author: pfandedd
 */

#include <sgpp/base/operation/hash/OperationMultipleEval.hpp>
#include <sgpp/base/exception/factory_exception.hpp>
#include <sgpp/globaldef.hpp>
#include <sgpp/datadriven/operation/hash/simple/DatadrivenOperationCommon.hpp>
#include "Configuration.hpp"
#include "OperationMultiEvalStreamingModOCLMaskMultiPlatform.hpp"

namespace SGPP {
namespace datadriven {

base::OperationMultipleEval* createStreamingModOCLMaskMultiPlatformConfigured(base::Grid& grid, base::DataMatrix& dataset,
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

    StreamingModOCLMaskMultiPlatform::Configuration::augmentDefaultParameters(*parameters);

    if ((*parameters)["INTERNAL_PRECISION"].get() == "float") {
        return new datadriven::OperationMultiEvalStreamingModOCLMaskMultiPlatform<float>(grid, dataset, manager, parameters);
    } else if ((*parameters)["INTERNAL_PRECISION"].get() == "double") {
        return new datadriven::OperationMultiEvalStreamingModOCLMaskMultiPlatform<double>(grid, dataset, manager, parameters);
    } else {
        throw base::factory_exception(
                "Error creating operation\"OperationMultiEvalStreamingModOCLMaskMultiPlatform\": invalid value for parameter \"INTERNAL_PRECISION\"");
    }

}

}
}
