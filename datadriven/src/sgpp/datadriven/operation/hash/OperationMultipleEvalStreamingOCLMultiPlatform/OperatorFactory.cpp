// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/datadriven/operation/hash/OperationMultipleEvalStreamingOCLMultiPlatform/Configuration.hpp>
#include <sgpp/datadriven/operation/hash/OperationMultipleEvalStreamingOCLMultiPlatform/OperationMultiEvalStreamingOCLMultiPlatform.hpp>
#include <sgpp/base/exception/factory_exception.hpp>
#include <sgpp/base/opencl/OCLOperationConfiguration.hpp>
#include <sgpp/base/operation/hash/OperationMultipleEval.hpp>
#include <sgpp/datadriven/operation/hash/DatadrivenOperationCommon.hpp>
#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace datadriven {

base::OperationMultipleEval* createStreamingOCLMultiPlatformConfigured(
    base::Grid& grid, base::DataMatrix& dataset,
    sgpp::datadriven::OperationMultipleEvalConfiguration& configuration) {
  std::shared_ptr<base::OCLManagerMultiPlatform> manager;

  std::shared_ptr<base::OCLOperationConfiguration> parameters;
  if (configuration.getParameters().operator bool()) {
    base::OCLOperationConfiguration* cloned =
        dynamic_cast<base::OCLOperationConfiguration*>(configuration.getParameters()->clone());
    parameters = std::shared_ptr<base::OCLOperationConfiguration>(cloned);
    manager = std::make_shared<base::OCLManagerMultiPlatform>(parameters);
  } else {
    manager = std::make_shared<base::OCLManagerMultiPlatform>();
    parameters = manager->getConfiguration();
  }

  StreamingOCLMultiPlatform::Configuration::augmentDefaultParameters(*parameters);

  //    std::string &firstPlatformName = (*parameters)["PLATFORMS"].keys()[0];
  //    std::string &firstDeviceName =
  //    (*parameters)["PLATFORMS"][firstPlatformName]["DEVICES"].keys()[0];
  //    json::Node &deviceNode =
  //    (*parameters)["PLATFORMS"][firstPlatformName]["DEVICES"][firstDeviceName];
  //    json::Node &firstDeviceConfig =
  //    deviceNode["[StreamingOCLMultiPlatform::Configuration::getKernelName()];
  //  std::cout << "INTERNAL_PRECISION: " << (*parameters)["INTERNAL_PRECISION"].get() << std::endl;
  if ((*parameters)["INTERNAL_PRECISION"].get().compare("float") == 0) {
    return new datadriven::StreamingOCLMultiPlatform::OperationMultiEvalStreamingOCLMultiPlatform<
        float>(grid, dataset, manager, parameters);
  } else if ((*parameters)["INTERNAL_PRECISION"].get().compare("double") == 0) {
    return new datadriven::StreamingOCLMultiPlatform::OperationMultiEvalStreamingOCLMultiPlatform<
        double>(grid, dataset, manager, parameters);
  } else {
    throw base::factory_exception(
        "Error creating "
        "operation\"OperationMultiEvalStreamingOCLMultiPlatform\": invalid "
        "value for parameter \"INTERNAL_PRECISION\"");
  }
}

}  // namespace datadriven
}  // namespace sgpp
