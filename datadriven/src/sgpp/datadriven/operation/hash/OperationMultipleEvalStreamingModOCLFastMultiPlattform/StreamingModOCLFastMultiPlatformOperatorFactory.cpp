// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <string>

#include "sgpp/base/operation/hash/OperationMultipleEval.hpp"
#include "sgpp/base/exception/factory_exception.hpp"
#include "sgpp/globaldef.hpp"
#include "StreamingModOCLFastMultiPlatformOperatorFactory.hpp"
#include "StreamingModOCLFastMultiPlatformConfiguration.hpp"
#include "sgpp/base/opencl/OCLOperationConfiguration.hpp"

namespace SGPP {
namespace datadriven {

base::OperationMultipleEval *createStreamingModOCLFastMultiPlatformConfigured(
    base::Grid &grid, base::DataMatrix &dataset,
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

  StreamingModOCLFastMultiPlatformConfiguration::augmentDefaultParameters(*parameters);

  std::string &firstPlatformName = (*parameters)["PLATFORMS"].keys()[0];
  std::string &firstDeviceName = (*parameters)["PLATFORMS"][firstPlatformName]["DEVICES"].keys()[0];
  json::Node &deviceNode =
      (*parameters)["PLATFORMS"][firstPlatformName]["DEVICES"][firstDeviceName];
  json::Node &firstDeviceConfig =
      deviceNode["KERNELS"][StreamingModOCLFastMultiPlatformConfiguration::getKernelName()];

  if (firstDeviceConfig["INTERNAL_PRECISION"].get() == "float") {
    return new datadriven::OperationMultiEvalStreamingModOCLFastMultiPlatform<float>(
        grid, dataset, manager, parameters, firstDeviceConfig);
  } else if (firstDeviceConfig["INTERNAL_PRECISION"].get() == "double") {
    return new datadriven::OperationMultiEvalStreamingModOCLFastMultiPlatform<double>(
        grid, dataset, manager, parameters, firstDeviceConfig);
  } else {
    throw base::factory_exception(
        "Error creating operation\"OperationMultiEvalStreamingModOCL\": invalid value for "
        "parameter \"INTERNAL_PRECISION\"");
  }
}
}  // namespace datadriven
}  // namespace SGPP
