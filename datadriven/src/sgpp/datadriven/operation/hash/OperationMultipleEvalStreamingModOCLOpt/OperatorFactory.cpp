// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/operation/hash/OperationMultipleEval.hpp>
#include <sgpp/base/exception/factory_exception.hpp>
#include <sgpp/globaldef.hpp>
#include <sgpp/datadriven/operation/hash/DatadrivenOperationCommon.hpp>
#include <sgpp/datadriven/operation/hash/OperationMultipleEvalStreamingModOCLOpt/Configuration.hpp>
#include <sgpp/datadriven/operation/hash/OperationMultipleEvalStreamingModOCLOpt/OperationMultiEvalStreamingModOCLOpt.hpp>

namespace sgpp {
namespace datadriven {

base::OperationMultipleEval* createStreamingModOCLOptConfigured(
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

  StreamingModOCLOpt::Configuration::augmentDefaultParameters(*parameters);

  if ((*parameters)["INTERNAL_PRECISION"].get() == "float") {
    return new datadriven::OperationMultiEvalStreamingModOCLOpt<float>(grid, dataset, manager,
                                                                       parameters);
  } else if ((*parameters)["INTERNAL_PRECISION"].get() == "double") {
    return new datadriven::OperationMultiEvalStreamingModOCLOpt<double>(grid, dataset, manager,
                                                                        parameters);
  } else {
    throw base::factory_exception(
        "Error creating "
        "operation\"OperationMultiEvalStreamingModOCLOpt\": "
        "invalid value for parameter \"INTERNAL_PRECISION\"");
  }
}

}  // namespace datadriven
}  // namespace sgpp
