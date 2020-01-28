// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/datadriven/operation/hash/OperationMultipleEvalStreamingBSplineOCL/StreamingBSplineOCLOperatorFactory.hpp>

#include <sgpp/base/operation/hash/OperationMultipleEval.hpp>
#include <sgpp/base/exception/factory_exception.hpp>
#include <sgpp/globaldef.hpp>
#include <sgpp/base/opencl/OCLOperationConfiguration.hpp>

namespace sgpp {
namespace datadriven {

base::OperationMultipleEval* createStreamingBSplineOCLConfigured(
    base::Grid& grid, base::DataMatrix& dataset,
    sgpp::datadriven::OperationMultipleEvalConfiguration& configuration) {
  std::shared_ptr<base::OCLOperationConfiguration> parameters;

  if (configuration.getParameters().operator bool()) {
    base::OCLOperationConfiguration* cloned =
        dynamic_cast<base::OCLOperationConfiguration*>(configuration.getParameters()->clone());
    parameters = std::shared_ptr<base::OCLOperationConfiguration>(cloned);
  } else {
    parameters = std::make_shared<base::OCLOperationConfiguration>("StreamingBSplineOCL.cfg");

    if ((*parameters).contains("KERNEL_USE_LOCAL_MEMORY") == false) {
      (*parameters).addTextAttr("KERNEL_USE_LOCAL_MEMORY", "true");
    }

    if ((*parameters).contains("KERNEL_DATA_BLOCK_SIZE") == false) {
      (*parameters).addTextAttr("KERNEL_DATA_BLOCK_SIZE", "1");
    }

    if ((*parameters).contains("LINEAR_LOAD_BALANCING_VERBOSE") == false) {
      (*parameters).addTextAttr("LINEAR_LOAD_BALANCING_VERBOSE", "false");
    }

    if ((*parameters).contains("KERNEL_TRANS_DATA_BLOCK_SIZE") == false) {
      (*parameters).addTextAttr("KERNEL_TRANS_DATA_BLOCK_SIZE", "1");
    }
  }

  if ((*parameters)["VERBOSE"].getBool()) {
    std::cout << "are optimizations on: " << (*parameters)["ENABLE_OPTIMIZATIONS"].getBool()
              << std::endl;
    std::cout << "is local memory on: " << (*parameters)["KERNEL_USE_LOCAL_MEMORY"].getBool()
              << std::endl;
    std::cout << "local size: " << (*parameters)["LOCAL_SIZE"].getUInt() << std::endl;
    std::cout << "internal precision: " << (*parameters)["INTERNAL_PRECISION"].get() << std::endl;
    std::cout << "platform is: " << (*parameters)["PLATFORM"].get() << std::endl;
    std::cout << "device type is: " << (*parameters)["DEVICE_TYPE"].get() << std::endl;
  }

  if ((*parameters)["INTERNAL_PRECISION"].get() == "float") {
    return new datadriven::OperationMultiEvalStreamingBSplineOCL<float>(grid, dataset, parameters);
  } else if ((*parameters)["INTERNAL_PRECISION"].get() == "double") {
    return new datadriven::OperationMultiEvalStreamingBSplineOCL<double>(grid, dataset, parameters);
  } else {
    throw base::factory_exception(
        "Error creating operation\"OperationMultiEvalStreamingBSplineOCL\": invalid value for "
        "parameter \"INTERNAL_PRECISION\"");
  }
}
}  // namespace datadriven
}  // namespace sgpp
