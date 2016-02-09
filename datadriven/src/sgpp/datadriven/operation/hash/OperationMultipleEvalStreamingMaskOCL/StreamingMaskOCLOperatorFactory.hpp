// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/base/operation/hash/OperationMultipleEval.hpp>
#include <sgpp/base/exception/factory_exception.hpp>
#include <sgpp/base/opencl/OpenCLConfigurationParameters.hpp>
#include <sgpp/datadriven/operation/hash/OperationMultipleEvalStreamingMaskOCL/OperationMultiEvalStreamingMaskOCL.hpp>

#include <sgpp/globaldef.hpp>

namespace SGPP {
namespace datadriven {

base::OperationMultipleEval* createStreamingMaskOCLConfigured(base::Grid& grid,
    base::DataMatrix& dataset) {

  std::map<std::string, std::string> defaultParameter;
  defaultParameter["STREAMING_MASK_OCL_USE_LOCAL_MEMORY"] = "true";
  defaultParameter["STREAMING_MASK_OCL_MAX_DIM_UNROLL"] = "10";
  defaultParameter["LINEAR_LOAD_BALANCING_VERBOSE"] = "false";

  base::OpenCLConfigurationParameters parameters("StreamingMaskOCL.cfg",
      defaultParameter);

  //  std::cout << "are optimizations on: " << parameters.getAsBoolean("STREAMING_OCL_ENABLE_OPTIMIZATIONS") << std::endl;
  //  std::cout << "is local memory on: " << parameters.getAsBoolean("STREAMING_OCL_USE_LOCAL_MEMORY") << std::endl;
  //  std::cout << "local size: " << parameters.getAsUnsigned("STREAMING_OCL_LOCAL_SIZE") << std::endl;
  //  std::cout << "max dim unroll: " << parameters.getAsUnsigned("STREAMING_OCL_MAX_DIM_UNROLL") << std::endl;
  //  std::cout << "internal precision: " << parameters["STREAMING_OCL_INTERNAL_PRECISION"] << std::endl;
  //  std::cout << "platform is: " << parameters["STREAMING_OCL_PLATFORM"] << std::endl;
  //  std::cout << "device type is: " << parameters["STREAMING_OCL_DEVICE_TYPE"] << std::endl;

  if (parameters["INTERNAL_PRECISION"] == "float") {
    return new datadriven::OperationMultiEvalStreamingMaskOCL<float>(grid, dataset,
           parameters);
  } else if (parameters["INTERNAL_PRECISION"] == "double") {
    return new datadriven::OperationMultiEvalStreamingMaskOCL<double>(grid, dataset,
           parameters);
  } else {
    throw base::factory_exception(
      "Error creating operation\"OperationMultiEvalStreamingMaskOCL\": invalid value for parameter \"STREAMING_MASK_OCL_INTERNAL_PRECISION\"");
  }
}

}
}
