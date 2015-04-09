/*
 * OCLOperatorFactory.hpp
 *
 *  Created on: Mar 25, 2015
 *      Author: pfandedd
 */

#include <sgpp/base/operation/hash/OperationMultipleEval.hpp>
#include <sgpp/base/exception/factory_exception.hpp>
#include <sgpp/base/opencl/OpenCLConfigurationParameters.hpp>
#include <sgpp/globaldef.hpp>
#include "OperationMultiEvalStreamingModOCL.hpp"

namespace SGPP {
namespace datadriven {

base::OperationMultipleEval *createStreamingModOCLConfigured(base::Grid& grid, base::DataMatrix& dataset) {

  std::map<std::string, std::string> defaultParameter;
  defaultParameter["KERNEL_USE_LOCAL_MEMORY"] = "true";
  defaultParameter["LINEAR_LOAD_BALANCING_VERBOSE"] = "false";

  base::OpenCLConfigurationParameters parameters("StreamingModOCL.cfg", defaultParameter);

//  std::cout << "are optimizations on: " << parameters.getAsBoolean("STREAMING_OCL_ENABLE_OPTIMIZATIONS") << std::endl;
//  std::cout << "is local memory on: " << parameters.getAsBoolean("STREAMING_OCL_USE_LOCAL_MEMORY") << std::endl;
//  std::cout << "local size: " << parameters.getAsUnsigned("STREAMING_OCL_LOCAL_SIZE") << std::endl;
//  std::cout << "max dim unroll: " << parameters.getAsUnsigned("STREAMING_OCL_MAX_DIM_UNROLL") << std::endl;
//  std::cout << "internal precision: " << parameters["STREAMING_OCL_INTERNAL_PRECISION"] << std::endl;
//  std::cout << "platform is: " << parameters["STREAMING_OCL_PLATFORM"] << std::endl;
//  std::cout << "device type is: " << parameters["STREAMING_OCL_DEVICE_TYPE"] << std::endl;

  if (parameters["INTERNAL_PRECISION"] == "float") {
    return new datadriven::OperationMultiEvalStreamingModOCL<float>(grid, dataset, parameters);
  } else if (parameters["INTERNAL_PRECISION"] == "double") {
    return new datadriven::OperationMultiEvalStreamingModOCL<double>(grid, dataset, parameters);
  } else {
    throw base::factory_exception(
        "Error creating operation\"OperationMultiEvalStreamingModOCL\": invalid value for parameter \"INTERNAL_PRECISION\"");
  }
}

}
}
