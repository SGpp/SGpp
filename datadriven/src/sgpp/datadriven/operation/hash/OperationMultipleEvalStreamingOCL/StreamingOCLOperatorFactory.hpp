/*
 * OCLOperatorFactory.hpp
 *
 *  Created on: Mar 25, 2015
 *      Author: pfandedd
 */

#pragma once

#include <sgpp/base/operation/hash/OperationMultipleEval.hpp>
#include <sgpp/base/exception/factory_exception.hpp>
#include <sgpp/base/tools/ConfigurationParameters.hpp>
#include "OperationMultiEvalStreamingOCL.hpp"
#include "StreamingOCLParameters.hpp"

#include <sgpp/globaldef.hpp>

namespace SGPP {
namespace datadriven {

base::OperationMultipleEval *createStreamingOCLConfigured(base::Grid& grid,
		base::DataMatrix& dataset) {

	std::map<std::string, std::string> defaultParameter;
	defaultParameter["STREAMING_OCL_ENABLE_OPTIMIZATIONS"] = "true";
	defaultParameter["STREAMING_OCL_USE_LOCAL_MEMORY"] = "true";
	defaultParameter["STREAMING_OCL_LOCAL_SIZE"] = "64";
	defaultParameter["STREAMING_OCL_MAX_DIM_UNROLL"] = "10";
	defaultParameter["STREAMING_OCL_INTERNAL_PRECISION"] = "float";

	base::ConfigurationParameters parameters("StreamingOCL.cfg", defaultParameter);

	std::cout << "are optimizations on: " << parameters.getAsBoolean("STREAMING_OCL_ENABLE_OPTIMIZATIONS") << std::endl;
	std::cout << "is local memory on: " << parameters.getAsBoolean("STREAMING_OCL_USE_LOCAL_MEMORY") << std::endl;
	std::cout << "local size: " << parameters.getAsUnsigned("STREAMING_OCL_LOCAL_SIZE") << std::endl;
	std::cout << "max dim unroll: " << parameters.getAsUnsigned("STREAMING_OCL_MAX_DIM_UNROLL") << std::endl;
	std::cout << "internal precision: " << parameters["STREAMING_OCL_INTERNAL_PRECISION"] << std::endl;

	if (parameters["STREAMING_OCL_INTERNAL_PRECISION"] == "float") {
		return new datadriven::OperationMultiEvalStreamingOCL<float>(grid, dataset, parameters);
	} else if (parameters["STREAMING_OCL_INTERNAL_PRECISION"] == "double") {
		return new datadriven::OperationMultiEvalStreamingOCL<double>(grid, dataset, parameters);
	} else {
		throw base::factory_exception("Error creating operation\"OperationMultiEvalStreamingOCL\": invalid value for parameter \"STREAMING_OCL_INTERNAL_PRECISION\"");
	}
}

}
}
