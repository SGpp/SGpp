/*
 * OCLOperatorFactory.hpp
 *
 *  Created on: Mar 25, 2015
 *      Author: pfandedd
 */

#pragma once

#include <sgpp/base/operation/hash/OperationMultipleEval.hpp>
#include <sgpp/base/exception/factory_exception.hpp>
#include <sgpp/globaldef.hpp>
#include <sgpp/base/opencl/OCLOperationConfiguration.hpp>
#include "../OperationMultipleEvalStreamingModOCLMaskMultiPlatform/OperationMultiEvalStreamingModOCLMaskMultiPlatform.hpp"

namespace SGPP {
namespace datadriven {

base::OperationMultipleEval* createStreamingModOCLMaskMultiPlatformConfigured(base::Grid& grid, base::DataMatrix& dataset,
SGPP::datadriven::OperationMultipleEvalConfiguration &configuration);

}
}
