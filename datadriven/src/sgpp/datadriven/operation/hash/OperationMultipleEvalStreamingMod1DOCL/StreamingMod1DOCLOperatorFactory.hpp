/*
 * OCLOperatorFactory.hpp
 *
 *  Created on: Mar 25, 2015
 *      Author: pfandedd
 */

#pragma once

#include <sgpp/base/operation/hash/OperationMultipleEval.hpp>
#include <sgpp/base/exception/factory_exception.hpp>
#include <sgpp/base/opencl/OpenCLConfigurationParameters.hpp>
#include <sgpp/globaldef.hpp>
#include "OperationMultiEvalStreamingMod1DOCL.hpp"

namespace SGPP {
namespace datadriven {

base::OperationMultipleEval *createStreamingMod1DOCLConfigured(base::Grid& grid, base::DataMatrix& dataset);

}
}
