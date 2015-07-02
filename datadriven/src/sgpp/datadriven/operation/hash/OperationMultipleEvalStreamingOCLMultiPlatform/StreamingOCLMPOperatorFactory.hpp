/*
 * OCLOperatorFactory.hpp
 *
 *  Created on: Mar 25, 2015
 *      Author: pfandedd
 */

#pragma once

#include <sgpp/globaldef.hpp>

#include <sgpp/base/operation/hash/OperationMultipleEval.hpp>

#include <sgpp/base/opencl/OCLConfigurationParameters.hpp>

namespace SGPP {
namespace datadriven {

base::OperationMultipleEval* createStreamingOCLMPConfigured(base::Grid& grid, base::DataMatrix& dataset,
        base::OCLConfigurationParameters *parameters = nullptr);

}
}
