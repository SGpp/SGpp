/*
 * OCLOperatorFactory.hpp
 *
 *  Created on: Mar 25, 2015
 *      Author: pfandedd
 */

#pragma once

#include <sgpp/base/operation/hash/OperationMultipleEval.hpp>
#include "OperationMultiEvalStreamingOCL.hpp"
#include "StreamingOCLParameters.hpp"

#include <sgpp/globaldef.hpp>

namespace SGPP {
namespace datadriven {

base::OperationMultipleEval *createStreamingOCLConfigured(base::Grid& grid,
		base::DataMatrix& dataset) {
	return new datadriven::OperationMultiEvalStreamingOCL<STREAMING_OCL_INTERNAL_PRECISION>(grid, dataset);
}

}
}
