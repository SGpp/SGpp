/*
 * OCLOperatorFactory.hpp
 *
 *  Created on: Mar 25, 2015
 *      Author: pfandedd
 */

#pragma once

#include <sgpp/globaldef.hpp>

#include <sgpp/base/operation/hash/OperationMultipleEval.hpp>

#include <sgpp/base/opencl/OCLOperationConfiguration.hpp>

namespace SGPP {
  namespace datadriven {

    base::OperationMultipleEval* createStreamingOCLConfigured(base::Grid& grid, base::DataMatrix& dataset,
        SGPP::datadriven::OperationMultipleEvalConfiguration& configuration);

  }
}
