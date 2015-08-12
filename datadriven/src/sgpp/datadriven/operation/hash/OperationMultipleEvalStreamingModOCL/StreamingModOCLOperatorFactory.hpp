/*
 * OCLOperatorFactory.hpp
 *
 *  Created on: Mar 25, 2015
 *      Author: pfandedd
 */

#pragma once

#include <sgpp/base/operation/hash/OperationMultipleEval.hpp>
#include <sgpp/base/exception/factory_exception.hpp>
#include <sgpp/base/opencl/OCLConfigurationParameters.hpp>
#include <sgpp/globaldef.hpp>
#include "OperationMultiEvalStreamingModOCL.hpp"

namespace SGPP {
  namespace datadriven {

    base::OperationMultipleEval* createStreamingModOCLConfigured(base::Grid& grid, base::DataMatrix& dataset, base::OCLConfigurationParameters parameters);

  }
}
