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
#include "OperationMultiEvalStreamingModOCLFast.hpp"

namespace SGPP {
  namespace datadriven {

    base::OperationMultipleEval* createStreamingModOCLFastConfigured(base::Grid& grid, base::DataMatrix& dataset, base::OCLConfigurationParameters parameters);

  }
}
