// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include "QuadratureOpFactory.hpp"
#include <sgpp/globaldef.hpp>

namespace SGPP {
  namespace op_factory {

    quadrature::OperationQuadratureMCAdvanced* createOperationQuadratureMCAdvanced(base::Grid& grid, size_t numberOfSamples) {
      return new quadrature::OperationQuadratureMCAdvanced(grid, numberOfSamples);
    }

  } /* namespace op_factory */
} /* namespace SGPP */
