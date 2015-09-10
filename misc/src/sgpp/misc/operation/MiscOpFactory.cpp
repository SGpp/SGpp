// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <cstring>
#include <sgpp/base/exception/factory_exception.hpp>

#include <sgpp/misc/operation/MiscOpFactory.hpp>
#include <sgpp/pde/operation/hash/OperationLaplaceEnhancedLinear.hpp>
#include <sgpp/pde/operation/hash/OperationLaplaceEnhancedLinearBoundary.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {

  namespace op_factory {
    base::OperationMatrix* createOperationLaplaceEnhanced(base::Grid& grid) {

      if (grid.getType() == base::GridType::Linear) {
        return new pde::OperationLaplaceEnhancedLinear(grid.getStorage());
      } else if (grid.getType() == base::GridType::LinearL0Boundary || grid.getType() == base::GridType::LinearBoundary) {
        return new pde::OperationLaplaceEnhancedLinearBoundary(grid.getStorage());
      } else {
        throw base::factory_exception("OperationLaplaceEnhanced is not implemented for this grid type.");
      }
    }

    base::OperationMatrix* createOperationLaplaceEnhanced(base::Grid& grid, SGPP::base::DataVector& coef) {

      if (grid.getType() == base::GridType::Linear) {
        return new pde::OperationLaplaceEnhancedLinear(grid.getStorage(), coef);
      } else if (grid.getType() == base::GridType::LinearL0Boundary || grid.getType() == base::GridType::LinearBoundary) {
        return new pde::OperationLaplaceEnhancedLinearBoundary(grid.getStorage(), coef);
      } else {
        throw base::factory_exception("OperationLaplaceEnhanced is not implemented for this grid type.");
      }
    }
  }
}
