/* *****************************************************************************
* Copyright (C) 2013 Technische Universitaet Muenchen                          *
* This file is part of the SG++ project. For conditions of distribution and    *
* use, please see the copyright notice at http://www5.in.tum.de/SGpp           *
***************************************************************************** */
// @author Alexander Heinecke (alexander.heinecke@mytum.de)

#include <cstring>
#include <sgpp/base/exception/factory_exception.hpp>

#include <sgpp/misc/operation/MiscOpFactory.hpp>
#include <sgpp/misc/pde/basis/linear/noboundary/operation/OperationLaplaceEnhancedLinear.hpp>
#include <sgpp/misc/pde/basis/linear/boundary/operation/OperationLaplaceEnhancedLinearBoundary.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {

  namespace op_factory {
    base::OperationMatrix* createOperationLaplaceEnhanced(base::Grid& grid) {

      if (strcmp(grid.getType(), "linear") == 0) {
        return new pde::OperationLaplaceEnhancedLinear(grid.getStorage());
      } else if (strcmp(grid.getType(), "linearBoundary") == 0 || strcmp(grid.getType(), "linearTrapezoidBoundary") == 0) {
        return new pde::OperationLaplaceEnhancedLinearBoundary(grid.getStorage());
      } else {
        throw base::factory_exception("OperationLaplaceEnhanced is not implemented for this grid type.");
      }
    }

    base::OperationMatrix* createOperationLaplaceEnhanced(base::Grid& grid, SGPP::base::DataVector& coef) {

      if (strcmp(grid.getType(), "linear") == 0) {
        return new pde::OperationLaplaceEnhancedLinear(grid.getStorage(), coef);
      } else if (strcmp(grid.getType(), "linearBoundary") == 0 || strcmp(grid.getType(), "linearTrapezoidBoundary") == 0) {
        return new pde::OperationLaplaceEnhancedLinearBoundary(grid.getStorage(), coef);
      } else {
        throw base::factory_exception("OperationLaplaceEnhanced is not implemented for this grid type.");
      }
    }
  }
}

