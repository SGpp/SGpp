// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/pde/operation/PdeOpFactory.hpp>

#include <cstring>

#include <sgpp/base/exception/factory_exception.hpp>

#include <sgpp/base/grid/type/PrewaveletGrid.hpp>

#include <sgpp/pde/operation/hash/OperationLaplaceLinear.hpp>
#include <sgpp/pde/operation/hash/OperationLaplaceLinearBoundary.hpp>
#include <sgpp/pde/operation/hash/OperationLaplaceModLinear.hpp>
#include <sgpp/pde/operation/hash/OperationLaplacePrewavelet.hpp>
#include <sgpp/pde/operation/hash/OperationLaplaceLinearStretched.hpp>
#include <sgpp/pde/operation/hash/OperationLaplaceLinearStretchedBoundary.hpp>

#include <sgpp/pde/operation/hash/OperationLTwoDotProductLinear.hpp>
#include <sgpp/pde/operation/hash/OperationMatrixLTwoDotExplicitLinear.hpp>
#include <sgpp/pde/operation/hash/OperationLTwoDotProductLinearBoundary.hpp>
#include <sgpp/pde/operation/hash/OperationMatrixLTwoDotExplicitLinearBoundary.hpp>
#include <sgpp/pde/operation/hash/OperationLTwoDotProductLinearStretched.hpp>
#include <sgpp/pde/operation/hash/OperationLTwoDotProductLinearStretchedBoundary.hpp>
#include <sgpp/pde/operation/hash/OperationMatrixLTwoDotExplicitPeriodic.hpp>
#include <sgpp/pde/operation/hash/OperationMatrixLTwoDotPeriodic.hpp>

#include <sgpp/pde/operation/hash/OperationLaplaceEnhancedLinear.hpp>
#include <sgpp/pde/operation/hash/OperationLaplaceEnhancedLinearBoundary.hpp>

#include <sgpp/globaldef.hpp>


namespace SGPP {

  namespace op_factory {

    base::OperationMatrix* createOperationLaplace(base::Grid& grid) {

      if (strcmp(grid.getType(), "linear") == 0) {
        return new pde::OperationLaplaceLinear(grid.getStorage());
      } else if (strcmp(grid.getType(), "linearL0Boundary") == 0 || strcmp(grid.getType(), "linearBoundary") == 0) {
        return new pde::OperationLaplaceLinearBoundary(grid.getStorage());
      } else if (strcmp(grid.getType(), "modlinear") == 0 ) {
        return new pde::OperationLaplaceModLinear(grid.getStorage());
      } else if (strcmp(grid.getType(), "prewavelet") == 0 ) {
        return new pde::OperationLaplacePrewavelet(grid.getStorage(),
               ((base::PrewaveletGrid*) &grid)->getShadowStorage());
      } else if (strcmp(grid.getType(), "linearStretched") == 0 ) {
        return new pde::OperationLaplaceLinearStretched(grid.getStorage());
      } else if (strcmp(grid.getType(), "LinearStretchedBoundary") == 0 ) {
        return new pde::OperationLaplaceLinearStretchedBoundary(grid.getStorage());
      } else {
        throw base::factory_exception("OperationLaplace is not implemented for this grid type.");
      }
    }

    base::OperationMatrix* createOperationLaplace(base::Grid& grid, SGPP::base::DataVector& coef) {

      if (strcmp(grid.getType(), "linear") == 0) {
        return new pde::OperationLaplaceLinear(grid.getStorage(), coef);
      } else if (strcmp(grid.getType(), "linearL0Boundary") == 0 || strcmp(grid.getType(), "linearBoundary") == 0) {
        return new pde::OperationLaplaceLinearBoundary(grid.getStorage(), coef);
      } else {
        throw base::factory_exception("OperationLaplace (with coefficients) is not implemented for this grid type.");
      }
    }

    base::OperationMatrix* createOperationLTwoDotProduct(base::Grid& grid) {

      if (strcmp(grid.getType(), "linear") == 0) {
        return new pde::OperationLTwoDotProductLinear(grid.getStorage());
      } else if (strcmp(grid.getType(), "linearL0Boundary") == 0 || strcmp(grid.getType(), "linearBoundary") == 0) {
        return new pde::OperationLTwoDotProductLinearBoundary(grid.getStorage());
      } else if (strcmp(grid.getType(), "linearStretched") == 0) {
        return new pde::OperationLTwoDotProductLinearStretched(grid.getStorage());
      } else if (strcmp(grid.getType(), "LinearStretchedBoundary") == 0) {
        return new pde::OperationLTwoDotProductLinearStretchedBoundary(grid.getStorage());
      } else if (strcmp(grid.getType(), "periodic") == 0) {
        return new pde::OperationMatrixLTwoDotPeriodic(grid.getStorage());
      } else
        throw base::factory_exception("OperationLTwoDotProduct is not implemented for this grid type.");
    }

    base::OperationMatrix* createOperationLTwoDotExplicit(base::Grid& grid) {
      if (strcmp(grid.getType(), "linear") == 0) {
        return new pde::OperationMatrixLTwoDotExplicitLinear(&grid);
      } else if (strcmp(grid.getType(), "linearL0Boundary") == 0) {
        return new pde::OperationMatrixLTwoDotExplicitLinearBoundary(&grid);
      } else if (strcmp(grid.getType(), "periodic") == 0) {
        return new pde::OperationMatrixLTwoDotExplicitPeriodic(&grid);
      } else
        throw base::factory_exception("OperationLTwoDotExplicit is not implemented for this grid type.");
    }

    base::OperationMatrix* createOperationLTwoDotExplicit(base::DataMatrix* m, base::Grid& grid) {
      if (strcmp(grid.getType(), "linear") == 0) {
        return new pde::OperationMatrixLTwoDotExplicitLinear(m, &grid);
      } else if (strcmp(grid.getType(), "linearL0Boundary") == 0) {
        return new pde::OperationMatrixLTwoDotExplicitLinearBoundary(m, &grid);
      } else if (strcmp(grid.getType(), "periodic") == 0) {
        return new pde::OperationMatrixLTwoDotExplicitPeriodic(m, &grid);
      } else
        throw base::factory_exception("OperationLTwoDotExplicit is not implemented for this grid type.");
    }

    base::OperationMatrix* createOperationLaplaceEnhanced(base::Grid& grid) {

      if (strcmp(grid.getType(), "linear") == 0) {
        return new pde::OperationLaplaceEnhancedLinear(grid.getStorage());
      } else if (strcmp(grid.getType(), "linearL0Boundary") == 0 || strcmp(grid.getType(), "linearBoundary") == 0) {
        return new pde::OperationLaplaceEnhancedLinearBoundary(grid.getStorage());
      } else {
        throw base::factory_exception("OperationLaplaceEnhanced is not implemented for this grid type.");
      }
    }

    base::OperationMatrix* createOperationLaplaceEnhanced(base::Grid& grid, SGPP::base::DataVector& coef) {

      if (strcmp(grid.getType(), "linear") == 0) {
        return new pde::OperationLaplaceEnhancedLinear(grid.getStorage(), coef);
      } else if (strcmp(grid.getType(), "linearL0Boundary") == 0 || strcmp(grid.getType(), "linearBoundary") == 0) {
        return new pde::OperationLaplaceEnhancedLinearBoundary(grid.getStorage(), coef);
      } else {
        throw base::factory_exception("OperationLaplaceEnhanced is not implemented for this grid type.");
      }
    }

  }
}
