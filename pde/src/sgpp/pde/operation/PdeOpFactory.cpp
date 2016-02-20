// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/pde/operation/PdeOpFactory.hpp>

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

#include <cstring>

namespace SGPP {

namespace op_factory {

base::OperationMatrix* createOperationLaplace(base::Grid& grid) {
  if (grid.getType() == base::GridType::Linear) {
    return new pde::OperationLaplaceLinear(&grid.getStorage());
  } else if (grid.getType() == base::GridType::LinearL0Boundary ||
             grid.getType() == base::GridType::LinearBoundary) {
    return new pde::OperationLaplaceLinearBoundary(&grid.getStorage());
  } else if (grid.getType() == base::GridType::ModLinear) {
    return new pde::OperationLaplaceModLinear(&grid.getStorage());
  } else if (grid.getType() == base::GridType::Prewavelet) {
    return new pde::OperationLaplacePrewavelet(&grid.getStorage(),
                                               ((base::PrewaveletGrid*)&grid)->getShadowStorage());
  } else if (grid.getType() == base::GridType::LinearStretched) {
    return new pde::OperationLaplaceLinearStretched(&grid.getStorage());
  } else if (grid.getType() == base::GridType::LinearStretchedBoundary) {
    return new pde::OperationLaplaceLinearStretchedBoundary(&grid.getStorage());
  } else {
    throw base::factory_exception("OperationLaplace is not implemented for this grid type.");
  }
}

base::OperationMatrix* createOperationLaplace(base::Grid& grid, SGPP::base::DataVector& coef) {
  if (grid.getType() == base::GridType::Linear) {
    return new pde::OperationLaplaceLinear(&grid.getStorage(), coef);
  } else if (grid.getType() == base::GridType::LinearL0Boundary ||
             grid.getType() == base::GridType::LinearBoundary) {
    return new pde::OperationLaplaceLinearBoundary(&grid.getStorage(), coef);
  } else {
    throw base::factory_exception(
        "OperationLaplace (with coefficients) is not implemented for this grid type.");
  }
}

base::OperationMatrix* createOperationLTwoDotProduct(base::Grid& grid) {
  if (grid.getType() == base::GridType::Linear) {
    return new pde::OperationLTwoDotProductLinear(&grid.getStorage());
  } else if (grid.getType() == base::GridType::LinearL0Boundary ||
             grid.getType() == base::GridType::LinearBoundary) {
    return new pde::OperationLTwoDotProductLinearBoundary(&grid.getStorage());
  } else if (grid.getType() == base::GridType::LinearStretched) {
    return new pde::OperationLTwoDotProductLinearStretched(&grid.getStorage());
  } else if (grid.getType() == base::GridType::LinearStretchedBoundary) {
    return new pde::OperationLTwoDotProductLinearStretchedBoundary(&grid.getStorage());
  } else if (grid.getType() == base::GridType::Periodic) {
    return new pde::OperationMatrixLTwoDotPeriodic(&grid.getStorage());
  } else {
    throw base::factory_exception("OperationLTwoDotProduct is not implemented for this grid type.");
  }
}

base::OperationMatrix* createOperationLTwoDotExplicit(base::Grid& grid) {
  if (grid.getType() == base::GridType::Linear) {
    return new pde::OperationMatrixLTwoDotExplicitLinear(&grid);
  } else if (grid.getType() == base::GridType::LinearL0Boundary) {
    return new pde::OperationMatrixLTwoDotExplicitLinearBoundary(&grid);
  } else if (grid.getType() == base::GridType::Periodic) {
    return new pde::OperationMatrixLTwoDotExplicitPeriodic(&grid);
  } else {
    throw base::factory_exception(
        "OperationLTwoDotExplicit is not implemented for this grid type.");
  }
}

base::OperationMatrix* createOperationLTwoDotExplicit(base::DataMatrix* m, base::Grid& grid) {
  if (grid.getType() == base::GridType::Linear) {
    return new pde::OperationMatrixLTwoDotExplicitLinear(m, &grid);
  } else if (grid.getType() == base::GridType::LinearL0Boundary) {
    return new pde::OperationMatrixLTwoDotExplicitLinearBoundary(m, &grid);
  } else if (grid.getType() == base::GridType::Periodic) {
    return new pde::OperationMatrixLTwoDotExplicitPeriodic(m, &grid);
  } else {
    throw base::factory_exception(
        "OperationLTwoDotExplicit is not implemented for this grid type.");
  }
}

base::OperationMatrix* createOperationLaplaceEnhanced(base::Grid& grid) {
  if (grid.getType() == base::GridType::Linear) {
    return new pde::OperationLaplaceEnhancedLinear(&grid.getStorage());
  } else if (grid.getType() == base::GridType::LinearL0Boundary ||
             grid.getType() == base::GridType::LinearBoundary) {
    return new pde::OperationLaplaceEnhancedLinearBoundary(&grid.getStorage());
  } else {
    throw base::factory_exception(
        "OperationLaplaceEnhanced is not implemented for this grid type.");
  }
}

base::OperationMatrix* createOperationLaplaceEnhanced(base::Grid& grid,
                                                      SGPP::base::DataVector& coef) {
  if (grid.getType() == base::GridType::Linear) {
    return new pde::OperationLaplaceEnhancedLinear(&grid.getStorage(), coef);
  } else if (grid.getType() == base::GridType::LinearL0Boundary ||
             grid.getType() == base::GridType::LinearBoundary) {
    return new pde::OperationLaplaceEnhancedLinearBoundary(&grid.getStorage(), coef);
  } else {
    throw base::factory_exception(
        "OperationLaplaceEnhanced is not implemented for this grid type.");
  }
}
}  // namespace op_factory
}  // namespace SGPP
