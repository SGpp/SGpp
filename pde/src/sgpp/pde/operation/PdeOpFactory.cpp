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
#include <sgpp/pde/operation/hash/OperationLaplacePoly.hpp>
#include <sgpp/pde/operation/hash/OperationLaplaceModPoly.hpp>
#include <sgpp/pde/operation/hash/OperationLaplacePolyBoundary.hpp>
#include <sgpp/pde/operation/hash/OperationLaplacePolyClenshawCurtis.hpp>
#include <sgpp/pde/operation/hash/OperationLaplacePolyClenshawCurtisBoundary.hpp>
#include <sgpp/pde/operation/hash/OperationLaplaceModPolyClenshawCurtis.hpp>
#include <sgpp/pde/operation/hash/OperationLaplaceBspline.hpp>
#include <sgpp/pde/operation/hash/OperationLaplaceModBspline.hpp>
#include <sgpp/pde/operation/hash/OperationLaplaceBsplineBoundary.hpp>
#include <sgpp/pde/operation/hash/OperationLaplaceBsplineClenshawCurtis.hpp>
#include <sgpp/pde/operation/hash/OperationLaplaceModBsplineClenshawCurtis.hpp>
#include <sgpp/pde/operation/hash/OperationLaplaceExplicitLinear.hpp>
#include <sgpp/pde/operation/hash/OperationLaplaceExplicitBspline.hpp>
#include <sgpp/pde/operation/hash/OperationLaplaceExplicitModBspline.hpp>

#include <sgpp/pde/operation/hash/OperationMatrixLTwoDotExplicitPeriodic.hpp>
#include <sgpp/pde/operation/hash/OperationMatrixLTwoDotExplicitLinear.hpp>
#include <sgpp/pde/operation/hash/OperationMatrixLTwoDotExplicitLinearBoundary.hpp>
#include <sgpp/pde/operation/hash/OperationMatrixLTwoDotExplicitModLinear.hpp>
#include <sgpp/pde/operation/hash/OperationMatrixLTwoDotExplicitKinkLinear.hpp>
#include <sgpp/pde/operation/hash/OperationMatrixLTwoDotExplicitBspline.hpp>
#include <sgpp/pde/operation/hash/OperationMatrixLTwoDotExplicitModBspline.hpp>
#include <sgpp/pde/operation/hash/OperationMatrixLTwoDotExplicitBsplineBoundary.hpp>
#include <sgpp/pde/operation/hash/OperationMatrixLTwoDotExplicitBsplineClenshawCurtis.hpp>
#include <sgpp/pde/operation/hash/OperationMatrixLTwoDotExplicitModBsplineClenshawCurtis.hpp>
#include <sgpp/pde/operation/hash/OperationMatrixLTwoDotExplicitPoly.hpp>
#include <sgpp/pde/operation/hash/OperationMatrixLTwoDotExplicitPolyBoundary.hpp>
#include <sgpp/pde/operation/hash/OperationMatrixLTwoDotExplicitModPoly.hpp>
#include <sgpp/pde/operation/hash/OperationMatrixLTwoDotExplicitPolyClenshawCurtis.hpp>
#include <sgpp/pde/operation/hash/OperationMatrixLTwoDotExplicitPolyClenshawCurtisBoundary.hpp>
#include <sgpp/pde/operation/hash/OperationMatrixLTwoDotExplicitModPolyClenshawCurtis.hpp>

#include <sgpp/pde/operation/hash/OperationMatrixLTwoDotPeriodic.hpp>
#include <sgpp/pde/operation/hash/OperationMatrixLTwoDotModLinear.hpp>
#include <sgpp/pde/operation/hash/OperationMatrixLTwoDotKinkLinear.hpp>
#include <sgpp/pde/operation/hash/OperationMatrixLTwoDotBspline.hpp>
#include <sgpp/pde/operation/hash/OperationMatrixLTwoDotBsplineBoundary.hpp>
#include <sgpp/pde/operation/hash/OperationMatrixLTwoDotModBspline.hpp>
#include <sgpp/pde/operation/hash/OperationMatrixLTwoDotBsplineClenshawCurtis.hpp>
#include <sgpp/pde/operation/hash/OperationMatrixLTwoDotModBsplineClenshawCurtis.hpp>
#include <sgpp/pde/operation/hash/OperationMatrixLTwoDotPoly.hpp>
#include <sgpp/pde/operation/hash/OperationMatrixLTwoDotPolyBoundary.hpp>
#include <sgpp/pde/operation/hash/OperationMatrixLTwoDotModPoly.hpp>
#include <sgpp/pde/operation/hash/OperationMatrixLTwoDotPolyClenshawCurtis.hpp>
#include <sgpp/pde/operation/hash/OperationMatrixLTwoDotPolyClenshawCurtisBoundary.hpp>
#include <sgpp/pde/operation/hash/OperationMatrixLTwoDotModPolyClenshawCurtis.hpp>

#include <sgpp/pde/operation/hash/OperationLTwoDotProductLinear.hpp>
#include <sgpp/pde/operation/hash/OperationLTwoDotProductLinearBoundary.hpp>
#include <sgpp/pde/operation/hash/OperationLTwoDotProductLinearStretched.hpp>
#include <sgpp/pde/operation/hash/OperationLTwoDotProductLinearStretchedBoundary.hpp>

#include <sgpp/pde/operation/hash/OperationLaplaceEnhancedLinear.hpp>
#include <sgpp/pde/operation/hash/OperationLaplaceEnhancedLinearBoundary.hpp>

#include <sgpp/globaldef.hpp>

#include <cstring>

namespace sgpp {

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
        &(dynamic_cast<base::PrewaveletGrid&>(grid).getShadowStorage()));
  } else if (grid.getType() == base::GridType::LinearStretched) {
    return new pde::OperationLaplaceLinearStretched(&grid.getStorage());
  } else if (grid.getType() == base::GridType::LinearStretchedBoundary) {
    return new pde::OperationLaplaceLinearStretchedBoundary(&grid.getStorage());
  } else if (grid.getType() == base::GridType::Poly) {
    return new pde::OperationLaplacePoly(&grid);
  } else if (grid.getType() == base::GridType::PolyBoundary) {
    return new pde::OperationLaplacePolyBoundary(&grid);
  } else if (grid.getType() == base::GridType::ModPoly) {
    return new pde::OperationLaplaceModPoly(&grid);
  } else if (grid.getType() == base::GridType::PolyClenshawCurtis) {
    return new pde::OperationLaplacePolyClenshawCurtis(&grid);
  } else if (grid.getType() == base::GridType::PolyClenshawCurtisBoundary) {
    return new pde::OperationLaplacePolyClenshawCurtisBoundary(&grid);
  } else if (grid.getType() == base::GridType::ModPolyClenshawCurtis) {
    return new pde::OperationLaplaceModPolyClenshawCurtis(&grid);
  } else if (grid.getType() == base::GridType::Bspline) {
    return new pde::OperationLaplaceBspline(&grid);
  } else if (grid.getType() == base::GridType::ModBspline) {
    return new pde::OperationLaplaceModBspline(&grid);
  } else if (grid.getType() == base::GridType::BsplineBoundary) {
    return new pde::OperationLaplaceBsplineBoundary(&grid);
  } else if (grid.getType() == base::GridType::BsplineClenshawCurtis) {
    return new pde::OperationLaplaceBsplineClenshawCurtis(&grid);
  } else if (grid.getType() == base::GridType::ModBsplineClenshawCurtis) {
    return new pde::OperationLaplaceModBsplineClenshawCurtis(&grid);
  } else {
    throw base::factory_exception("OperationLaplace is not implemented for this grid type.");
  }
}

base::OperationMatrix* createOperationLaplace(base::Grid& grid, sgpp::base::DataVector& coef) {
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

base::OperationMatrix* createOperationLaplaceExplicit(base::Grid& grid) {
  if (grid.getType() == base::GridType::Linear) {
    return new pde::OperationLaplaceExplicitLinear(&grid.getStorage());
  } else if (grid.getType() == base::GridType::Bspline) {
    return new pde::OperationLaplaceExplicitBspline(&grid);
  } else if (grid.getType() == base::GridType::ModBspline) {
    return new pde::OperationLaplaceExplicitModBspline(&grid);
  } else {
    throw base::factory_exception(
        "OperationLaplaceExplicit is not implemented for this grid type.");
  }
}

base::OperationMatrix* createOperationLaplaceExplicit(base::DataMatrix* m, base::Grid& grid) {
  if (grid.getType() == base::GridType::Linear) {
    return new pde::OperationLaplaceExplicitLinear(m, &grid.getStorage());
  } else if (grid.getType() == base::GridType::Bspline) {
    return new pde::OperationLaplaceExplicitBspline(m, &grid);
  } else if (grid.getType() == base::GridType::ModBspline) {
    return new pde::OperationLaplaceExplicitModBspline(m, &grid);
  } else {
    throw base::factory_exception(
        "OperationLaplaceExplicit is not implemented for this grid type.");
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
  } else if (grid.getType() == base::GridType::ModLinear) {
    return new pde::OperationMatrixLTwoDotModLinear(&grid);
  } else if (grid.getType() == base::GridType::KinkLinear) {
    return new pde::OperationMatrixLTwoDotKinkLinear(&grid);
  } else if (grid.getType() == base::GridType::Poly) {
    return new pde::OperationMatrixLTwoDotPoly(&grid);
  } else if (grid.getType() == base::GridType::PolyBoundary) {
    return new pde::OperationMatrixLTwoDotPolyBoundary(&grid);
  } else if (grid.getType() == base::GridType::ModPoly) {
    return new pde::OperationMatrixLTwoDotModPoly(&grid);
  } else if (grid.getType() == base::GridType::PolyClenshawCurtis) {
    return new pde::OperationMatrixLTwoDotPolyClenshawCurtis(&grid);
  } else if (grid.getType() == base::GridType::PolyClenshawCurtisBoundary) {
    return new pde::OperationMatrixLTwoDotPolyClenshawCurtisBoundary(&grid);
  } else if (grid.getType() == base::GridType::ModPolyClenshawCurtis) {
    return new pde::OperationMatrixLTwoDotModPolyClenshawCurtis(&grid);
  } else if (grid.getType() == base::GridType::Bspline) {
    return new pde::OperationMatrixLTwoDotBspline(&grid);
  } else if (grid.getType() == base::GridType::BsplineBoundary) {
    return new pde::OperationMatrixLTwoDotBsplineBoundary(&grid);
  } else if (grid.getType() == base::GridType::ModBspline) {
    return new pde::OperationMatrixLTwoDotModBspline(&grid);
  } else if (grid.getType() == base::GridType::BsplineClenshawCurtis) {
    return new pde::OperationMatrixLTwoDotBsplineClenshawCurtis(&grid);
  } else if (grid.getType() == base::GridType::ModBsplineClenshawCurtis) {
    return new pde::OperationMatrixLTwoDotModBsplineClenshawCurtis(&grid);
  } else {
    throw base::factory_exception("OperationLTwoDotProduct is not implemented for this grid type.");
  }
}

base::OperationMatrix* createOperationLTwoDotExplicit(base::Grid& grid) {
  if (grid.getType() == base::GridType::Linear) {
    return new pde::OperationMatrixLTwoDotExplicitLinear(&grid);
  } else if (grid.getType() == base::GridType::LinearL0Boundary ||
             grid.getType() == base::GridType::LinearBoundary) {
    return new pde::OperationMatrixLTwoDotExplicitLinearBoundary(&grid);
  } else if (grid.getType() == base::GridType::ModLinear) {
    return new pde::OperationMatrixLTwoDotExplicitModLinear(&grid);
  } else if (grid.getType() == base::GridType::KinkLinear) {
    return new pde::OperationMatrixLTwoDotExplicitKinkLinear(&grid);
  } else if (grid.getType() == base::GridType::Periodic) {
    return new pde::OperationMatrixLTwoDotExplicitPeriodic(&grid);
  } else if (grid.getType() == base::GridType::Bspline) {
    return new pde::OperationMatrixLTwoDotExplicitBspline(&grid);
  } else if (grid.getType() == base::GridType::BsplineBoundary) {
    return new pde::OperationMatrixLTwoDotExplicitBsplineBoundary(&grid);
  } else if (grid.getType() == base::GridType::ModBspline) {
    return new pde::OperationMatrixLTwoDotExplicitModBspline(&grid);
  } else if (grid.getType() == base::GridType::BsplineClenshawCurtis) {
    return new pde::OperationMatrixLTwoDotExplicitBsplineClenshawCurtis(&grid);
  } else if (grid.getType() == base::GridType::ModBsplineClenshawCurtis) {
    return new pde::OperationMatrixLTwoDotExplicitModBsplineClenshawCurtis(&grid);
  } else if (grid.getType() == base::GridType::Poly) {
    return new pde::OperationMatrixLTwoDotExplicitPoly(&grid);
  } else if (grid.getType() == base::GridType::PolyBoundary) {
    return new pde::OperationMatrixLTwoDotExplicitPolyBoundary(&grid);
  } else if (grid.getType() == base::GridType::ModPoly) {
    return new pde::OperationMatrixLTwoDotExplicitModPoly(&grid);
  } else if (grid.getType() == base::GridType::PolyClenshawCurtis) {
    return new pde::OperationMatrixLTwoDotExplicitPolyClenshawCurtis(&grid);
  } else if (grid.getType() == base::GridType::PolyClenshawCurtisBoundary) {
    return new pde::OperationMatrixLTwoDotExplicitPolyClenshawCurtisBoundary(&grid);
  } else if (grid.getType() == base::GridType::ModPolyClenshawCurtis) {
    return new pde::OperationMatrixLTwoDotExplicitModPolyClenshawCurtis(&grid);
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
  } else if (grid.getType() == base::GridType::ModLinear) {
    return new pde::OperationMatrixLTwoDotExplicitModLinear(m, &grid);
  } else if (grid.getType() == base::GridType::KinkLinear) {
    return new pde::OperationMatrixLTwoDotExplicitKinkLinear(m, &grid);
  } else if (grid.getType() == base::GridType::Periodic) {
    return new pde::OperationMatrixLTwoDotExplicitPeriodic(m, &grid);
  } else if (grid.getType() == base::GridType::Bspline) {
    return new pde::OperationMatrixLTwoDotExplicitBspline(m, &grid);
  } else if (grid.getType() == base::GridType::BsplineBoundary) {
    return new pde::OperationMatrixLTwoDotExplicitBsplineBoundary(m, &grid);
  } else if (grid.getType() == base::GridType::ModBspline) {
    return new pde::OperationMatrixLTwoDotExplicitModBspline(m, &grid);
  } else if (grid.getType() == base::GridType::BsplineClenshawCurtis) {
    return new pde::OperationMatrixLTwoDotExplicitBsplineClenshawCurtis(m, &grid);
  } else if (grid.getType() == base::GridType::ModBsplineClenshawCurtis) {
    return new pde::OperationMatrixLTwoDotExplicitModBsplineClenshawCurtis(m, &grid);
  } else if (grid.getType() == base::GridType::Poly) {
    return new pde::OperationMatrixLTwoDotExplicitPoly(m, &grid);
  } else if (grid.getType() == base::GridType::PolyBoundary) {
    return new pde::OperationMatrixLTwoDotExplicitPolyBoundary(m, &grid);
  } else if (grid.getType() == base::GridType::ModPoly) {
    return new pde::OperationMatrixLTwoDotExplicitModPoly(m, &grid);
  } else if (grid.getType() == base::GridType::PolyClenshawCurtis) {
    return new pde::OperationMatrixLTwoDotExplicitPolyClenshawCurtis(m, &grid);
  } else if (grid.getType() == base::GridType::PolyClenshawCurtisBoundary) {
    return new pde::OperationMatrixLTwoDotExplicitPolyClenshawCurtisBoundary(m, &grid);
  } else if (grid.getType() == base::GridType::ModPolyClenshawCurtis) {
    return new pde::OperationMatrixLTwoDotExplicitModPolyClenshawCurtis(m, &grid);
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
                                                      sgpp::base::DataVector& coef) {
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
}  // namespace sgpp
