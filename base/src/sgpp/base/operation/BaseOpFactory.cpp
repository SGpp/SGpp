// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/operation/BaseOpFactory.hpp>

#include <sgpp/base/exception/factory_exception.hpp>

#include <sgpp/base/grid/type/BsplineBoundaryGrid.hpp>
#include <sgpp/base/grid/type/BsplineClenshawCurtisGrid.hpp>
#include <sgpp/base/grid/type/BsplineGrid.hpp>
#include <sgpp/base/grid/type/FundamentalSplineGrid.hpp>
#include <sgpp/base/grid/type/ModBsplineClenshawCurtisGrid.hpp>
#include <sgpp/base/grid/type/ModBsplineGrid.hpp>
#include <sgpp/base/grid/type/ModFundamentalSplineGrid.hpp>
#include <sgpp/base/grid/type/ModPolyGrid.hpp>
#include <sgpp/base/grid/type/PolyBoundaryGrid.hpp>
#include <sgpp/base/grid/type/PolyGrid.hpp>
#include <sgpp/base/grid/type/PolyClenshawCurtisBoundaryGrid.hpp>
#include <sgpp/base/grid/type/PolyClenshawCurtisGrid.hpp>

#include <sgpp/base/grid/type/PrewaveletGrid.hpp>

#include <sgpp/base/operation/hash/OperationHierarchisationFundamentalSpline.hpp>
#include <sgpp/base/operation/hash/OperationHierarchisationLinear.hpp>
#include <sgpp/base/operation/hash/OperationHierarchisationLinearBoundary.hpp>
#include <sgpp/base/operation/hash/OperationHierarchisationLinearStretched.hpp>
#include <sgpp/base/operation/hash/OperationHierarchisationLinearStretchedBoundary.hpp>
#include <sgpp/base/operation/hash/OperationHierarchisationModBspline.hpp>
#include <sgpp/base/operation/hash/OperationHierarchisationModFundamentalSpline.hpp>
#include <sgpp/base/operation/hash/OperationHierarchisationModLinear.hpp>
#include <sgpp/base/operation/hash/OperationHierarchisationModPoly.hpp>
#include <sgpp/base/operation/hash/OperationHierarchisationModWavelet.hpp>
#include <sgpp/base/operation/hash/OperationHierarchisationPoly.hpp>
#include <sgpp/base/operation/hash/OperationHierarchisationPolyBoundary.hpp>
#include <sgpp/base/operation/hash/OperationHierarchisationPolyClenshawCurtis.hpp>
#include <sgpp/base/operation/hash/OperationHierarchisationPolyClenshawCurtisBoundary.hpp>
#include <sgpp/base/operation/hash/OperationHierarchisationPrewavelet.hpp>
#include <sgpp/base/operation/hash/OperationStencilHierarchisationLinear.hpp>
#include <sgpp/base/operation/hash/OperationStencilHierarchisationModLinear.hpp>
#include <sgpp/base/operation/hash/OperationHierarchisationLinearClenshawCurtis.hpp>
#include <sgpp/base/operation/hash/OperationHierarchisationLinearClenshawCurtisBoundary.hpp>

#include <sgpp/base/operation/hash/OperationFirstMomentLinear.hpp>
#include <sgpp/base/operation/hash/OperationQuadrature.hpp>
#include <sgpp/base/operation/hash/OperationQuadratureLinear.hpp>
#include <sgpp/base/operation/hash/OperationQuadratureLinearBoundary.hpp>
#include <sgpp/base/operation/hash/OperationQuadraturePoly.hpp>
#include <sgpp/base/operation/hash/OperationQuadratureModPoly.hpp>
#include <sgpp/base/operation/hash/OperationQuadraturePolyBoundary.hpp>
#include <sgpp/base/operation/hash/OperationQuadraturePolyClenshawCurtis.hpp>
#include <sgpp/base/operation/hash/OperationQuadratureBspline.hpp>
#include <sgpp/base/operation/hash/OperationQuadratureBsplineBoundary.hpp>
#include <sgpp/base/operation/hash/OperationQuadratureBsplineClenshawCurtis.hpp>
#include <sgpp/base/operation/hash/OperationQuadratureModBspline.hpp>
#include <sgpp/base/operation/hash/OperationQuadratureFundamentalSpline.hpp>
#include <sgpp/base/operation/hash/OperationQuadratureModFundamentalSpline.hpp>

#include <sgpp/base/operation/hash/OperationSecondMomentLinear.hpp>

#include <sgpp/base/operation/hash/OperationConvertPrewavelet.hpp>

#include <sgpp/base/operation/hash/OperationEvalLinear.hpp>
#include <sgpp/base/operation/hash/OperationEvalLinearBoundary.hpp>
#include <sgpp/base/operation/hash/OperationEvalLinearStretched.hpp>
#include <sgpp/base/operation/hash/OperationEvalLinearStretchedBoundary.hpp>
#include <sgpp/base/operation/hash/OperationEvalModLinear.hpp>
#include <sgpp/base/operation/hash/OperationEvalModPoly.hpp>
#include <sgpp/base/operation/hash/OperationEvalPeriodic.hpp>
#include <sgpp/base/operation/hash/OperationEvalPoly.hpp>
#include <sgpp/base/operation/hash/OperationEvalPolyBoundary.hpp>
#include <sgpp/base/operation/hash/OperationEvalPrewavelet.hpp>
#include <sgpp/base/operation/hash/OperationEvalLinearClenshawCurtisNaive.hpp>
#include <sgpp/base/operation/hash/OperationEvalLinearClenshawCurtisBoundaryNaive.hpp>

#include <sgpp/base/operation/hash/OperationMultipleEvalLinear.hpp>
#include <sgpp/base/operation/hash/OperationMultipleEvalLinearBoundary.hpp>
#include <sgpp/base/operation/hash/OperationMultipleEvalLinearStretched.hpp>
#include <sgpp/base/operation/hash/OperationMultipleEvalLinearStretchedBoundary.hpp>
#include <sgpp/base/operation/hash/OperationMultipleEvalModLinear.hpp>
#include <sgpp/base/operation/hash/OperationMultipleEvalModPoly.hpp>
#include <sgpp/base/operation/hash/OperationMultipleEvalPeriodic.hpp>
#include <sgpp/base/operation/hash/OperationMultipleEvalPoly.hpp>
#include <sgpp/base/operation/hash/OperationMultipleEvalPolyBoundary.hpp>
#include <sgpp/base/operation/hash/OperationMultipleEvalPolyClenshawCurtis.hpp>
#include <sgpp/base/operation/hash/OperationMultipleEvalPolyClenshawCurtisBoundary.hpp>
#include <sgpp/base/operation/hash/OperationMultipleEvalPrewavelet.hpp>
#include <sgpp/base/operation/hash/OperationMultipleEvalLinearClenshawCurtis.hpp>
#include <sgpp/base/operation/hash/OperationMultipleEvalLinearClenshawCurtisBoundary.hpp>

#include <sgpp/base/operation/hash/OperationMultipleEvalBsplineNaive.hpp>
#include <sgpp/base/operation/hash/OperationMultipleEvalModBsplineNaive.hpp>

#include <sgpp/base/operation/hash/OperationEvalBsplineNaive.hpp>
#include <sgpp/base/operation/hash/OperationEvalBsplineBoundaryNaive.hpp>
#include <sgpp/base/operation/hash/OperationEvalBsplineClenshawCurtisNaive.hpp>
#include <sgpp/base/operation/hash/OperationEvalFundamentalSplineNaive.hpp>
#include <sgpp/base/operation/hash/OperationEvalLinearNaive.hpp>
#include <sgpp/base/operation/hash/OperationEvalLinearBoundaryNaive.hpp>
#include <sgpp/base/operation/hash/OperationEvalModBsplineNaive.hpp>
#include <sgpp/base/operation/hash/OperationEvalModBsplineClenshawCurtisNaive.hpp>
#include <sgpp/base/operation/hash/OperationEvalModFundamentalSplineNaive.hpp>
#include <sgpp/base/operation/hash/OperationEvalModLinearNaive.hpp>
#include <sgpp/base/operation/hash/OperationEvalModPolyNaive.hpp>
#include <sgpp/base/operation/hash/OperationEvalModWaveletNaive.hpp>
#include <sgpp/base/operation/hash/OperationEvalPolyNaive.hpp>
#include <sgpp/base/operation/hash/OperationEvalPolyBoundaryNaive.hpp>
#include <sgpp/base/operation/hash/OperationEvalPolyClenshawCurtisNaive.hpp>
#include <sgpp/base/operation/hash/OperationEvalPolyClenshawCurtisBoundaryNaive.hpp>
#include <sgpp/base/operation/hash/OperationEvalWaveletNaive.hpp>
#include <sgpp/base/operation/hash/OperationEvalWaveletBoundaryNaive.hpp>

#include <sgpp/base/operation/hash/OperationEvalGradientBsplineNaive.hpp>
#include <sgpp/base/operation/hash/OperationEvalGradientBsplineBoundaryNaive.hpp>
#include <sgpp/base/operation/hash/OperationEvalGradientBsplineClenshawCurtisNaive.hpp>
#include <sgpp/base/operation/hash/OperationEvalGradientFundamentalSplineNaive.hpp>
#include <sgpp/base/operation/hash/OperationEvalGradientModBsplineNaive.hpp>
#include <sgpp/base/operation/hash/OperationEvalGradientModBsplineClenshawCurtisNaive.hpp>
#include <sgpp/base/operation/hash/OperationEvalGradientModFundamentalSplineNaive.hpp>
#include <sgpp/base/operation/hash/OperationEvalGradientModWaveletNaive.hpp>
#include <sgpp/base/operation/hash/OperationEvalGradientWaveletNaive.hpp>
#include <sgpp/base/operation/hash/OperationEvalGradientWaveletBoundaryNaive.hpp>

#include <sgpp/base/operation/hash/OperationEvalHessianBsplineNaive.hpp>
#include <sgpp/base/operation/hash/OperationEvalHessianBsplineBoundaryNaive.hpp>
#include <sgpp/base/operation/hash/OperationEvalHessianBsplineClenshawCurtisNaive.hpp>
#include <sgpp/base/operation/hash/OperationEvalHessianFundamentalSplineNaive.hpp>
#include <sgpp/base/operation/hash/OperationEvalHessianModBsplineNaive.hpp>
#include <sgpp/base/operation/hash/OperationEvalHessianModBsplineClenshawCurtisNaive.hpp>
#include <sgpp/base/operation/hash/OperationEvalHessianModFundamentalSplineNaive.hpp>
#include <sgpp/base/operation/hash/OperationEvalHessianModWaveletNaive.hpp>
#include <sgpp/base/operation/hash/OperationEvalHessianWaveletNaive.hpp>
#include <sgpp/base/operation/hash/OperationEvalHessianWaveletBoundaryNaive.hpp>

#include <sgpp/base/operation/hash/OperationEvalPartialDerivativeBsplineNaive.hpp>
#include <sgpp/base/operation/hash/OperationEvalPartialDerivativeBsplineBoundaryNaive.hpp>
#include <sgpp/base/operation/hash/OperationEvalPartialDerivativeBsplineClenshawCurtisNaive.hpp>
#include <sgpp/base/operation/hash/OperationEvalPartialDerivativeFundamentalSplineNaive.hpp>
#include <sgpp/base/operation/hash/OperationEvalPartialDerivativeModBsplineNaive.hpp>
#include <sgpp/base/operation/hash/OperationEvalPartialDerivativeModBsplineClenshawCurtisNaive.hpp>
#include <sgpp/base/operation/hash/OperationEvalPartialDerivativeModFundamentalSplineNaive.hpp>
#include <sgpp/base/operation/hash/OperationEvalPartialDerivativeModWaveletNaive.hpp>
#include <sgpp/base/operation/hash/OperationEvalPartialDerivativeWaveletNaive.hpp>
#include <sgpp/base/operation/hash/OperationEvalPartialDerivativeWaveletBoundaryNaive.hpp>

#include <sgpp/globaldef.hpp>

#include <cstring>

namespace sgpp {

namespace op_factory {

base::OperationHierarchisation* createOperationHierarchisation(base::Grid& grid) {
  if (grid.getType() == base::GridType::Linear) {
    return new base::OperationHierarchisationLinear(grid.getStorage());
  } else if (grid.getType() == base::GridType::LinearStencil) {
    return new base::OperationStencilHierarchisationLinear(grid.getStorage());
  } else if (grid.getType() == base::GridType::ModLinearStencil) {
    return new base::OperationStencilHierarchisationModLinear(grid.getStorage());
  } else if (grid.getType() == base::GridType::ModLinear) {
    return new base::OperationHierarchisationModLinear(grid.getStorage());
  } else if (grid.getType() == base::GridType::LinearClenshawCurtis) {
    return new base::OperationHierarchisationLinearClenshawCurtis(grid.getStorage());
  } else if (grid.getType() == base::GridType::LinearClenshawCurtisBoundary) {
    return new base::OperationHierarchisationLinearClenshawCurtisBoundary(grid.getStorage());
  } else if (grid.getType() == base::GridType::LinearL0Boundary ||
             grid.getType() == base::GridType::LinearBoundary ||
             grid.getType() == base::GridType::LinearTruncatedBoundary ||
             grid.getType() == base::GridType::SquareRoot) {
    return new base::OperationHierarchisationLinearBoundary(grid.getStorage());
  } else if (grid.getType() == base::GridType::LinearStretched) {
    return new base::OperationHierarchisationLinearStretched(grid.getStorage());
  } else if (grid.getType() == base::GridType::LinearStretchedBoundary) {
    return new base::OperationHierarchisationLinearStretchedBoundary(grid.getStorage());
  } else if (grid.getType() == base::GridType::Poly) {
    return new base::OperationHierarchisationPoly(
        grid.getStorage(), dynamic_cast<base::PolyGrid*>(&grid)->getDegree());
  } else if (grid.getType() == base::GridType::PolyBoundary) {
    return new base::OperationHierarchisationPolyBoundary(
        grid.getStorage(), dynamic_cast<base::PolyBoundaryGrid*>(&grid)->getDegree());
  } else if (grid.getType() == base::GridType::ModPoly) {
    return new base::OperationHierarchisationModPoly(
        grid.getStorage(), dynamic_cast<base::ModPolyGrid*>(&grid)->getDegree());
  } else if (grid.getType() == base::GridType::PolyClenshawCurtis) {
    return new base::OperationHierarchisationPolyClenshawCurtis(
        grid.getStorage(), dynamic_cast<base::PolyClenshawCurtisGrid*>(&grid)->getDegree());
  } else if (grid.getType() == base::GridType::PolyClenshawCurtisBoundary) {
    return new base::OperationHierarchisationPolyClenshawCurtisBoundary(
        grid.getStorage(), dynamic_cast<base::PolyClenshawCurtisBoundaryGrid*>(&grid)->getDegree());
  } else if (grid.getType() == base::GridType::Prewavelet) {
    return new base::OperationHierarchisationPrewavelet(
        grid.getStorage(), dynamic_cast<base::PrewaveletGrid*>(&grid)->getShadowStorage());
  } else if (grid.getType() == base::GridType::FundamentalSpline) {
    return new base::OperationHierarchisationFundamentalSpline(
        dynamic_cast<base::FundamentalSplineGrid*>(&grid));
  } else if (grid.getType() == base::GridType::ModFundamentalSpline) {
    return new base::OperationHierarchisationModFundamentalSpline(
        dynamic_cast<base::ModFundamentalSplineGrid*>(&grid));
  } else {
    throw base::factory_exception(
        "createOperationHierarchisation is not implemented for this grid type.");
  }
}

base::OperationQuadrature* createOperationQuadrature(base::Grid& grid) {
  if (grid.getType() == base::GridType::Linear) {
    return new base::OperationQuadratureLinear(grid.getStorage());
  } else if (grid.getType() == base::GridType::LinearL0Boundary ||
             grid.getType() == base::GridType::LinearBoundary) {
    return new base::OperationQuadratureLinearBoundary(grid.getStorage());
  } else if (grid.getType() == base::GridType::Poly) {
    return new base::OperationQuadraturePoly(grid.getStorage(),
                                             dynamic_cast<base::PolyGrid*>(&grid)->getDegree());
  } else if (grid.getType() == base::GridType::ModPoly) {
    return new base::OperationQuadratureModPoly(
        grid.getStorage(), dynamic_cast<base::ModPolyGrid*>(&grid)->getDegree());
  } else if (grid.getType() == base::GridType::PolyBoundary) {
    return new base::OperationQuadraturePolyBoundary(
        grid.getStorage(), dynamic_cast<base::PolyBoundaryGrid*>(&grid)->getDegree());
  } else if (grid.getType() == base::GridType::PolyClenshawCurtisBoundary) {
    return new base::OperationQuadraturePolyClenshawCurtis(
        grid.getStorage(), dynamic_cast<base::PolyClenshawCurtisBoundaryGrid*>(&grid)->getDegree());
  } else if (grid.getType() == base::GridType::Bspline) {
    return new base::OperationQuadratureBspline(
        grid.getStorage(), dynamic_cast<base::BsplineGrid*>(&grid)->getDegree());
  } else if (grid.getType() == base::GridType::BsplineBoundary) {
    return new base::OperationQuadratureBsplineBoundary(
        grid.getStorage(), dynamic_cast<base::BsplineBoundaryGrid*>(&grid)->getDegree());
  } else if (grid.getType() == base::GridType::BsplineClenshawCurtis) {
    return new base::OperationQuadratureBsplineClenshawCurtis(
        grid.getStorage(), dynamic_cast<base::BsplineClenshawCurtisGrid*>(&grid)->getDegree());
  } else if (grid.getType() == base::GridType::ModBspline) {
    return new base::OperationQuadratureModBspline(
        grid.getStorage(), dynamic_cast<base::ModBsplineGrid*>(&grid)->getDegree());
  } else if (grid.getType() == base::GridType::FundamentalSpline) {
    return new base::OperationQuadratureFundamentalSpline(
        grid.getStorage(), dynamic_cast<base::FundamentalSplineGrid*>(&grid)->getDegree());
  } else if (grid.getType() == base::GridType::ModFundamentalSpline) {
    return new base::OperationQuadratureModFundamentalSpline(
        grid.getStorage(), dynamic_cast<base::ModFundamentalSplineGrid*>(&grid)->getDegree());
  } else {
    throw base::factory_exception(
        "createOperationQuadrature is not implemented for this grid type.");
  }
}

base::OperationFirstMoment* createOperationFirstMoment(base::Grid& grid) {
  if (grid.getType() == base::GridType::Linear) {
    return new base::OperationFirstMomentLinear(grid.getStorage());
  } else {
    throw base::factory_exception(
        "createOperationFirstMoment is not implemented for this grid type.");
  }
}

base::OperationSecondMoment* createOperationSecondMoment(base::Grid& grid) {
  if (grid.getType() == base::GridType::Linear) {
    return new base::OperationSecondMomentLinear(grid.getStorage());
  } else {
    throw base::factory_exception(
        "createOperationSecondMoment is not implemented for this grid type.");
  }
}

base::OperationConvert* createOperationConvert(base::Grid& grid) {
  if (grid.getType() == base::GridType::Prewavelet) {
    return new base::OperationConvertPrewavelet(grid.getStorage(),
                                                ((base::PrewaveletGrid*)&grid)->getShadowStorage());
  } else {
    throw base::factory_exception("createOperationConvert is not implemented for this grid type.");
  }
}

base::OperationMatrix* createOperationIdentity(base::Grid& grid) {
  return new base::OperationIdentity();
}

base::OperationEval* createOperationEval(base::Grid& grid) {
  if (grid.getType() == base::GridType::Linear) {
    return new base::OperationEvalLinear(grid.getStorage());
  } else if (grid.getType() == base::GridType::LinearL0Boundary ||
             grid.getType() == base::GridType::LinearBoundary ||
             grid.getType() == base::GridType::LinearTruncatedBoundary ||
             grid.getType() == base::GridType::SquareRoot) {
    return new base::OperationEvalLinearBoundary(grid.getStorage());
  } else if (grid.getType() == base::GridType::ModLinear) {
    return new base::OperationEvalModLinear(grid.getStorage());
  } else if (grid.getType() == base::GridType::Poly) {
    return new base::OperationEvalPoly(grid.getStorage(),
                                       dynamic_cast<base::PolyGrid*>(&grid)->getDegree());
  } else if (grid.getType() == base::GridType::PolyBoundary) {
    return new base::OperationEvalPolyBoundary(
        grid.getStorage(), dynamic_cast<base::PolyBoundaryGrid*>(&grid)->getDegree());
  } else if (grid.getType() == base::GridType::ModPoly) {
    return new base::OperationEvalModPoly(grid.getStorage(),
                                          dynamic_cast<base::ModPolyGrid*>(&grid)->getDegree());
  } else if (grid.getType() == base::GridType::Prewavelet) {
    return new base::OperationEvalPrewavelet(grid.getStorage());
  } else if (grid.getType() == base::GridType::LinearStretched) {
    return new base::OperationEvalLinearStretched(grid.getStorage());
  } else if (grid.getType() == base::GridType::LinearStretchedBoundary) {
    return new base::OperationEvalLinearStretchedBoundary(grid.getStorage());
  } else if (grid.getType() == base::GridType::Periodic) {
    return new base::OperationEvalPeriodic(grid.getStorage());
  } else {
    throw base::factory_exception(
        "createOperationEval is not implemented for this grid type. "
        "Try createOperationEvalNaive instead.");
  }
}

base::OperationMultipleEval* createOperationMultipleEval(base::Grid& grid,
                                                         base::DataMatrix& dataset) {
  if (grid.getType() == base::GridType::Linear) {
    return new base::OperationMultipleEvalLinear(grid, dataset);
  } else if (grid.getType() == base::GridType::LinearL0Boundary ||
             grid.getType() == base::GridType::LinearBoundary) {
    return new base::OperationMultipleEvalLinearBoundary(grid, dataset);
  } else if (grid.getType() == base::GridType::ModLinear) {
    return new base::OperationMultipleEvalModLinear(grid, dataset);
  } else if (grid.getType() == base::GridType::Poly) {
    return new base::OperationMultipleEvalPoly(
        grid, dynamic_cast<base::PolyGrid*>(&grid)->getDegree(), dataset);
  } else if (grid.getType() == base::GridType::PolyBoundary) {
    return new base::OperationMultipleEvalPolyBoundary(
        grid, dynamic_cast<base::PolyBoundaryGrid*>(&grid)->getDegree(), dataset);
  } else if (grid.getType() == base::GridType::ModPoly) {
    return new base::OperationMultipleEvalModPoly(
        grid, dynamic_cast<base::ModPolyGrid*>(&grid)->getDegree(), dataset);
  } else if (grid.getType() == base::GridType::PolyClenshawCurtis) {
    return new base::OperationMultipleEvalPolyClenshawCurtis(
        grid, dynamic_cast<base::PolyClenshawCurtisGrid*>(&grid)->getDegree(), dataset);
  } else if (grid.getType() == base::GridType::PolyClenshawCurtisBoundary) {
    return new base::OperationMultipleEvalPolyClenshawCurtisBoundary(
        grid, dynamic_cast<base::PolyClenshawCurtisBoundaryGrid*>(&grid)->getDegree(), dataset);
  } else if (grid.getType() == base::GridType::LinearClenshawCurtis) {
    return new base::OperationMultipleEvalLinearClenshawCurtis(grid, dataset);
  } else if (grid.getType() == base::GridType::LinearClenshawCurtisBoundary) {
    return new base::OperationMultipleEvalLinearClenshawCurtisBoundary(grid, dataset);
  } else if (grid.getType() == base::GridType::Prewavelet) {
    return new base::OperationMultipleEvalPrewavelet(grid, dataset);
  } else if (grid.getType() == base::GridType::LinearStretched) {
    return new base::OperationMultipleEvalLinearStretched(grid, dataset);
  } else if (grid.getType() == base::GridType::LinearStretchedBoundary) {
    return new base::OperationMultipleEvalLinearStretchedBoundary(grid, dataset);
  } else if (grid.getType() == base::GridType::Periodic) {
    return new base::OperationMultipleEvalPeriodic(grid, dataset);
  } else {
    throw base::factory_exception(
        "createOperationMultipleEval is not implemented for this grid type.");
  }
}

base::OperationMultipleEval* createOperationMultipleEvalNaive(base::Grid& grid,
                                                              base::DataMatrix& dataset) {
  if (grid.getType() == base::GridType::Bspline) {
    return new base::OperationMultipleEvalBsplineNaive(
        grid, dynamic_cast<base::BsplineGrid*>(&grid)->getDegree(), dataset);
  } else if (grid.getType() == base::GridType::ModBspline) {
    return new base::OperationMultipleEvalModBsplineNaive(
        grid, dynamic_cast<base::ModBsplineGrid*>(&grid)->getDegree(), dataset);
  } else {
    throw base::factory_exception(
        "createOperationMultipleEvalNaive is not implemented for this grid type.");
  }
}

base::OperationEval* createOperationEvalNaive(base::Grid& grid) {
  if (grid.getType() == base::GridType::Linear) {
    return new base::OperationEvalLinearNaive(grid.getStorage());
  } else if (grid.getType() == base::GridType::ModLinear) {
    return new base::OperationEvalModLinearNaive(grid.getStorage());
  } else if (grid.getType() == base::GridType::LinearBoundary) {
    return new base::OperationEvalLinearBoundaryNaive(grid.getStorage());
  } else if (grid.getType() == base::GridType::LinearClenshawCurtis) {
    return new base::OperationEvalLinearClenshawCurtisNaive(grid.getStorage());
  } else if (grid.getType() == base::GridType::LinearClenshawCurtisBoundary) {
    return new base::OperationEvalLinearClenshawCurtisBoundaryNaive(grid.getStorage());
  } else if (grid.getType() == base::GridType::Bspline) {
    return new base::OperationEvalBsplineNaive(grid.getStorage(),
                                               dynamic_cast<base::BsplineGrid&>(grid).getDegree());
  } else if (grid.getType() == base::GridType::ModBspline) {
    return new base::OperationEvalModBsplineNaive(
        grid.getStorage(), dynamic_cast<base::ModBsplineGrid&>(grid).getDegree());
  } else if (grid.getType() == base::GridType::ModBsplineClenshawCurtis) {
    return new base::OperationEvalModBsplineClenshawCurtisNaive(
        grid.getStorage(), dynamic_cast<base::ModBsplineClenshawCurtisGrid&>(grid).getDegree());
  } else if (grid.getType() == base::GridType::BsplineBoundary) {
    return new base::OperationEvalBsplineBoundaryNaive(
        grid.getStorage(), dynamic_cast<base::BsplineBoundaryGrid&>(grid).getDegree());
  } else if (grid.getType() == base::GridType::BsplineClenshawCurtis) {
    return new base::OperationEvalBsplineClenshawCurtisNaive(
        grid.getStorage(), dynamic_cast<base::BsplineClenshawCurtisGrid&>(grid).getDegree());
  } else if (grid.getType() == base::GridType::Wavelet) {
    return new base::OperationEvalWaveletNaive(grid.getStorage());
  } else if (grid.getType() == base::GridType::ModWavelet) {
    return new base::OperationEvalModWaveletNaive(grid.getStorage());
  } else if (grid.getType() == base::GridType::WaveletBoundary) {
    return new base::OperationEvalWaveletBoundaryNaive(grid.getStorage());
  } else if (grid.getType() == base::GridType::FundamentalSpline) {
    return new base::OperationEvalFundamentalSplineNaive(
        grid.getStorage(), dynamic_cast<base::FundamentalSplineGrid&>(grid).getDegree());
  } else if (grid.getType() == base::GridType::ModFundamentalSpline) {
    return new base::OperationEvalModFundamentalSplineNaive(
        grid.getStorage(), dynamic_cast<base::ModFundamentalSplineGrid&>(grid).getDegree());
  } else if (grid.getType() == base::GridType::Poly) {
    return new base::OperationEvalPolyNaive(grid.getStorage(),
                                            dynamic_cast<base::PolyGrid*>(&grid)->getDegree());
  } else if (grid.getType() == base::GridType::PolyBoundary) {
    return new base::OperationEvalPolyBoundaryNaive(
        grid.getStorage(), dynamic_cast<base::PolyBoundaryGrid*>(&grid)->getDegree());
  } else if (grid.getType() == base::GridType::ModPoly) {
    return new base::OperationEvalModPolyNaive(
        grid.getStorage(), dynamic_cast<base::ModPolyGrid*>(&grid)->getDegree());
  } else if (grid.getType() == base::GridType::PolyClenshawCurtisBoundary) {
    return new base::OperationEvalPolyClenshawCurtisBoundaryNaive(
        grid.getStorage(), dynamic_cast<base::PolyClenshawCurtisBoundaryGrid*>(&grid)->getDegree());
  } else if (grid.getType() == base::GridType::PolyClenshawCurtis) {
    return new base::OperationEvalPolyClenshawCurtisNaive(
        grid.getStorage(), dynamic_cast<base::PolyClenshawCurtisGrid*>(&grid)->getDegree());
  } else {
    throw base::factory_exception(
        "createOperationEvalNaive is not implemented for this grid type.");
  }
}

base::OperationEvalGradient* createOperationEvalGradientNaive(base::Grid& grid) {
  if (grid.getType() == base::GridType::Bspline) {
    return new base::OperationEvalGradientBsplineNaive(
        grid.getStorage(), dynamic_cast<base::BsplineGrid&>(grid).getDegree());
  } else if (grid.getType() == base::GridType::ModBspline) {
    return new base::OperationEvalGradientModBsplineNaive(
        grid.getStorage(), dynamic_cast<base::ModBsplineGrid&>(grid).getDegree());
  } else if (grid.getType() == base::GridType::ModBsplineClenshawCurtis) {
    return new base::OperationEvalGradientModBsplineClenshawCurtisNaive(
        grid.getStorage(), dynamic_cast<base::ModBsplineClenshawCurtisGrid&>(grid).getDegree());
  } else if (grid.getType() == base::GridType::BsplineBoundary) {
    return new base::OperationEvalGradientBsplineBoundaryNaive(
        grid.getStorage(), dynamic_cast<base::BsplineBoundaryGrid&>(grid).getDegree());
  } else if (grid.getType() == base::GridType::BsplineClenshawCurtis) {
    return new base::OperationEvalGradientBsplineClenshawCurtisNaive(
        grid.getStorage(), dynamic_cast<base::BsplineClenshawCurtisGrid&>(grid).getDegree());
  } else if (grid.getType() == base::GridType::Wavelet) {
    return new base::OperationEvalGradientWaveletNaive(grid.getStorage());
  } else if (grid.getType() == base::GridType::ModWavelet) {
    return new base::OperationEvalGradientModWaveletNaive(grid.getStorage());
  } else if (grid.getType() == base::GridType::WaveletBoundary) {
    return new base::OperationEvalGradientWaveletBoundaryNaive(grid.getStorage());
  } else if (grid.getType() == base::GridType::FundamentalSpline) {
    return new base::OperationEvalGradientFundamentalSplineNaive(
        grid.getStorage(), dynamic_cast<base::FundamentalSplineGrid&>(grid).getDegree());
  } else if (grid.getType() == base::GridType::ModFundamentalSpline) {
    return new base::OperationEvalGradientModFundamentalSplineNaive(
        grid.getStorage(), dynamic_cast<base::ModFundamentalSplineGrid&>(grid).getDegree());
  } else {
    throw base::factory_exception(
        "createOperationEvalGradient is not implemented for this grid type.");
  }
}

base::OperationEvalHessian* createOperationEvalHessianNaive(base::Grid& grid) {
  if (grid.getType() == base::GridType::Bspline) {
    return new base::OperationEvalHessianBsplineNaive(
        grid.getStorage(), dynamic_cast<base::BsplineGrid&>(grid).getDegree());
  } else if (grid.getType() == base::GridType::ModBspline) {
    return new base::OperationEvalHessianModBsplineNaive(
        grid.getStorage(), dynamic_cast<base::ModBsplineGrid&>(grid).getDegree());
  } else if (grid.getType() == base::GridType::ModBsplineClenshawCurtis) {
    return new base::OperationEvalHessianModBsplineClenshawCurtisNaive(
        grid.getStorage(), dynamic_cast<base::ModBsplineClenshawCurtisGrid&>(grid).getDegree());
  } else if (grid.getType() == base::GridType::BsplineBoundary) {
    return new base::OperationEvalHessianBsplineBoundaryNaive(
        grid.getStorage(), dynamic_cast<base::BsplineBoundaryGrid&>(grid).getDegree());
  } else if (grid.getType() == base::GridType::BsplineClenshawCurtis) {
    return new base::OperationEvalHessianBsplineClenshawCurtisNaive(
        grid.getStorage(), dynamic_cast<base::BsplineClenshawCurtisGrid&>(grid).getDegree());
  } else if (grid.getType() == base::GridType::Wavelet) {
    return new base::OperationEvalHessianWaveletNaive(grid.getStorage());
  } else if (grid.getType() == base::GridType::ModWavelet) {
    return new base::OperationEvalHessianModWaveletNaive(grid.getStorage());
  } else if (grid.getType() == base::GridType::WaveletBoundary) {
    return new base::OperationEvalHessianWaveletBoundaryNaive(grid.getStorage());
  } else if (grid.getType() == base::GridType::FundamentalSpline) {
    return new base::OperationEvalHessianFundamentalSplineNaive(
        grid.getStorage(), dynamic_cast<base::FundamentalSplineGrid&>(grid).getDegree());
  } else if (grid.getType() == base::GridType::ModFundamentalSpline) {
    return new base::OperationEvalHessianModFundamentalSplineNaive(
        grid.getStorage(), dynamic_cast<base::ModFundamentalSplineGrid&>(grid).getDegree());
  } else {
    throw base::factory_exception(
        "createOperationEvalHessian is not implemented for this grid type.");
  }
}

base::OperationEvalPartialDerivative* createOperationEvalPartialDerivativeNaive(base::Grid& grid) {
  if (grid.getType() == base::GridType::Bspline) {
    return new base::OperationEvalPartialDerivativeBsplineNaive(
        grid.getStorage(), dynamic_cast<base::BsplineGrid&>(grid).getDegree());
  } else if (grid.getType() == base::GridType::ModBspline) {
    return new base::OperationEvalPartialDerivativeModBsplineNaive(
        grid.getStorage(), dynamic_cast<base::ModBsplineGrid&>(grid).getDegree());
  } else if (grid.getType() == base::GridType::ModBsplineClenshawCurtis) {
    return new base::OperationEvalPartialDerivativeModBsplineClenshawCurtisNaive(
        grid.getStorage(), dynamic_cast<base::ModBsplineClenshawCurtisGrid&>(grid).getDegree());
  } else if (grid.getType() == base::GridType::BsplineBoundary) {
    return new base::OperationEvalPartialDerivativeBsplineBoundaryNaive(
        grid.getStorage(), dynamic_cast<base::BsplineBoundaryGrid&>(grid).getDegree());
  } else if (grid.getType() == base::GridType::BsplineClenshawCurtis) {
    return new base::OperationEvalPartialDerivativeBsplineClenshawCurtisNaive(
        grid.getStorage(), dynamic_cast<base::BsplineClenshawCurtisGrid&>(grid).getDegree());
  } else if (grid.getType() == base::GridType::Wavelet) {
    return new base::OperationEvalPartialDerivativeWaveletNaive(grid.getStorage());
  } else if (grid.getType() == base::GridType::ModWavelet) {
    return new base::OperationEvalPartialDerivativeModWaveletNaive(grid.getStorage());
  } else if (grid.getType() == base::GridType::WaveletBoundary) {
    return new base::OperationEvalPartialDerivativeWaveletBoundaryNaive(grid.getStorage());
  } else if (grid.getType() == base::GridType::FundamentalSpline) {
    return new base::OperationEvalPartialDerivativeFundamentalSplineNaive(
        grid.getStorage(), dynamic_cast<base::FundamentalSplineGrid&>(grid).getDegree());
  } else if (grid.getType() == base::GridType::ModFundamentalSpline) {
    return new base::OperationEvalPartialDerivativeModFundamentalSplineNaive(
        grid.getStorage(), dynamic_cast<base::ModFundamentalSplineGrid&>(grid).getDegree());
  } else {
    throw base::factory_exception(
        "createOperationEvalPartialDerivative is not implemented for "
        "this grid type.");
  }
}

}  // namespace op_factory
}  // namespace sgpp
