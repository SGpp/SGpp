// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/operation/BaseOpFactory.hpp>

#include <sgpp/base/exception/factory_exception.hpp>

#include <sgpp/base/grid/type/PolyGrid.hpp>
#include <sgpp/base/grid/type/PolyBoundaryGrid.hpp>
#include <sgpp/base/grid/type/ModPolyGrid.hpp>
#include <sgpp/base/grid/type/PrewaveletGrid.hpp>
#include <sgpp/base/grid/type/BsplineGrid.hpp>
#include <sgpp/base/grid/type/BsplineBoundaryGrid.hpp>
#include <sgpp/base/grid/type/BsplineClenshawCurtisGrid.hpp>
#include <sgpp/base/grid/type/ModBsplineGrid.hpp>
#include <sgpp/base/grid/type/ModBsplineClenshawCurtisGrid.hpp>
#include <sgpp/base/grid/type/FundamentalSplineGrid.hpp>
#include <sgpp/base/grid/type/ModFundamentalSplineGrid.hpp>

#include <sgpp/base/operation/hash/OperationStencilHierarchisationLinear.hpp>
#include <sgpp/base/operation/hash/OperationStencilHierarchisationModLinear.hpp>
#include <sgpp/base/operation/hash/OperationHierarchisationLinear.hpp>
#include <sgpp/base/operation/hash/OperationHierarchisationModLinear.hpp>
#include <sgpp/base/operation/hash/OperationHierarchisationLinearBoundary.hpp>
#include <sgpp/base/operation/hash/OperationHierarchisationLinearStretched.hpp>
#include <sgpp/base/operation/hash/OperationHierarchisationLinearStretchedBoundary.hpp>
#include <sgpp/base/operation/hash/OperationHierarchisationPoly.hpp>
#include <sgpp/base/operation/hash/OperationHierarchisationPolyBoundary.hpp>
#include <sgpp/base/operation/hash/OperationHierarchisationModPoly.hpp>
#include <sgpp/base/operation/hash/OperationHierarchisationPrewavelet.hpp>
#include <sgpp/base/operation/hash/OperationHierarchisationModBspline.hpp>
#include <sgpp/base/operation/hash/OperationHierarchisationModWavelet.hpp>
#include <sgpp/base/operation/hash/OperationHierarchisationFundamentalSpline.hpp>
#include <sgpp/base/operation/hash/OperationHierarchisationModFundamentalSpline.hpp>

#include <sgpp/base/operation/hash/OperationQuadrature.hpp>
#include <sgpp/base/operation/hash/OperationQuadratureLinear.hpp>
#include <sgpp/base/operation/hash/OperationFirstMomentLinear.hpp>
#include <sgpp/base/operation/hash/OperationSecondMomentLinear.hpp>
#include <sgpp/base/operation/hash/OperationQuadratureLinearBoundary.hpp>
#include <sgpp/base/operation/hash/OperationQuadraturePoly.hpp>
#include <sgpp/base/operation/hash/OperationQuadraturePolyBoundary.hpp>

#include <sgpp/base/operation/hash/OperationConvertPrewavelet.hpp>

#include <sgpp/base/operation/hash/OperationEvalLinear.hpp>
#include <sgpp/base/operation/hash/OperationEvalLinearBoundary.hpp>
#include <sgpp/base/operation/hash/OperationEvalModLinear.hpp>
#include <sgpp/base/operation/hash/OperationEvalPoly.hpp>
#include <sgpp/base/operation/hash/OperationEvalPolyBoundary.hpp>
#include <sgpp/base/operation/hash/OperationEvalModPoly.hpp>
#include <sgpp/base/operation/hash/OperationEvalModBspline.hpp>
#include <sgpp/base/operation/hash/OperationEvalModWavelet.hpp>
#include <sgpp/base/operation/hash/OperationEvalPrewavelet.hpp>
#include <sgpp/base/operation/hash/OperationEvalLinearStretched.hpp>
#include <sgpp/base/operation/hash/OperationEvalLinearStretchedBoundary.hpp>
#include <sgpp/base/operation/hash/OperationEvalPeriodic.hpp>

#include <sgpp/base/operation/hash/OperationMultipleEvalLinear.hpp>
#include <sgpp/base/operation/hash/OperationMultipleEvalLinearBoundary.hpp>
#include <sgpp/base/operation/hash/OperationMultipleEvalModLinear.hpp>
#include <sgpp/base/operation/hash/OperationMultipleEvalPoly.hpp>
#include <sgpp/base/operation/hash/OperationMultipleEvalPolyBoundary.hpp>
#include <sgpp/base/operation/hash/OperationMultipleEvalModPoly.hpp>
#include <sgpp/base/operation/hash/OperationMultipleEvalModBspline.hpp>
#include <sgpp/base/operation/hash/OperationMultipleEvalModWavelet.hpp>
#include <sgpp/base/operation/hash/OperationMultipleEvalPrewavelet.hpp>
#include <sgpp/base/operation/hash/OperationMultipleEvalLinearStretched.hpp>
#include <sgpp/base/operation/hash/OperationMultipleEvalLinearStretchedBoundary.hpp>
#include <sgpp/base/operation/hash/OperationMultipleEvalPeriodic.hpp>

#include <sgpp/base/operation/hash/OperationNaiveEvalBsplineBoundary.hpp>
#include <sgpp/base/operation/hash/OperationNaiveEvalBsplineClenshawCurtis.hpp>
#include <sgpp/base/operation/hash/OperationNaiveEvalModBspline.hpp>
#include <sgpp/base/operation/hash/OperationNaiveEvalModBsplineClenshawCurtis.hpp>
#include <sgpp/base/operation/hash/OperationNaiveEvalBspline.hpp>
#include <sgpp/base/operation/hash/OperationNaiveEvalLinearBoundary.hpp>
#include <sgpp/base/operation/hash/OperationNaiveEvalLinearClenshawCurtis.hpp>
#include <sgpp/base/operation/hash/OperationNaiveEvalModLinear.hpp>
#include <sgpp/base/operation/hash/OperationNaiveEvalLinear.hpp>
#include <sgpp/base/operation/hash/OperationNaiveEvalWaveletBoundary.hpp>
#include <sgpp/base/operation/hash/OperationNaiveEvalModWavelet.hpp>
#include <sgpp/base/operation/hash/OperationNaiveEvalWavelet.hpp>
#include <sgpp/base/operation/hash/OperationNaiveEvalFundamentalSpline.hpp>
#include <sgpp/base/operation/hash/OperationNaiveEvalModFundamentalSpline.hpp>
#include <sgpp/base/operation/hash/OperationNaiveEvalPoly.hpp>
#include <sgpp/base/operation/hash/OperationNaiveEvalPolyBoundary.hpp>

#include <sgpp/base/operation/hash/OperationNaiveEvalGradientBsplineBoundary.hpp>
#include <sgpp/base/operation/hash/OperationNaiveEvalGradientBsplineClenshawCurtis.hpp>
#include <sgpp/base/operation/hash/OperationNaiveEvalGradientModBspline.hpp>
#include <sgpp/base/operation/hash/OperationNaiveEvalGradientModBsplineClenshawCurtis.hpp>
#include <sgpp/base/operation/hash/OperationNaiveEvalGradientBspline.hpp>
#include <sgpp/base/operation/hash/OperationNaiveEvalGradientWaveletBoundary.hpp>
#include <sgpp/base/operation/hash/OperationNaiveEvalGradientModWavelet.hpp>
#include <sgpp/base/operation/hash/OperationNaiveEvalGradientWavelet.hpp>
#include <sgpp/base/operation/hash/OperationNaiveEvalGradientFundamentalSpline.hpp>
#include <sgpp/base/operation/hash/OperationNaiveEvalGradientModFundamentalSpline.hpp>

#include <sgpp/base/operation/hash/OperationNaiveEvalHessianBsplineBoundary.hpp>
#include <sgpp/base/operation/hash/OperationNaiveEvalHessianBsplineClenshawCurtis.hpp>
#include <sgpp/base/operation/hash/OperationNaiveEvalHessianModBspline.hpp>
#include <sgpp/base/operation/hash/OperationNaiveEvalHessianModBsplineClenshawCurtis.hpp>
#include <sgpp/base/operation/hash/OperationNaiveEvalHessianBspline.hpp>
#include <sgpp/base/operation/hash/OperationNaiveEvalHessianWaveletBoundary.hpp>
#include <sgpp/base/operation/hash/OperationNaiveEvalHessianModWavelet.hpp>
#include <sgpp/base/operation/hash/OperationNaiveEvalHessianWavelet.hpp>
#include <sgpp/base/operation/hash/OperationNaiveEvalHessianFundamentalSpline.hpp>
#include <sgpp/base/operation/hash/OperationNaiveEvalHessianModFundamentalSpline.hpp>

#include <sgpp/base/operation/hash/OperationNaiveEvalPartialDerivativeBsplineBoundary.hpp>
#include <sgpp/base/operation/hash/OperationNaiveEvalPartialDerivativeBsplineClenshawCurtis.hpp>
#include <sgpp/base/operation/hash/OperationNaiveEvalPartialDerivativeModBspline.hpp>
#include <sgpp/base/operation/hash/OperationNaiveEvalPartialDerivativeModBsplineClenshawCurtis.hpp>
#include <sgpp/base/operation/hash/OperationNaiveEvalPartialDerivativeBspline.hpp>
#include <sgpp/base/operation/hash/OperationNaiveEvalPartialDerivativeWaveletBoundary.hpp>
#include <sgpp/base/operation/hash/OperationNaiveEvalPartialDerivativeModWavelet.hpp>
#include <sgpp/base/operation/hash/OperationNaiveEvalPartialDerivativeWavelet.hpp>
#include <sgpp/base/operation/hash/OperationNaiveEvalPartialDerivativeFundamentalSpline.hpp>
#include <sgpp/base/operation/hash/OperationNaiveEvalPartialDerivativeModFundamentalSpline.hpp>

#include <sgpp/globaldef.hpp>

#include <cstring>

namespace SGPP {

namespace op_factory {

std::unique_ptr<base::OperationHierarchisation> createOperationHierarchisation(
  base::Grid& grid) {
  if (grid.getType() == base::GridType::Linear) {
    return std::unique_ptr<base::OperationHierarchisation>(
        new base::OperationHierarchisationLinear(grid.getStorage()));
  } else if (grid.getType() == base::GridType::LinearStencil) {
    return std::unique_ptr<base::OperationHierarchisation>(
        new base::OperationStencilHierarchisationLinear(
             grid.getStorage()));
  } else if (grid.getType() == base::GridType::ModLinearStencil) {
    return std::unique_ptr<base::OperationHierarchisation>(
        new base::OperationStencilHierarchisationModLinear(
             grid.getStorage()));
  } else if (grid.getType() == base::GridType::ModLinear) {
    return std::unique_ptr<base::OperationHierarchisation>(
        new base::OperationHierarchisationModLinear(grid.getStorage()));
  } else if (grid.getType() == base::GridType::LinearL0Boundary
             || grid.getType() == base::GridType::LinearBoundary
             || grid.getType() == base::GridType::LinearTruncatedBoundary
             || grid.getType() == base::GridType::SquareRoot) {
    return std::unique_ptr<base::OperationHierarchisation>(
        new base::OperationHierarchisationLinearBoundary(
             grid.getStorage()));
  } else if (grid.getType() == base::GridType::LinearStretched) {
    return std::unique_ptr<base::OperationHierarchisation>(
        new base::OperationHierarchisationLinearStretched(
             grid.getStorage()));
  } else if (grid.getType() == base::GridType::LinearStretchedBoundary) {
    return std::unique_ptr<base::OperationHierarchisation>(
        new base::OperationHierarchisationLinearStretchedBoundary(
             grid.getStorage()));
  } else if (grid.getType() == base::GridType::Poly) {
    return std::unique_ptr<base::OperationHierarchisation>(
        new base::OperationHierarchisationPoly(grid.getStorage(),
           dynamic_cast<base::PolyGrid*>(&grid)->getDegree()));
  } else if (grid.getType() == base::GridType::PolyBoundary) {
    return std::unique_ptr<base::OperationHierarchisation>(
        new base::OperationHierarchisationPolyBoundary(grid.getStorage(),
           dynamic_cast<base::PolyBoundaryGrid*>(&grid)->getDegree()));
  } else if (grid.getType() == base::GridType::ModPoly) {
    return std::unique_ptr<base::OperationHierarchisation>(
        new base::OperationHierarchisationModPoly(grid.getStorage(),
           dynamic_cast<base::ModPolyGrid*>(&grid)->getDegree()));
  } else if (grid.getType() == base::GridType::Prewavelet) {
    return std::unique_ptr<base::OperationHierarchisation>(
        new base::OperationHierarchisationPrewavelet(grid.getStorage(),
           dynamic_cast<base::PrewaveletGrid*>(&grid)->getShadowStorage()));
  } else if (grid.getType() == base::GridType::FundamentalSpline) {
    return std::unique_ptr<base::OperationHierarchisation>(
        new base::OperationHierarchisationFundamentalSpline(
             dynamic_cast<base::FundamentalSplineGrid*>(&grid)));
  } else if (grid.getType() == base::GridType::ModFundamentalSpline) {
    return std::unique_ptr<base::OperationHierarchisation>(
        new base::OperationHierarchisationModFundamentalSpline(
             dynamic_cast<base::ModFundamentalSplineGrid*>(&grid)));
  } else {
    throw base::factory_exception(
      "OperationHierarchisation is not implemented for this grid type.");
  }
}

std::unique_ptr<base::OperationQuadrature> createOperationQuadrature(base::Grid& grid) {
  if (grid.getType() == base::GridType::Linear) {
    return std::unique_ptr<base::OperationQuadrature>(
        new base::OperationQuadratureLinear(grid.getStorage()));
  } else if (grid.getType() == base::GridType::LinearL0Boundary
             || grid.getType() == base::GridType::LinearBoundary) {
    return std::unique_ptr<base::OperationQuadrature>(
        new base::OperationQuadratureLinearBoundary(grid.getStorage()));
  } else if (grid.getType() == base::GridType::Poly) {
    return std::unique_ptr<base::OperationQuadrature>(
        new base::OperationQuadraturePoly(grid.getStorage(),
           dynamic_cast<base::PolyGrid*>(&grid)->getDegree()));
  } else if (grid.getType() == base::GridType::PolyBoundary) {
    return std::unique_ptr<base::OperationQuadrature>(
        new base::OperationQuadraturePolyBoundary(grid.getStorage(),
           dynamic_cast<base::PolyBoundaryGrid*>(&grid)->getDegree()));
  } else {
    throw base::factory_exception(
      "OperationQuadrature is not implemented for this grid type.");
  }
}

std::unique_ptr<base::OperationFirstMoment> createOperationFirstMoment(base::Grid& grid) {
  if (grid.getType() == base::GridType::Linear) {
    return std::unique_ptr<base::OperationFirstMoment>(
        new base::OperationFirstMomentLinear(grid.getStorage()));
  } else {
    throw base::factory_exception(
      "OperationFirstMoment is not implemented for this grid type.");
  }
}

std::unique_ptr<base::OperationSecondMoment> createOperationSecondMoment(base::Grid& grid) {
  if (grid.getType() == base::GridType::Linear) {
    return std::unique_ptr<base::OperationSecondMoment>(
        new base::OperationSecondMomentLinear(grid.getStorage()));
  } else {
    throw base::factory_exception(
      "OperationSecondMoment is not implemented for this grid type.");
  }
}

std::unique_ptr<base::OperationConvert> createOperationConvert(base::Grid& grid) {
  if (grid.getType() == base::GridType::Prewavelet) {
    return std::unique_ptr<base::OperationConvert>(
        new base::OperationConvertPrewavelet(grid.getStorage(),
           ((base::PrewaveletGrid*) &grid)->getShadowStorage()));
  } else {
    throw base::factory_exception(
      "OperationConvert is not implemented for this grid type.");
  }
}

std::unique_ptr<base::OperationMatrix> createOperationIdentity(base::Grid& grid) {
  return std::unique_ptr<base::OperationMatrix>(
        new base::OperationIdentity());
}

std::unique_ptr<base::OperationEval> createOperationEval(base::Grid& grid) {
  if (grid.getType() == base::GridType::Linear) {
    return std::unique_ptr<base::OperationEval>(
        new base::OperationEvalLinear(grid.getStorage()));
  } else if (grid.getType() == base::GridType::LinearL0Boundary
             || grid.getType() == base::GridType::LinearBoundary
             || grid.getType() == base::GridType::LinearTruncatedBoundary
             || grid.getType() == base::GridType::SquareRoot) {
    return std::unique_ptr<base::OperationEval>(
        new base::OperationEvalLinearBoundary(grid.getStorage()));
  } else if (grid.getType() == base::GridType::ModLinear) {
    return std::unique_ptr<base::OperationEval>(
        new base::OperationEvalModLinear(grid.getStorage()));
  } else if (grid.getType() == base::GridType::Poly) {
    return std::unique_ptr<base::OperationEval>(
        new base::OperationEvalPoly(
             grid.getStorage(),
             dynamic_cast<base::PolyGrid*>(&grid)->getDegree()));
  } else if (grid.getType() == base::GridType::PolyBoundary) {
    return std::unique_ptr<base::OperationEval>(
        new base::OperationEvalPolyBoundary(grid.getStorage(),
           dynamic_cast<base::PolyBoundaryGrid*>(&grid)->getDegree()));
  } else if (grid.getType() == base::GridType::ModPoly) {
    return std::unique_ptr<base::OperationEval>(
        new base::OperationEvalModPoly(
             grid.getStorage(),
             dynamic_cast<base::ModPolyGrid*>(&grid)->getDegree()));
  } else if (grid.getType() == base::GridType::ModBspline) {
    return std::unique_ptr<base::OperationEval>(
        new base::OperationEvalModBspline(grid.getStorage(),
           dynamic_cast<base::ModBsplineGrid*>(&grid)->getDegree()));
  } else if (grid.getType() == base::GridType::ModWavelet) {
    return std::unique_ptr<base::OperationEval>(
        new base::OperationEvalModWavelet(grid.getStorage()));
  } else if (grid.getType() == base::GridType::Prewavelet) {
    return std::unique_ptr<base::OperationEval>(
        new base::OperationEvalPrewavelet(grid.getStorage()));
  } else if (grid.getType() == base::GridType::LinearStretched) {
    return std::unique_ptr<base::OperationEval>(
        new base::OperationEvalLinearStretched(grid.getStorage()));
  } else if (grid.getType() == base::GridType::LinearStretchedBoundary) {
    return std::unique_ptr<base::OperationEval>(
        new base::OperationEvalLinearStretchedBoundary(grid.getStorage()));
  } else if (grid.getType() == base::GridType::Periodic) {
    return std::unique_ptr<base::OperationEval>(
        new base::OperationEvalPeriodic(grid.getStorage()));
  } else {
    throw base::factory_exception(
      "OperationEval is not implemented for this grid type.");
  }
}

std::unique_ptr<base::OperationMultipleEval> createOperationMultipleEval(base::Grid& grid,
    base::DataMatrix& dataset) {
  if (grid.getType() == base::GridType::Linear) {
    return std::unique_ptr<base::OperationMultipleEval>(
        new base::OperationMultipleEvalLinear(grid, dataset));
  } else if (grid.getType() == base::GridType::LinearL0Boundary
             || grid.getType() == base::GridType::LinearBoundary) {
    return std::unique_ptr<base::OperationMultipleEval>(
        new base::OperationMultipleEvalLinearBoundary(grid, dataset));
  } else if (grid.getType() == base::GridType::ModLinear) {
    return std::unique_ptr<base::OperationMultipleEval>(
        new base::OperationMultipleEvalModLinear(grid, dataset));
  } else if (grid.getType() == base::GridType::Poly) {
    return std::unique_ptr<base::OperationMultipleEval>(
        new base::OperationMultipleEvalPoly(grid,
           dynamic_cast<base::PolyGrid*>(&grid)->getDegree(), dataset));
  } else if (grid.getType() == base::GridType::PolyBoundary) {
    return std::unique_ptr<base::OperationMultipleEval>(
        new base::OperationMultipleEvalPolyBoundary(grid,
           dynamic_cast<base::PolyBoundaryGrid*>(&grid)->getDegree(),
           dataset));
  } else if (grid.getType() == base::GridType::ModPoly) {
    return std::unique_ptr<base::OperationMultipleEval>(
        new base::OperationMultipleEvalModPoly(grid,
           dynamic_cast<base::ModPolyGrid*>(&grid)->getDegree(), dataset));
  } else if (grid.getType() == base::GridType::ModBspline) {
    return std::unique_ptr<base::OperationMultipleEval>(
        new base::OperationMultipleEvalModBspline(grid,
           dynamic_cast<base::ModBsplineGrid*>(&grid)->getDegree(),
           dataset));
  } else if (grid.getType() == base::GridType::ModWavelet) {
    return std::unique_ptr<base::OperationMultipleEval>(
        new base::OperationMultipleEvalModWavelet(grid, dataset));
  } else if (grid.getType() == base::GridType::Prewavelet) {
    return std::unique_ptr<base::OperationMultipleEval>(
        new base::OperationMultipleEvalPrewavelet(grid, dataset));
  } else if (grid.getType() == base::GridType::LinearStretched) {
    return std::unique_ptr<base::OperationMultipleEval>(
        new base::OperationMultipleEvalLinearStretched(grid, dataset));
  } else if (grid.getType() == base::GridType::LinearStretchedBoundary) {
    return std::unique_ptr<base::OperationMultipleEval>(
        new base::OperationMultipleEvalLinearStretchedBoundary(grid,
           dataset));
  } else if (grid.getType() == base::GridType::Periodic) {
    return std::unique_ptr<base::OperationMultipleEval>(
        new base::OperationMultipleEvalPeriodic(grid, dataset));
  } else {
    throw base::factory_exception(
      "OperationMultipleEval is not implemented for this grid type.");
  }
}

std::unique_ptr<base::OperationNaiveEval> createOperationNaiveEval(base::Grid& grid) {
  if (grid.getType() == base::GridType::Linear) {
    return std::unique_ptr<base::OperationNaiveEval>(
        new base::OperationNaiveEvalLinear(grid.getStorage()));
  } else if (grid.getType() == base::GridType::ModLinear) {
    return std::unique_ptr<base::OperationNaiveEval>(
        new base::OperationNaiveEvalModLinear(grid.getStorage()));
  } else if (grid.getType() == base::GridType::LinearBoundary) {
    return std::unique_ptr<base::OperationNaiveEval>(
        new base::OperationNaiveEvalLinearBoundary(grid.getStorage()));
  } else if (grid.getType() == base::GridType::LinearClenshawCurtis) {
    return std::unique_ptr<base::OperationNaiveEval>(
        new base::OperationNaiveEvalLinearClenshawCurtis(
             grid.getStorage()));
  } else if (grid.getType() == base::GridType::Bspline) {
    return std::unique_ptr<base::OperationNaiveEval>(
        new base::OperationNaiveEvalBspline(grid.getStorage(),
           dynamic_cast<base::BsplineGrid&>(grid).getDegree()));
  } else if (grid.getType() == base::GridType::ModBspline) {
    return std::unique_ptr<base::OperationNaiveEval>(
        new base::OperationNaiveEvalModBspline(grid.getStorage(),
           dynamic_cast<base::ModBsplineGrid&>(grid).getDegree()));
  } else if (grid.getType() == base::GridType::ModBsplineClenshawCurtis) {
    return std::unique_ptr<base::OperationNaiveEval>(
        new base::OperationNaiveEvalModBsplineClenshawCurtis(
             grid.getStorage(),
             dynamic_cast<base::ModBsplineClenshawCurtisGrid&>(
               grid).getDegree()));
  } else if (grid.getType() == base::GridType::BsplineBoundary) {
    return std::unique_ptr<base::OperationNaiveEval>(
        new base::OperationNaiveEvalBsplineBoundary(grid.getStorage(),
           dynamic_cast<base::BsplineBoundaryGrid&>(grid).getDegree()));
  } else if (grid.getType() == base::GridType::BsplineClenshawCurtis) {
    return std::unique_ptr<base::OperationNaiveEval>(
        new base::OperationNaiveEvalBsplineClenshawCurtis(
             grid.getStorage(),
             dynamic_cast<base::BsplineClenshawCurtisGrid&>(grid).getDegree()));
  } else if (grid.getType() == base::GridType::Wavelet) {
    return std::unique_ptr<base::OperationNaiveEval>(
        new base::OperationNaiveEvalWavelet(grid.getStorage()));
  } else if (grid.getType() == base::GridType::ModWavelet) {
    return std::unique_ptr<base::OperationNaiveEval>(
        new base::OperationNaiveEvalModWavelet(grid.getStorage()));
  } else if (grid.getType() == base::GridType::WaveletBoundary) {
    return std::unique_ptr<base::OperationNaiveEval>(
        new base::OperationNaiveEvalWaveletBoundary(grid.getStorage()));
  } else if (grid.getType() == base::GridType::FundamentalSpline) {
    return std::unique_ptr<base::OperationNaiveEval>(
        new base::OperationNaiveEvalFundamentalSpline(grid.getStorage(),
           dynamic_cast<base::FundamentalSplineGrid&>(grid).getDegree()));
  } else if (grid.getType() == base::GridType::ModFundamentalSpline) {
    return std::unique_ptr<base::OperationNaiveEval>(
        new base::OperationNaiveEvalModFundamentalSpline(
             grid.getStorage(),
             dynamic_cast<base::ModFundamentalSplineGrid&>(grid).getDegree()));
  } else if (grid.getType() == base::GridType::Poly) {
    return std::unique_ptr<base::OperationNaiveEval>(
        new base::OperationNaiveEvalPoly(
             grid.getStorage(),
             dynamic_cast<base::PolyGrid*>(&grid)->getDegree()));
  } else if (grid.getType() == base::GridType::PolyBoundary) {
    return std::unique_ptr<base::OperationNaiveEval>(
        new base::OperationNaiveEvalPolyBoundary(grid.getStorage(),
           dynamic_cast<base::PolyBoundaryGrid*>(&grid)->getDegree()));
  } else {
    throw base::factory_exception(
      "OperationNaiveEval is not implemented for this grid type.");
  }
}

std::unique_ptr<base::OperationNaiveEvalGradient> createOperationNaiveEvalGradient(
  base::Grid& grid) {
  if (grid.getType() == base::GridType::Bspline) {
    return std::unique_ptr<base::OperationNaiveEvalGradient>(
        new base::OperationNaiveEvalGradientBspline(grid.getStorage(),
           dynamic_cast<base::BsplineGrid&>(grid).getDegree()));
  } else if (grid.getType() == base::GridType::ModBspline) {
    return std::unique_ptr<base::OperationNaiveEvalGradient>(
        new base::OperationNaiveEvalGradientModBspline(grid.getStorage(),
           dynamic_cast<base::ModBsplineGrid&>(grid).getDegree()));
  } else if (grid.getType() == base::GridType::ModBsplineClenshawCurtis) {
    return std::unique_ptr<base::OperationNaiveEvalGradient>(
        new base::OperationNaiveEvalGradientModBsplineClenshawCurtis(
             grid.getStorage(),
             dynamic_cast<base::ModBsplineClenshawCurtisGrid&>(
               grid).getDegree()));
  } else if (grid.getType() == base::GridType::BsplineBoundary) {
    return std::unique_ptr<base::OperationNaiveEvalGradient>(
        new base::OperationNaiveEvalGradientBsplineBoundary(
             grid.getStorage(),
             dynamic_cast<base::BsplineBoundaryGrid&>(grid).getDegree()));
  } else if (grid.getType() == base::GridType::BsplineClenshawCurtis) {
    return std::unique_ptr<base::OperationNaiveEvalGradient>(
        new base::OperationNaiveEvalGradientBsplineClenshawCurtis(
             grid.getStorage(),
             dynamic_cast<base::BsplineClenshawCurtisGrid&>(grid).getDegree()));
  } else if (grid.getType() == base::GridType::Wavelet) {
    return std::unique_ptr<base::OperationNaiveEvalGradient>(
        new base::OperationNaiveEvalGradientWavelet(grid.getStorage()));
  } else if (grid.getType() == base::GridType::ModWavelet) {
    return std::unique_ptr<base::OperationNaiveEvalGradient>(
        new base::OperationNaiveEvalGradientModWavelet(grid.getStorage()));
  } else if (grid.getType() == base::GridType::WaveletBoundary) {
    return std::unique_ptr<base::OperationNaiveEvalGradient>(
        new base::OperationNaiveEvalGradientWaveletBoundary(
             grid.getStorage()));
  } else if (grid.getType() == base::GridType::FundamentalSpline) {
    return std::unique_ptr<base::OperationNaiveEvalGradient>(
        new base::OperationNaiveEvalGradientFundamentalSpline(
             grid.getStorage(),
             dynamic_cast<base::FundamentalSplineGrid&>(grid).getDegree()));
  } else if (grid.getType() == base::GridType::ModFundamentalSpline) {
    return std::unique_ptr<base::OperationNaiveEvalGradient>(
        new base::OperationNaiveEvalGradientModFundamentalSpline(
             grid.getStorage(),
             dynamic_cast<base::ModFundamentalSplineGrid&>(grid).getDegree()));
  } else {
    throw base::factory_exception(
      "OperationNaiveEvalGradient is not implemented for this grid type.");
  }
}

std::unique_ptr<base::OperationNaiveEvalHessian> createOperationNaiveEvalHessian(
  base::Grid& grid) {
  if (grid.getType() == base::GridType::Bspline) {
    return std::unique_ptr<base::OperationNaiveEvalHessian>(
        new base::OperationNaiveEvalHessianBspline(grid.getStorage(),
           dynamic_cast<base::BsplineGrid&>(grid).getDegree()));
  } else if (grid.getType() == base::GridType::ModBspline) {
    return std::unique_ptr<base::OperationNaiveEvalHessian>(
        new base::OperationNaiveEvalHessianModBspline(grid.getStorage(),
           dynamic_cast<base::ModBsplineGrid&>(grid).getDegree()));
  } else if (grid.getType() == base::GridType::ModBsplineClenshawCurtis) {
    return std::unique_ptr<base::OperationNaiveEvalHessian>(
        new base::OperationNaiveEvalHessianModBsplineClenshawCurtis(
             grid.getStorage(),
             dynamic_cast<base::ModBsplineClenshawCurtisGrid&>(
               grid).getDegree()));
  } else if (grid.getType() == base::GridType::BsplineBoundary) {
    return std::unique_ptr<base::OperationNaiveEvalHessian>(
        new base::OperationNaiveEvalHessianBsplineBoundary(
             grid.getStorage(),
             dynamic_cast<base::BsplineBoundaryGrid&>(grid).getDegree()));
  } else if (grid.getType() == base::GridType::BsplineClenshawCurtis) {
    return std::unique_ptr<base::OperationNaiveEvalHessian>(
        new base::OperationNaiveEvalHessianBsplineClenshawCurtis(
             grid.getStorage(),
             dynamic_cast<base::BsplineClenshawCurtisGrid&>(grid).getDegree()));
  } else if (grid.getType() == base::GridType::Wavelet) {
    return std::unique_ptr<base::OperationNaiveEvalHessian>(
        new base::OperationNaiveEvalHessianWavelet(grid.getStorage()));
  } else if (grid.getType() == base::GridType::ModWavelet) {
    return std::unique_ptr<base::OperationNaiveEvalHessian>(
        new base::OperationNaiveEvalHessianModWavelet(grid.getStorage()));
  } else if (grid.getType() == base::GridType::WaveletBoundary) {
    return std::unique_ptr<base::OperationNaiveEvalHessian>(
        new base::OperationNaiveEvalHessianWaveletBoundary(
             grid.getStorage()));
  } else if (grid.getType() == base::GridType::FundamentalSpline) {
    return std::unique_ptr<base::OperationNaiveEvalHessian>(
        new base::OperationNaiveEvalHessianFundamentalSpline(
             grid.getStorage(),
             dynamic_cast<base::FundamentalSplineGrid&>(grid).getDegree()));
  } else if (grid.getType() == base::GridType::ModFundamentalSpline) {
    return std::unique_ptr<base::OperationNaiveEvalHessian>(
        new base::OperationNaiveEvalHessianModFundamentalSpline(
             grid.getStorage(),
             dynamic_cast<base::ModFundamentalSplineGrid&>(grid).getDegree()));
  } else {
    throw base::factory_exception(
      "OperationNaiveEvalHessian is not implemented for this grid type.");
  }
}

std::unique_ptr<base::OperationNaiveEvalPartialDerivative>
createOperationNaiveEvalPartialDerivative(
  base::Grid& grid) {
  if (grid.getType() == base::GridType::Bspline) {
    return std::unique_ptr<base::OperationNaiveEvalPartialDerivative>(
        new base::OperationNaiveEvalPartialDerivativeBspline(
             grid.getStorage(),
             dynamic_cast<base::BsplineGrid&>(grid).getDegree()));
  } else if (grid.getType() == base::GridType::ModBspline) {
    return std::unique_ptr<base::OperationNaiveEvalPartialDerivative>(
        new base::OperationNaiveEvalPartialDerivativeModBspline(
             grid.getStorage(),
             dynamic_cast<base::ModBsplineGrid&>(grid).getDegree()));
  } else if (grid.getType() == base::GridType::ModBsplineClenshawCurtis) {
    return std::unique_ptr<base::OperationNaiveEvalPartialDerivative>(
        new
           base::OperationNaiveEvalPartialDerivativeModBsplineClenshawCurtis(
             grid.getStorage(),
             dynamic_cast<base::ModBsplineClenshawCurtisGrid&>(
               grid).getDegree()));
  } else if (grid.getType() == base::GridType::BsplineBoundary) {
    return std::unique_ptr<base::OperationNaiveEvalPartialDerivative>(
        new base::OperationNaiveEvalPartialDerivativeBsplineBoundary(
             grid.getStorage(),
             dynamic_cast<base::BsplineBoundaryGrid&>(grid).getDegree()));
  } else if (grid.getType() == base::GridType::BsplineClenshawCurtis) {
    return std::unique_ptr<base::OperationNaiveEvalPartialDerivative>(
        new base::OperationNaiveEvalPartialDerivativeBsplineClenshawCurtis(
             grid.getStorage(),
             dynamic_cast<base::BsplineClenshawCurtisGrid&>(grid).getDegree()));
  } else if (grid.getType() == base::GridType::Wavelet) {
    return std::unique_ptr<base::OperationNaiveEvalPartialDerivative>(
        new base::OperationNaiveEvalPartialDerivativeWavelet(
             grid.getStorage()));
  } else if (grid.getType() == base::GridType::ModWavelet) {
    return std::unique_ptr<base::OperationNaiveEvalPartialDerivative>(
        new base::OperationNaiveEvalPartialDerivativeModWavelet(
             grid.getStorage()));
  } else if (grid.getType() == base::GridType::WaveletBoundary) {
    return std::unique_ptr<base::OperationNaiveEvalPartialDerivative>(
        new base::OperationNaiveEvalPartialDerivativeWaveletBoundary(
             grid.getStorage()));
  } else if (grid.getType() == base::GridType::FundamentalSpline) {
    return std::unique_ptr<base::OperationNaiveEvalPartialDerivative>(
        new base::OperationNaiveEvalPartialDerivativeFundamentalSpline(
             grid.getStorage(),
             dynamic_cast<base::FundamentalSplineGrid&>(grid).getDegree()));
  } else if (grid.getType() == base::GridType::ModFundamentalSpline) {
    return std::unique_ptr<base::OperationNaiveEvalPartialDerivative>(
        new base::OperationNaiveEvalPartialDerivativeModFundamentalSpline(
             grid.getStorage(),
             dynamic_cast<base::ModFundamentalSplineGrid&>(grid).getDegree()));
  } else {
    throw base::factory_exception(
      "OperationNaiveEvalPartialDerivative is not implemented for "
      "this grid type.");
  }
}

}  // namespace op_factory
}  // namespace SGPP
