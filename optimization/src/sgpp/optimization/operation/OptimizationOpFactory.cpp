// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/globaldef.hpp>

#include <sgpp/optimization/operation/OptimizationOpFactory.hpp>
#include <sgpp/base/exception/factory_exception.hpp>

#include <sgpp/base/grid/type/LinearGrid.hpp>
#include <sgpp/base/grid/type/LinearBoundaryGrid.hpp>
#include <sgpp/base/grid/type/LinearClenshawCurtisGrid.hpp>
#include <sgpp/base/grid/type/ModLinearGrid.hpp>
#include <sgpp/base/grid/type/BsplineGrid.hpp>
#include <sgpp/base/grid/type/BsplineBoundaryGrid.hpp>
#include <sgpp/base/grid/type/BsplineClenshawCurtisGrid.hpp>
#include <sgpp/base/grid/type/ModBsplineGrid.hpp>
#include <sgpp/base/grid/type/ModBsplineClenshawCurtisGrid.hpp>
#include <sgpp/base/grid/type/WaveletGrid.hpp>
#include <sgpp/base/grid/type/WaveletBoundaryGrid.hpp>
#include <sgpp/base/grid/type/ModWaveletGrid.hpp>
#include <sgpp/base/grid/type/FundamentalNotAKnotSplineBoundaryGrid.hpp>
#include <sgpp/base/grid/type/FundamentalSplineGrid.hpp>
#include <sgpp/base/grid/type/FundamentalSplineBoundaryGrid.hpp>
#include <sgpp/base/grid/type/ModFundamentalSplineGrid.hpp>

#include <sgpp/optimization/operation/hash/OperationMultipleHierarchisationBsplineBoundary.hpp>
#include <sgpp/optimization/operation/hash/OperationMultipleHierarchisationBsplineClenshawCurtis.hpp>
#include <sgpp/optimization/operation/hash/OperationMultipleHierarchisationModBspline.hpp>
#include <sgpp/optimization/operation/hash/OperationMultipleHierarchisationModBsplineClenshawCurtis.hpp>
#include <sgpp/optimization/operation/hash/OperationMultipleHierarchisationBspline.hpp>
#include <sgpp/optimization/operation/hash/OperationMultipleHierarchisationLinearBoundary.hpp>
#include <sgpp/optimization/operation/hash/OperationMultipleHierarchisationLinearClenshawCurtis.hpp>
#include <sgpp/optimization/operation/hash/OperationMultipleHierarchisationModLinear.hpp>
#include <sgpp/optimization/operation/hash/OperationMultipleHierarchisationLinear.hpp>
#include <sgpp/optimization/operation/hash/OperationMultipleHierarchisationWaveletBoundary.hpp>
#include <sgpp/optimization/operation/hash/OperationMultipleHierarchisationModWavelet.hpp>
#include <sgpp/optimization/operation/hash/OperationMultipleHierarchisationWavelet.hpp>
#include <sgpp/optimization/operation/hash/OperationMultipleHierarchisationFundamentalNotAKnotSplineBoundary.hpp>
#include <sgpp/optimization/operation/hash/OperationMultipleHierarchisationFundamentalSpline.hpp>
#include <sgpp/optimization/operation/hash/OperationMultipleHierarchisationFundamentalSplineBoundary.hpp>
#include <sgpp/optimization/operation/hash/OperationMultipleHierarchisationModFundamentalSpline.hpp>
#include <sgpp/optimization/operation/hash/OperationMultipleHierarchisationLagrangeNotAKnotSplineBoundary.hpp>
#include <sgpp/optimization/operation/hash/OperationMultipleHierarchisationLagrangeSplineBoundary.hpp>
#include <sgpp/optimization/operation/hash/OperationMultipleHierarchisationModLagrangeNotAKnotSpline.hpp>
#include <sgpp/optimization/operation/hash/OperationMultipleHierarchisationModNotAKnotBspline.hpp>
#include <sgpp/optimization/operation/hash/OperationMultipleHierarchisationNaturalBsplineBoundary.hpp>
#include <sgpp/optimization/operation/hash/OperationMultipleHierarchisationNotAKnotBsplineBoundary.hpp>

#include <cstring>

namespace sgpp {
namespace op_factory {

optimization::OperationMultipleHierarchisation*
createOperationMultipleHierarchisation(base::Grid& grid) {
  if (grid.getType() == base::GridType::Linear) {
    return new optimization::OperationMultipleHierarchisationLinear(
        dynamic_cast<base::LinearGrid&>(grid));
  } else if (grid.getType() == base::GridType::LinearBoundary) {
    return new optimization::OperationMultipleHierarchisationLinearBoundary(
        dynamic_cast<base::LinearBoundaryGrid&>(grid));
  } else if (grid.getType() == base::GridType::LinearClenshawCurtis) {
    return new optimization::OperationMultipleHierarchisationLinearClenshawCurtis(
        dynamic_cast<base::LinearClenshawCurtisGrid&>(grid));
  } else if (grid.getType() == base::GridType::ModLinear) {
    return new optimization::OperationMultipleHierarchisationModLinear(
        dynamic_cast<base::ModLinearGrid&>(grid));
  } else if (grid.getType() == base::GridType::Bspline) {
    return new optimization::OperationMultipleHierarchisationBspline(
        dynamic_cast<base::BsplineGrid&>(grid));
  } else if (grid.getType() == base::GridType::BsplineBoundary) {
    return new optimization::OperationMultipleHierarchisationBsplineBoundary(
        dynamic_cast<base::BsplineBoundaryGrid&>(grid));
  } else if (grid.getType() == base::GridType::BsplineClenshawCurtis) {
    return new optimization::OperationMultipleHierarchisationBsplineClenshawCurtis(
        dynamic_cast<base::BsplineClenshawCurtisGrid&>(grid));
  } else if (grid.getType() == base::GridType::ModBspline) {
    return new optimization::OperationMultipleHierarchisationModBspline(
        dynamic_cast<base::ModBsplineGrid&>(grid));
  } else if (grid.getType() == base::GridType::ModBsplineClenshawCurtis) {
    return new optimization::OperationMultipleHierarchisationModBsplineClenshawCurtis(
        dynamic_cast<base::ModBsplineClenshawCurtisGrid&>(grid));
  } else if (grid.getType() == base::GridType::Wavelet) {
    return new optimization::OperationMultipleHierarchisationWavelet(
        dynamic_cast<base::WaveletGrid&>(grid));
  } else if (grid.getType() == base::GridType::WaveletBoundary) {
    return new optimization::OperationMultipleHierarchisationWaveletBoundary(
        dynamic_cast<base::WaveletBoundaryGrid&>(grid));
  } else if (grid.getType() == base::GridType::ModWavelet) {
    return new optimization::OperationMultipleHierarchisationModWavelet(
        dynamic_cast<base::ModWaveletGrid&>(grid));
  } else if (grid.getType() == base::GridType::FundamentalNotAKnotSplineBoundary) {
    return new optimization::OperationMultipleHierarchisationFundamentalNotAKnotSplineBoundary(
        dynamic_cast<base::FundamentalNotAKnotSplineBoundaryGrid&>(grid));
  } else if (grid.getType() == base::GridType::FundamentalSpline) {
    return new optimization::OperationMultipleHierarchisationFundamentalSpline(
        dynamic_cast<base::FundamentalSplineGrid&>(grid));
  } else if (grid.getType() == base::GridType::FundamentalSplineBoundary) {
    return new optimization::OperationMultipleHierarchisationFundamentalSplineBoundary(
        dynamic_cast<base::FundamentalSplineBoundaryGrid&>(grid));
  } else if (grid.getType() == base::GridType::ModFundamentalSpline) {
    return new optimization::OperationMultipleHierarchisationModFundamentalSpline(
        dynamic_cast<base::ModFundamentalSplineGrid&>(grid));
  } else if (grid.getType() == base::GridType::LagrangeNotAKnotSplineBoundary) {
    return new optimization::OperationMultipleHierarchisationLagrangeNotAKnotSplineBoundary(
        dynamic_cast<base::LagrangeNotAKnotSplineBoundaryGrid&>(grid));
  } else if (grid.getType() == base::GridType::LagrangeSplineBoundary) {
    return new optimization::OperationMultipleHierarchisationLagrangeSplineBoundary(
        dynamic_cast<base::LagrangeSplineBoundaryGrid&>(grid));
  } else if (grid.getType() == base::GridType::ModLagrangeNotAKnotSpline) {
    return new optimization::OperationMultipleHierarchisationModLagrangeNotAKnotSpline(
        dynamic_cast<base::ModLagrangeNotAKnotSplineGrid&>(grid));
  } else if (grid.getType() == base::GridType::ModNotAKnotBspline) {
    return new optimization::OperationMultipleHierarchisationModNotAKnotBspline(
        dynamic_cast<base::ModNotAKnotBsplineGrid&>(grid));
  } else if (grid.getType() == base::GridType::NaturalBsplineBoundary) {
    return new optimization::OperationMultipleHierarchisationNaturalBsplineBoundary(
        dynamic_cast<base::NaturalBsplineBoundaryGrid&>(grid));
  } else if (grid.getType() == base::GridType::NotAKnotBsplineBoundary) {
    return new optimization::OperationMultipleHierarchisationNotAKnotBsplineBoundary(
        dynamic_cast<base::NotAKnotBsplineBoundaryGrid&>(grid));
  } else {
    throw base::factory_exception(
        "OperationMultipleHierarchisation is not implemented for this grid type.");
  }
}
}  // namespace op_factory
}  // namespace sgpp
