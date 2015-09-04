// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <cstring>

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
#include <sgpp/base/grid/type/FundamentalSplineGrid.hpp>
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
#include <sgpp/optimization/operation/hash/OperationMultipleHierarchisationFundamentalSpline.hpp>
#include <sgpp/optimization/operation/hash/OperationMultipleHierarchisationModFundamentalSpline.hpp>

namespace SGPP {
  namespace op_factory {

    optimization::OperationMultipleHierarchisation*
    createOperationMultipleHierarchisation(base::Grid& grid) {
      if (std::strcmp(grid.getType(), "linear") == 0) {
        return new optimization::OperationMultipleHierarchisationLinear(
                 dynamic_cast<base::LinearGrid&>(grid));
      } else if (std::strcmp(grid.getType(), "linearBoundary") == 0) {
        return new optimization::OperationMultipleHierarchisationLinearBoundary(
                 dynamic_cast<base::LinearBoundaryGrid&>(grid));
      } else if (std::strcmp(grid.getType(), "linearClenshawCurtis") == 0) {
        return new optimization::OperationMultipleHierarchisationLinearClenshawCurtis(
                 dynamic_cast<base::LinearClenshawCurtisGrid&>(grid));
      } else if (std::strcmp(grid.getType(), "modlinear") == 0) {
        return new optimization::OperationMultipleHierarchisationModLinear(
                 dynamic_cast<base::ModLinearGrid&>(grid));
      } else if (std::strcmp(grid.getType(), "bspline") == 0) {
        return new optimization::OperationMultipleHierarchisationBspline(
                 dynamic_cast<base::BsplineGrid&>(grid));
      } else if (std::strcmp(grid.getType(), "bsplineBoundary") == 0) {
        return new optimization::OperationMultipleHierarchisationBsplineBoundary(
                 dynamic_cast<base::BsplineBoundaryGrid&>(grid));
      } else if (std::strcmp(grid.getType(), "bsplineClenshawCurtis") == 0) {
        return new optimization::OperationMultipleHierarchisationBsplineClenshawCurtis(
                 dynamic_cast<base::BsplineClenshawCurtisGrid&>(grid));
      } else if (std::strcmp(grid.getType(), "modBspline") == 0) {
        return new optimization::OperationMultipleHierarchisationModBspline(
                 dynamic_cast<base::ModBsplineGrid&>(grid));
      } else if (std::strcmp(grid.getType(), "modBsplineClenshawCurtis") == 0) {
        return new optimization::OperationMultipleHierarchisationModBsplineClenshawCurtis(
                 dynamic_cast<base::ModBsplineClenshawCurtisGrid&>(grid));
      } else if (std::strcmp(grid.getType(), "wavelet") == 0) {
        return new optimization::OperationMultipleHierarchisationWavelet(
                 dynamic_cast<base::WaveletGrid&>(grid));
      } else if (std::strcmp(grid.getType(), "waveletBoundary") == 0) {
        return new optimization::OperationMultipleHierarchisationWaveletBoundary(
                 dynamic_cast<base::WaveletBoundaryGrid&>(grid));
      } else if (std::strcmp(grid.getType(), "modWavelet") == 0) {
        return new optimization::OperationMultipleHierarchisationModWavelet(
                 dynamic_cast<base::ModWaveletGrid&>(grid));
      } else if (std::strcmp(grid.getType(), "fundamentalSpline") == 0) {
        return new optimization::OperationMultipleHierarchisationFundamentalSpline(
                 dynamic_cast<base::FundamentalSplineGrid&>(grid));
      } else if (std::strcmp(grid.getType(), "modFundamentalSpline") == 0) {
        return new optimization::OperationMultipleHierarchisationModFundamentalSpline(
                 dynamic_cast<base::ModFundamentalSplineGrid&>(grid));
      } else {
        throw base::factory_exception(
          "OperationMultipleHierarchisation is not implemented for this grid type.");
      }
    }

  }
}
