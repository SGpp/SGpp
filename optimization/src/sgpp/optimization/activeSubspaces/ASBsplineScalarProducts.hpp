// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once
// #ifdef USE_EIGEN

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/grid/type/NakBsplineModifiedGrid.hpp>
#include <sgpp/base/operation/hash/common/basis/NakBsplineBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/NakBsplineBoundaryBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/NakBsplineExtendedBasis.hpp>
#include <sgpp/base/operation/hash/common/basis/NakBsplineModifiedBasis.hpp>
#include <sgpp/optimization/activeSubspaces/GaussQuadrature.hpp>

#include <map>
#include <tuple>
#include <utility>

namespace sgpp {
namespace optimization {

/**
 * Calculates and stores the scalar products of B-spline functions
 */
class ASBsplineScalarProducts {
 public:
  /**
   * Constructor
   *
   * @param gridType          type of the grid for the interpolant
   * @param degree            degree for the B-spline basis functions
   */
  ASBsplineScalarProducts(sgpp::base::GridType gridType, size_t numDim, size_t degree,
                          size_t quadOrder)
      : gridType(gridType), degree(degree), quadOrder(quadOrder) {
    if (gridType == sgpp::base::GridType::NakBspline) {
      basis = std::make_unique<sgpp::base::SNakBsplineBase>(degree);
    } else if (gridType == sgpp::base::GridType::NakBsplineBoundary) {
      basis = std::make_unique<sgpp::base::SNakBsplineBoundaryBase>(degree);
    } else if (gridType == sgpp::base::GridType::NakBsplineModified) {
      basis = std::make_unique<sgpp::base::SNakBsplineModifiedBase>(degree);
    } else if (gridType == sgpp::base::GridType::NakBsplineExtended) {
      basis = std::make_unique<sgpp::base::SNakBsplineExtendedBase>(degree);
    } else {
      throw sgpp::base::generation_exception("ASMatrixNakBspline: gridType not supported.");
    }
    base::GaussLegendreQuadRule1D gauss;
    gauss.getLevelPointsAndWeightsNormalized(quadOrder, coordinates, weights);
  }

  /**
   * calculates the one diensional integral \int f*g dx where f and g are B-spline basis functions
   * or first derivatives of B-spline basis functions
   *
   * @param level1 level of the first B-spline
   * @param index1 index of the first B-spline
   * @param dx1 evaluate B-spline if False, evaluate d/dx B-spline if True
   * @param level2 level of the second B-spline
   * @param index2 index of the second B-spline
   * @param coordinates coordinates for the Gauss quadrature
   * @param weights weights for the Gauss quadrature
   * @param dx2 evaluate B-spline if False, evaluate d/dx B-spline if True
   *
   * @return  integral (derivative of) first basis function * (derivative of) second basis function
   */
  double univariateScalarProduct(size_t level1, size_t index1, bool dx1, size_t level2,
                                 size_t index2, bool dx2);

  /**
   * used to get the support segments of a b-splien basis functions. Needed for Gauss quadrature
   *
   * @param level	level of the B-spline basis function
   * @param index	index of the B-spline basis function
   *
   * @return the indices of the segments of the B-spline basis functions support
   */
  sgpp::base::DataVector nakBSplineSupport(size_t level, size_t index);

 private:
  sgpp::base::GridType gridType;
  size_t degree;
  size_t quadOrder;
  sgpp::base::DataVector coordinates;
  sgpp::base::DataVector weights;
  std::unique_ptr<sgpp::base::SBasis> basis;
  typedef std::tuple<size_t, size_t, bool, size_t, size_t, bool> asMatrixHashType;
  std::map<asMatrixHashType, double> innerProducts;
};

}  // namespace optimization
}  // namespace sgpp

// #endif /* USE_EIGEN */
