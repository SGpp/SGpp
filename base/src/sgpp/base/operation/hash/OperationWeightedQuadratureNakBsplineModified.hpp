// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <memory>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/operation/hash/OperationWeightedQuadrature.hpp>
#include <sgpp/globaldef.hpp>
#include "common/basis/NakBsplineModifiedBasis.hpp"

namespace sgpp {
namespace base {

/**
 * Quadrature on sparse grid, not a knot B-spline grid created by transformation of a not a knot
 * B-spline combigrid
 */
class OperationWeightedQuadratureNakBsplineModified : public OperationWeightedQuadrature {
 public:
  /**
   * Constructor of OperationQuadratureNakBsplineModified
   *
   * @param storage Pointer to the grid's GridStorage object
   * @param degree the B-spline degree
   */
  OperationWeightedQuadratureNakBsplineModified(GridStorage& storage, size_t degree,
                                                size_t quadOrder)
      : storage(storage), base(degree), quadOrder(quadOrder) {
    base::DataVector temp_quadCoordinates, temp_quadWeights;
    base::GaussLegendreQuadRule1D gauss;
    gauss.getLevelPointsAndWeightsNormalized(quadOrder, temp_quadCoordinates, temp_quadWeights);
    quadCoordinates = std::make_shared<sgpp::base::DataVector>(temp_quadCoordinates);
    quadWeights = std::make_shared<sgpp::base::DataVector>(temp_quadWeights);
  }

  ~OperationWeightedQuadratureNakBsplineModified() override {}

  /**
   * Quadrature for not a knot B-spline basis functions w.r.t. a probability density function
   *
   * @param alpha   	Coefficient vector for current grid
   * @param pdf			probability density function
   * @parm quadOrder	order for the gauss Legendre quadrature
   */
  double doWeightedQuadrature(DataVector& alpha, sgpp::base::DistributionsVector pdfs);

 protected:
  // Pointer to the grid's GridStorage object
  GridStorage& storage;
  /// NakBsplineBoundaryCombigrid Basis object
  SNakBsplineModifiedBase base;
  /// quadrature rule order
  size_t quadOrder;
  /// quadrature rule coordinates
  std::shared_ptr<base::DataVector> quadCoordinates;
  /// quadrature rule weights
  std::shared_ptr<base::DataVector> quadWeights;
};

}  // namespace base
}  // namespace sgpp