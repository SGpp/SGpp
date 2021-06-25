// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <omp.h>

#include <memory>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/operation/hash/OperationWeightedQuadrature.hpp>
#include <sgpp/base/operation/hash/common/basis/NakPBsplineBasis.hpp>
#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace base {

/**
 * Quadrature on sparse grid
 */
class OperationWeightedQuadratureNakPBspline : public OperationWeightedQuadrature {
 public:
  /**
   * Constructor of OperationQuadratureNakPBspline
   *
   * @param storage Pointer to the grid's GridStorage object
   * @param degree the B-spline degree
   * @param quadOrder quadrature order
   */
  OperationWeightedQuadratureNakPBspline(GridStorage& storage, size_t degree, size_t quadOrder)
      : storage(storage), base(degree), quadOrder(quadOrder) {
    base::DataVector temp_quadCoordinates, temp_quadWeights;
    base::GaussLegendreQuadRule1D gauss;
    gauss.getLevelPointsAndWeightsNormalized(quadOrder, temp_quadCoordinates, temp_quadWeights);
    quadCoordinates = std::make_shared<sgpp::base::DataVector>(temp_quadCoordinates);
    quadWeights = std::make_shared<sgpp::base::DataVector>(temp_quadWeights);
  }

  ~OperationWeightedQuadratureNakPBspline() override {}

  /**
   * Quadrature for not a knot B-spline basis functions w.r.t. a probability density function
   *
   * @param alpha   	Coefficient vector for current grid
   * @param pdfs			probability density functions
   */
  double doWeightedQuadrature(DataVector& alpha, sgpp::base::DistributionsVector pdfs) override;

 protected:
  // Pointer to the grid's GridStorage object
  GridStorage& storage;
  /// basis object
  SNakPBsplineBase base;
  /// quadrature rule order
  size_t quadOrder;
  /// quadrature rule coordinates
  std::shared_ptr<base::DataVector> quadCoordinates;
  /// quadrature rule weights
  std::shared_ptr<base::DataVector> quadWeights;
};

}  // namespace base
}  // namespace sgpp
