// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/operation/hash/OperationQuadrature.hpp>
#include <sgpp/combigrid/functions/WeightFunctionsCollection.hpp>

#include <sgpp/globaldef.hpp>
#include "../../../../../../base/src/sgpp/base/operation/hash/common/basis/NakBsplineModifiedBasis.hpp"

namespace sgpp {
namespace combigrid {

/**
 * Weighted Quadrature on sparse grid, not a knot B-spline grid created by transformation of a not a
 * knot B-spline combigrid In UQ the weights are probability density functions.
 * It is designed analogously to the operations in BaseOpFactory
 */
class OperationWeightedQuadratureNakBsplineModified : public sgpp::base::OperationQuadrature {
 public:
  /**
   * Constructor of OperationWeightedQuadratureNakBsplineModified
   *
   * @param storage Pointer to the grid's GridStorage object
   * @param degree the B-spline degree
   */
  OperationWeightedQuadratureNakBsplineModified(sgpp::base::GridStorage& storage, size_t degree,
                                                WeightFunctionsCollection weightFunctionsCollection,
                                                sgpp::base::DataVector bounds)
      : storage(storage),
        base(degree),
        weightFunctionsCollection(weightFunctionsCollection),
        bounds(bounds) {}

  ~OperationWeightedQuadratureNakBsplineModified() override {}

  /**
   * Quadrature for weighted not a knot B-spline basis functions
   *
   * @param alpha Coefficient vector for current grid
   */
  double doQuadrature(sgpp::base::DataVector& alpha) override;

 protected:
  // Pointer to the grid's GridStorage object
  sgpp::base::GridStorage& storage;
  /// NakBsplineBoundaryCombigrid Basis object
  sgpp::base::SNakBsplineModifiedBase base;
  WeightFunctionsCollection weightFunctionsCollection;
  sgpp::base::DataVector bounds;
};

}  // namespace combigrid
}  // namespace sgpp
