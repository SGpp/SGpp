// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <omp.h>

#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/operation/hash/OperationWeightedQuadrature.hpp>
#include <sgpp/base/operation/hash/common/basis/NakBsplineExtendedBasis.hpp>
#include <sgpp/globaldef.hpp>

namespace sgpp {
namespace base {

/**
 * Quadrature on sparse grid, not a knot B-spline grid created by transformation of a not a knot
 * B-spline combigrid
 */
class OperationWeightedQuadratureNakBsplineExtended : public OperationWeightedQuadrature {
 public:
  /**
   * Constructor of OperationQuadratureNakBsplineExtended
   *
   * @param storage Pointer to the grid's GridStorage object
   * @param degree the B-spline degree
   */
  OperationWeightedQuadratureNakBsplineExtended(GridStorage& storage, size_t degree)
      : storage(storage), base(degree) {}

  ~OperationWeightedQuadratureNakBsplineExtended() override {}

  /**
   * Quadrature for not a knot B-spline basis functions w.r.t. a probability density function
   *
   * @param alpha   	Coefficient vector for current grid
   * @param pdf			probability density function
   * @parm quadOrder	order for the gauss Legendre quadrature
   */
  double doWeightedQuadrature(DataVector& alpha, sgpp::base::DistributionsVector pdfs,
                              size_t quadOrder);

 protected:
  // Pointer to the grid's GridStorage object
  GridStorage& storage;
  /// NakBsplineBoundaryCombigrid Basis object
  SNakBsplineExtendedBase base;
};

}  // namespace base
}  // namespace sgpp
