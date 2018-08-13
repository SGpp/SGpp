// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifndef SGPP_OPTIMIZATION_OPERATION_HASH_OPMULTHIERFUNDAMENTALNOTAKNOTSPLINEBOUNDARY_HPP
#define SGPP_OPTIMIZATION_OPERATION_HASH_OPMULTHIERFUNDAMENTALNOTAKNOTSPLINEBOUNDARY_HPP

#include <sgpp/globaldef.hpp>

#include <sgpp/optimization/operation/hash/OperationMultipleHierarchisation.hpp>
#include <sgpp/base/grid/type/FundamentalNotAKnotSplineBoundaryGrid.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/operation/hash/OperationHierarchisationFundamentalNotAKnotSplineBoundary.hpp>

namespace sgpp {
namespace optimization {

/**
 * Hierarchisation operation for B-spline basis functions on
 * Noboundary grids.
 */
class OperationMultipleHierarchisationFundamentalNotAKnotSplineBoundary :
    public OperationMultipleHierarchisation {
 public:
  /**
   * Constructor.
   *
   * @param grid      grid
   */
  explicit OperationMultipleHierarchisationFundamentalNotAKnotSplineBoundary(
      base::FundamentalNotAKnotSplineBoundaryGrid& grid);

  /**
   * Destructor.
   */
  ~OperationMultipleHierarchisationFundamentalNotAKnotSplineBoundary() override;

  /**
   * @param[in,out] nodeValues before: vector of function values at
   *                           the grid points,
   *                           after: vector of hierarchical coefficients
   * @return                   whether hierarchisation was successful
   */
  bool doHierarchisation(base::DataVector& nodeValues) override;

  /**
   * @param[in,out] alpha before: vector of hierarchical coefficients,
   *                      after: vector of function values at
   *                      the grid points
   */
  void doDehierarchisation(base::DataVector& alpha) override;

  /**
   * @param[in,out] nodeValues before: matrix of function values at
   *                           the grid points,
   *                           after: matrix of hierarchical coefficients
   * @return                   whether hierarchisation was successful
   */
  bool doHierarchisation(base::DataMatrix& nodeValues) override;

  /**
   * @param[in,out] alpha before: matrix of hierarchical coefficients,
   *                      after: matrix of function values at
   *                      the grid points
   */
  void doDehierarchisation(base::DataMatrix& alpha) override;

 protected:
  /// storage of the sparse grid
  base::FundamentalNotAKnotSplineBoundaryGrid& grid;
  /// hierarchization operation
  base::OperationHierarchisationFundamentalNotAKnotSplineBoundary op;
};
}  // namespace optimization
}  // namespace sgpp

#endif /* SGPP_OPTIMIZATION_OPERATION_HASH_OPMULTHIERFUNDAMENTALNOTAKNOTSPLINEBOUNDARY_HPP */
