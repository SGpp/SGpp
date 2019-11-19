// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/globaldef.hpp>
#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/operation/hash/OperationEval.hpp>
#include <sgpp/combigrid/grid/FullGrid.hpp>

namespace sgpp {
namespace combigrid {

/**
 * Operation for evaluating a full grid function (linear combination of full grid basis functions).
 */
class OperationEvalFullGrid : public base::OperationEval {
 public:
  /**
   * Default constructor.
   */
  OperationEvalFullGrid();

  /**
   * Constructor.
   *
   * @param grid  full grid
   */
  explicit OperationEvalFullGrid(const FullGrid& grid);

  /**
   * Virtual destructor.
   */
  ~OperationEvalFullGrid() override;

  /**
   * Evaluate a full grid function at one point.
   *
   * @param surpluses   coefficients for the full grid basis functions (may be nodal/hierarchical)
   * @param point       point at which to evaluate the full grid function
   * @return value of the full grid function at the given point
   */
  double eval(const base::DataVector& surpluses, const base::DataVector& point) override;

  /**
   * Evaluate a full grid function at multiple points.
   *
   * @param[in] surpluses   coefficients for the full grid basis functions (may be nodal/hierarchical)
   * @param[in] points      points at which to evaluate the full grid function
   *                        (every row corresponds to one point)
   * @param[out] result     values of the full grid function at the given points
   */
  virtual void multiEval(const base::DataVector& surpluses, const base::DataMatrix& points,
      base::DataVector& result);

  /**
   * @return full grid
   */
  const FullGrid& getGrid() const;

  /**
   * @param grid  full grid
   */
  void setGrid(const FullGrid& grid);

 protected:
  /// full grid
  FullGrid grid;
};

}  // namespace combigrid
}  // namespace sgpp
