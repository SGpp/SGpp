// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/globaldef.hpp>
#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/combigrid/grid/CombinationGrid.hpp>

#include <vector>

namespace sgpp {
namespace combigrid {

/**
 * Operation for evaluating a combination grid function
 * (linear combination of linear combinations of full grid basis functions, where the coefficients
 * of the outer linear combination are given by the combination grid).
 */
class OperationEvalCombinationGrid {
 public:
  /**
   * Constructor.
   *
   * @param grid  combination grid
   */
  explicit OperationEvalCombinationGrid(const CombinationGrid& grid);

  /**
   * Evaluate a combination grid function at one point.
   *
   * @param surpluses   coefficients for the basis functions (may be nodal/hierarchical),
   *                    every vector corresponds to one full grid (the order of DataVector
   *                    entries is given by IndexVectorRange)
   * @param point       point at which to evaluate the combination grid function
   * @return value of the combination grid function at the given point
   */
  double eval(const std::vector<base::DataVector>& surpluses,
      const base::DataVector& point);

  /**
   * Evaluate a combination grid function at multiple points.
   *
   * @param[in] surpluses   coefficients for the basis functions (may be nodal/hierarchical),
   *                        every vector corresponds to one full grid (the order of DataVector
   *                        entries is given by IndexVectorRange)
   * @param[in] points      points at which to evaluate the combination grid function
   *                        (every row corresponds to one point)
   * @param[out] result     values of the combination grid function at the given points
   */
  void multiEval(const std::vector<base::DataVector>& surpluses,
      const base::DataMatrix& points, base::DataVector& result);

  /**
   * @return combination grid
   */
  const CombinationGrid& getGrid() const;

  /**
   * @param grid  combination grid
   */
  void setGrid(const CombinationGrid& grid);

 protected:
  /// combination grid
  CombinationGrid grid;
};

}  // namespace combigrid
}  // namespace sgpp
