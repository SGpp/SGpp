// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/globaldef.hpp>

#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/function/scalar/ScalarFunction.hpp>
#include <sgpp/base/grid/Grid.hpp>

#include <cstddef>
namespace sgpp {
namespace optimization {

/**
 * Abstract base class for iterative grid generation methods.
 */
class IterativeGridGenerator {
 public:
  /**
   * Constructor.
   * Do not destruct the grid before this object!
   *
   * @param f     objective function
   * @param grid  grid (should be empty)
   * @param N     maximal number of grid points
   */
  IterativeGridGenerator(base::ScalarFunction& f, base::Grid& grid, size_t N);

  /**
   * Destructor.
   */
  virtual ~IterativeGridGenerator();

  /**
   * Pure virtual method for iterative grid generation.
   *
   * @return true on success, otherwise false
   */
  virtual bool generate() = 0;

  /**
   * @return underlying grid
   */
  base::Grid& getGrid() const;

  /**
   * @return vector of function values at the grid points
   */
  const base::DataVector& getFunctionValues() const;

  /**
   * Print a grid (grid points and function values).
   *
   * @param gridGen       grid to be printed
   */
  void printIterativeGridGenerator() const;

 protected:
  /// objective function
  base::ScalarFunction& f;
  /// underlying grid
  base::Grid& grid;
  /// maximal number of grid points
  size_t N;
  /// vector of function values at the grid points
  base::DataVector functionValues;

  /**
   * Removes grid points with indices
   * [oldGridSize, oldGridSize + 1, ..., grid.getSize() - 1]
   * from the grid.
   *
   * @param oldGridSize   number of grid points after removal
   */
  void undoRefinement(size_t oldGridSize);

  /**
   * Evaluates the objective function at grid points with indices
   * [oldGridSize, oldGridSize + 1, ..., grid.getSize() - 1]
   * and saves values in functionValues.
   *
   * @param oldGridSize   number of grid points already evaluated
   */
  void evalFunction(size_t oldGridSize = 0);
};
}  // namespace optimization
}  // namespace sgpp
