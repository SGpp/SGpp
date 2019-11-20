// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/globaldef.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/combigrid/LevelIndexTypes.hpp>

namespace sgpp {
namespace combigrid {

/**
 * Operation on a pole of a full grid. A pole is a one-dimensional sub-grid (1D entries of
 * the index vector/coordinates are fixed in all dimensions but one).
 */
class OperationPole {
 public:
  /**
   * Default constructor.
   */
  OperationPole() {
  }

  /**
   * Virtual destructor.
   */
  virtual ~OperationPole() {
  }

  /**
   * Apply the operator on data. This operates in-place on the whole data vector (which is assumed
   * to hold data for all full grid points), so this should only be called for 1D grids.
   *
   * @param[in,out] values    data vector for all full grid points
   *                          (the order is given by IndexVectorRange)
   * @param[in] level         level of the full grid
   * @param[in] hasBoundary   whether the full grid has points on the boundary
   */
  virtual void apply(base::DataVector& values, level_t level, bool hasBoundary = true) {
    apply(values, 0, 1, values.size(), level, hasBoundary);
  }

  /**
   * Apply the operator on data. This operates in-place on a subset of the given data vector
   * (which is assumed to hold data for all full grid points), where the specified subset
   * exactly corresponds to the points of the pole.
   *
   * @param[in,out] values    data vector for all full grid points
   *                          (the order is given by IndexVectorRange)
   * @param[in] start         sequence number of the first grid point of the pole
   * @param[in] step          difference of sequence numbers of two subsequent grid points of
   *                          the pole
   * @param[in] count         number of grid points of the pole
   * @param[in] level         level of the full grid
   * @param[in] hasBoundary   whether the full grid has points on the boundary
   */
  virtual void apply(base::DataVector& values, size_t start, size_t step, size_t count,
      level_t level, bool hasBoundary = true) = 0;
};

}  // namespace combigrid
}  // namespace sgpp
