// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/globaldef.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/combigrid/LevelIndexTypes.hpp>
#include <sgpp/combigrid/operation/OperationPole.hpp>

namespace sgpp {
namespace combigrid {

/**
 * Operation of dehierarchising values on a 1D pole of a full grid with linear basis functions.
 */
class OperationPoleDehierarchisationLinear : public OperationPole {
 public:
  /**
   * Default constructor.
   */
  OperationPoleDehierarchisationLinear();

  /**
   * Virtual destructor.
   */
  ~OperationPoleDehierarchisationLinear() override;

  /**
   * Apply the operator on data.
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
  void apply(base::DataVector& values, size_t start, size_t step, size_t count,
      level_t level, bool hasBoundary = true) override;
};

}  // namespace combigrid
}  // namespace sgpp
