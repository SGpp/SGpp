// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/globaldef.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/combigrid/grid/CombinationGrid.hpp>
#include <sgpp/combigrid/operation/OperationPole.hpp>

#include <vector>

namespace sgpp {
namespace combigrid {

/**
 * Operation for applying 1D OperationPole operators on all poles of all full grids of some
 * combination grid, in all dimensions via the unidirectional principle (UP).
 */
class OperationUPCombinationGrid {
 public:
  /**
   * Constructor.
   *
   * @param grid            combination grid
   * @param operationPole   vector of unique_ptr to 1D pole operators
   *                        (do not destruct before this object)
   */
  OperationUPCombinationGrid(const CombinationGrid& grid,
      const std::vector<std::unique_ptr<OperationPole>>& operationPole);

  /**
   * Constructor.
   *
   * @param grid            combination grid
   * @param operationPole   vector of pointers to 1D pole operators
   *                        (do not delete before this object)
   */
  OperationUPCombinationGrid(const CombinationGrid& grid,
      const std::vector<OperationPole*> operationPole);

  /**
   * Constructor for the special case where the same OperationPole should be used for all
   * dimensions.
   *
   * @param grid            combination grid
   * @param operationPole   1D pole operator (do not destruct before this object)
   */
  OperationUPCombinationGrid(const CombinationGrid& grid, OperationPole& operationPole);

  /**
   * Apply the unidirectional principle in-place.
   *
   * @param[in,out] values  vector of vectors with values on the full grids, every vector
   *                        corresponds to one full grid of the combination grid, every vector
   *                        has the same size as the number of grid points of the respective
   *                        full grid (the order of DataVector entries is given by
   *                        IndexVectorRange)
   */
  void apply(std::vector<base::DataVector>& values);

  /**
   * @return combination grid
   */
  const CombinationGrid& getGrid() const;

  /**
   * @param grid  combination grid
   */
  void setGrid(const CombinationGrid& grid);

  /**
   * @return vector of pointers to 1D pole operators (do not delete before this object)
   */
  const std::vector<OperationPole*>& getOperationPole() const;

  /**
   * @param operationPole   vector of pointers to 1D pole operators
   *                        (do not delete before this object)
   */
  void setOperationPole(const std::vector<OperationPole*>& operationPole);

 protected:
  /// combination grid
  CombinationGrid grid;
  /// vector of pointers to 1D pole operators (do not delete before this object)
  std::vector<OperationPole*> operationPole;
};

}  // namespace combigrid
}  // namespace sgpp
