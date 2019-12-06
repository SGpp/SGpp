// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/globaldef.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/combigrid/grid/FullGrid.hpp>
#include <sgpp/combigrid/operation/OperationPole.hpp>

#include <memory>
#include <vector>

namespace sgpp {
namespace combigrid {

/**
 * Operation for applying 1D OperationPole operators on all poles of a full grid in all dimensions
 * via the unidirectional principle (UP).
 */
class OperationUPFullGrid {
 public:
  /**
   * Constructor.
   *
   * @param grid            full grid
   * @param operationPole   vector of unique_ptr to 1D pole operators
   *                        (do not destruct before this object)
   */
  OperationUPFullGrid(const FullGrid& grid,
      const std::vector<std::unique_ptr<OperationPole>>& operationPole);

  /**
   * Constructor.
   *
   * @param grid            full grid
   * @param operationPole   vector of pointers to 1D pole operators
   *                        (do not delete before this object)
   */
  OperationUPFullGrid(const FullGrid& grid, const std::vector<OperationPole*>& operationPole);

  /**
   * Constructor for the special case where the same OperationPole should be used for all
   * dimensions.
   *
   * @param grid            full grid
   * @param operationPole   1D pole operator (do not destruct before this object)
   */
  OperationUPFullGrid(const FullGrid& grid, OperationPole& operationPole);

  /**
   * Apply the unidirectional principle in-place.
   *
   * @param[in,out] values  data vector, same size as the number of grid points of the full grid
   *                        (the order is given by IndexVectorRange)
   */
  void apply(base::DataVector& values);

  /**
   * @return full grid
   */
  const FullGrid& getGrid() const;

  /**
   * @param grid  full grid
   */
  void setGrid(const FullGrid& grid);

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
  /// full grid
  FullGrid grid;
  /// vector of pointers to 1D pole operators (do not delete before this object)
  std::vector<OperationPole*> operationPole;
};

}  // namespace combigrid
}  // namespace sgpp
