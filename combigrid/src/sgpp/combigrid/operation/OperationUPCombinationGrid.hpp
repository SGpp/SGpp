// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/globaldef.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/combigrid/grid/CombinationGrid.hpp>
#include <sgpp/combigrid/grid/FullGrid.hpp>
#include <sgpp/combigrid/operation/OperationPole.hpp>
#include <sgpp/combigrid/operation/OperationUPFullGrid.hpp>

#include <vector>

namespace sgpp {
namespace combigrid {

class OperationUPCombinationGrid {
 public:
  OperationUPCombinationGrid(const CombinationGrid& grid,
      const std::vector<std::unique_ptr<OperationPole>>& operationPole) :
      grid(grid), operationPole() {
    for (const std::unique_ptr<OperationPole>& operationPole1d : operationPole) {
      this->operationPole.push_back(operationPole1d.get());
    }
  }

  OperationUPCombinationGrid(const CombinationGrid& grid,
      const std::vector<OperationPole*> operationPole) :
      grid(grid), operationPole(operationPole) {
  }

  OperationUPCombinationGrid(const CombinationGrid& grid, OperationPole& operationPole) :
      grid(grid), operationPole(grid.getDimension(), &operationPole) {
  }

  virtual ~OperationUPCombinationGrid() {
  }

  virtual void apply(std::vector<base::DataVector>& values) {
    const std::vector<FullGrid>& fullGrids = grid.getFullGrids();

    if (fullGrids.empty()) {
      return;
    }

    OperationUPFullGrid operationUPFullGrid(fullGrids[0], operationPole);

    for (size_t i = 0; i < values.size(); i++) {
      operationUPFullGrid.setGrid(fullGrids[i]);
      operationUPFullGrid.apply(values[i]);
    }
  }

 protected:
  CombinationGrid grid;
  std::vector<OperationPole*> operationPole;
};

}  // namespace combigrid
}  // namespace sgpp
