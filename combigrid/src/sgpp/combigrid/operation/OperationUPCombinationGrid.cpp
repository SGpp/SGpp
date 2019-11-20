// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/globaldef.hpp>
#include <sgpp/combigrid/grid/FullGrid.hpp>
#include <sgpp/combigrid/operation/OperationUPCombinationGrid.hpp>
#include <sgpp/combigrid/operation/OperationUPFullGrid.hpp>

#include <vector>

namespace sgpp {
namespace combigrid {

OperationUPCombinationGrid::OperationUPCombinationGrid(const CombinationGrid& grid,
    const std::vector<std::unique_ptr<OperationPole>>& operationPole) :
    grid(grid), operationPole() {
  for (const std::unique_ptr<OperationPole>& operationPole1d : operationPole) {
    this->operationPole.push_back(operationPole1d.get());
  }
}

OperationUPCombinationGrid::OperationUPCombinationGrid(const CombinationGrid& grid,
    const std::vector<OperationPole*> operationPole) :
    grid(grid), operationPole(operationPole) {
}

OperationUPCombinationGrid::OperationUPCombinationGrid(const CombinationGrid& grid,
    OperationPole& operationPole) : grid(grid), operationPole(grid.getDimension(), &operationPole) {
}

void OperationUPCombinationGrid::apply(std::vector<base::DataVector>& values) {
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

const CombinationGrid& OperationUPCombinationGrid::getGrid() const {
  return grid;
}

void OperationUPCombinationGrid::setGrid(const CombinationGrid& grid) {
  this->grid = grid;
}

const std::vector<OperationPole*>& OperationUPCombinationGrid::getOperationPole() const {
  return operationPole;
}

void OperationUPCombinationGrid::setOperationPole(
    const std::vector<OperationPole*>& operationPole) {
  this->operationPole = operationPole;
}

}  // namespace combigrid
}  // namespace sgpp
