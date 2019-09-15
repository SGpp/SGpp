// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/globaldef.hpp>
#include <sgpp/combigrid/grid/FullGrid.hpp>
#include <sgpp/combigrid/operation/OperationEvalCombinationGrid.hpp>
#include <sgpp/combigrid/operation/OperationEvalFullGrid.hpp>

#include <vector>

namespace sgpp {
namespace combigrid {

OperationEvalCombinationGrid::OperationEvalCombinationGrid(const CombinationGrid& grid) :
    grid(grid) {
}

double OperationEvalCombinationGrid::eval(const std::vector<base::DataVector>& surpluses,
    const base::DataVector& point) {
  const std::vector<FullGrid>& fullGrids = grid.getFullGrids();
  base::DataVector values(fullGrids.size());
  OperationEvalFullGrid operationEvalFullGrid;

  for (size_t i = 0; i < fullGrids.size(); i++) {
    operationEvalFullGrid.setGrid(fullGrids[i]);
    values[i] = operationEvalFullGrid.eval(surpluses[i], point);
  }

  return grid.combineValues(values);
}

void OperationEvalCombinationGrid::eval(const std::vector<base::DataVector>& surpluses,
    const base::DataMatrix& points, base::DataVector& result) {
  const std::vector<FullGrid>& fullGrids = grid.getFullGrids();
  base::DataMatrix values(points.getNrows(), fullGrids.size());
  base::DataVector curValues(points.getNrows());
  OperationEvalFullGrid operationEvalFullGrid;

  for (size_t i = 0; i < fullGrids.size(); i++) {
    operationEvalFullGrid.setGrid(fullGrids[i]);
    operationEvalFullGrid.eval(surpluses[i], points, curValues);
    values.setColumn(i, curValues);
  }

  grid.combineValues(values, result);
}

}  // namespace combigrid
}  // namespace sgpp
