// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/globaldef.hpp>
#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/combigrid/CombinationGrid.hpp>
#include <sgpp/combigrid/FullGrid.hpp>
#include <sgpp/combigrid/OperationEvalFullGrid.hpp>

#include <vector>

namespace sgpp {
namespace combigrid {

class OperationEvalCombinationGrid {
 public:
  explicit OperationEvalCombinationGrid(const CombinationGrid& grid) : grid(grid) {
  }

  virtual ~OperationEvalCombinationGrid() {
  }

  virtual double eval(const std::vector<base::DataVector>& surpluses,
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

  virtual void eval(const std::vector<base::DataVector>& surpluses,
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

 protected:
  CombinationGrid grid;
};

}  // namespace combigrid
}  // namespace sgpp
