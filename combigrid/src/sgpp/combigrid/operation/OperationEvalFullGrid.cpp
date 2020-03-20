// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/globaldef.hpp>
#include <sgpp/combigrid/LevelIndexTypes.hpp>
#include <sgpp/combigrid/operation/OperationEvalFullGrid.hpp>
#include <sgpp/combigrid/tools/IndexVectorRange.hpp>

namespace sgpp {
namespace combigrid {

OperationEvalFullGrid::OperationEvalFullGrid() : grid() {
}

OperationEvalFullGrid::OperationEvalFullGrid(const FullGrid& grid) : grid(grid) {
}

OperationEvalFullGrid::~OperationEvalFullGrid() {
}

double OperationEvalFullGrid::eval(const base::DataVector& surpluses,
    const base::DataVector& point) {
  const LevelVector& level = grid.getLevel();
  const HeterogeneousBasis& basis = grid.getBasis();
  size_t i = 0;
  double result = 0.0;

  for (const IndexVector& index : IndexVectorRange(grid)) {
    result += surpluses[i] * basis.eval(level, index, point);
    i++;
  }

  return result;
}

void OperationEvalFullGrid::multiEval(const base::DataVector& surpluses,
    const base::DataMatrix& points, base::DataVector& result) {
  const LevelVector& level = grid.getLevel();
  const HeterogeneousBasis& basis = grid.getBasis();
  const size_t n = points.getNrows();
  base::DataVector point(points.getNcols());
  size_t i = 0;
  result.resize(n);
  result.setAll(0.0);

  for (const IndexVector& index : IndexVectorRange(grid)) {
    for (size_t j = 0; j < n; j++) {
      points.getRow(j, point);
      result[j] += surpluses[i] * basis.eval(level, index, point);
    }

    i++;
  }
}

const FullGrid& OperationEvalFullGrid::getGrid() const {
  return grid;
}

void OperationEvalFullGrid::setGrid(const FullGrid& grid) {
  this->grid = grid;
}

}  // namespace combigrid
}  // namespace sgpp
