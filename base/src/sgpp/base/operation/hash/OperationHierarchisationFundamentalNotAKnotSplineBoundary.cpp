// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/globaldef.hpp>
#include <sgpp/base/operation/hash/OperationHierarchisationFundamentalNotAKnotSplineBoundary.hpp>
#include <sgpp/base/operation/hash/common/algorithm_bfs/HierarchisationFundamentalNotAKnotSplineBoundary.hpp>
#include <sgpp/base/operation/hash/common/algorithm_bfs/DehierarchisationFundamentalNotAKnotSplineBoundary.hpp>
#include <sgpp/base/algorithm/BreadthFirstSearch.hpp>

namespace sgpp {
namespace base {

OperationHierarchisationFundamentalNotAKnotSplineBoundary::
OperationHierarchisationFundamentalNotAKnotSplineBoundary(
    FundamentalNotAKnotSplineBoundaryGrid* grid) :
  grid(grid) {
}

OperationHierarchisationFundamentalNotAKnotSplineBoundary::
~OperationHierarchisationFundamentalNotAKnotSplineBoundary() {
}

void OperationHierarchisationFundamentalNotAKnotSplineBoundary::doHierarchisation(
  DataVector& node_values) {
  HierarchisationFundamentalNotAKnotSplineBoundary func(grid);
  BreadthFirstSearch<HierarchisationFundamentalNotAKnotSplineBoundary>
  bfs(func, grid->getStorage());
  bfs.execute(node_values, node_values);
}

void OperationHierarchisationFundamentalNotAKnotSplineBoundary::doDehierarchisation(
  DataVector& alpha) {
  DehierarchisationFundamentalNotAKnotSplineBoundary func(grid);
  BreadthFirstSearch<DehierarchisationFundamentalNotAKnotSplineBoundary>
  bfs(func, grid->getStorage());
  DataVector alphaCopy(alpha);
  bfs.execute(alphaCopy, alpha);
}

void OperationHierarchisationFundamentalNotAKnotSplineBoundary::doHierarchisation(
  DataMatrix& node_values) {
  HierarchisationFundamentalNotAKnotSplineBoundary func(grid);
  BreadthFirstSearch<HierarchisationFundamentalNotAKnotSplineBoundary>
  bfs(func, grid->getStorage());
  bfs.execute(node_values, node_values);
}

void OperationHierarchisationFundamentalNotAKnotSplineBoundary::doDehierarchisation(
  DataMatrix& alpha) {
  DehierarchisationFundamentalNotAKnotSplineBoundary func(grid);
  BreadthFirstSearch<DehierarchisationFundamentalNotAKnotSplineBoundary>
  bfs(func, grid->getStorage());
  DataMatrix alphaCopy(alpha);
  bfs.execute(alphaCopy, alpha);
}

}  // namespace base
}  // namespace sgpp
