// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/globaldef.hpp>
#include <sgpp/base/operation/hash/OperationHierarchisationFundamentalSplineBoundary.hpp>
#include <sgpp/base/operation/hash/common/algorithm_bfs/HierarchisationFundamentalSplineBoundary.hpp>
#include <sgpp/base/operation/hash/common/algorithm_bfs/DehierarchisationFundamentalSplineBoundary.hpp>
#include <sgpp/base/algorithm/BreadthFirstSearch.hpp>

namespace sgpp {
namespace base {

OperationHierarchisationFundamentalSplineBoundary::
OperationHierarchisationFundamentalSplineBoundary(FundamentalSplineBoundaryGrid* grid) :
  grid(grid) {
}

OperationHierarchisationFundamentalSplineBoundary::
~OperationHierarchisationFundamentalSplineBoundary() {
}

void OperationHierarchisationFundamentalSplineBoundary::doHierarchisation(
  DataVector& node_values) {
  HierarchisationFundamentalSplineBoundary func(grid);
  BreadthFirstSearch<HierarchisationFundamentalSplineBoundary>
  bfs(func, grid->getStorage());
  bfs.execute(node_values, node_values);
}

void OperationHierarchisationFundamentalSplineBoundary::doDehierarchisation(
  DataVector& alpha) {
  DehierarchisationFundamentalSplineBoundary func(grid);
  BreadthFirstSearch<DehierarchisationFundamentalSplineBoundary>
  bfs(func, grid->getStorage());
  DataVector alphaCopy(alpha);
  bfs.execute(alphaCopy, alpha);
}

void OperationHierarchisationFundamentalSplineBoundary::doHierarchisation(
  DataMatrix& node_values) {
  HierarchisationFundamentalSplineBoundary func(grid);
  BreadthFirstSearch<HierarchisationFundamentalSplineBoundary>
  bfs(func, grid->getStorage());
  bfs.execute(node_values, node_values);
}

void OperationHierarchisationFundamentalSplineBoundary::doDehierarchisation(
  DataMatrix& alpha) {
  DehierarchisationFundamentalSplineBoundary func(grid);
  BreadthFirstSearch<DehierarchisationFundamentalSplineBoundary>
  bfs(func, grid->getStorage());
  DataMatrix alphaCopy(alpha);
  bfs.execute(alphaCopy, alpha);
}

}  // namespace base
}  // namespace sgpp
