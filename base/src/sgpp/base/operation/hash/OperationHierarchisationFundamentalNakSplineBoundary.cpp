// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/globaldef.hpp>
#include <sgpp/base/operation/hash/OperationHierarchisationFundamentalNakSplineBoundary.hpp>
#include <sgpp/base/operation/hash/common/algorithm_bfs/HierarchisationFundamentalNakSplineBoundary.hpp>
#include <sgpp/base/operation/hash/common/algorithm_bfs/DehierarchisationFundamentalNakSplineBoundary.hpp>
#include <sgpp/base/algorithm/BreadthFirstSearch.hpp>

namespace sgpp {
namespace base {

OperationHierarchisationFundamentalNakSplineBoundary::
OperationHierarchisationFundamentalNakSplineBoundary(
    FundamentalNakSplineBoundaryGrid* grid) :
  grid(grid) {
}

OperationHierarchisationFundamentalNakSplineBoundary::
~OperationHierarchisationFundamentalNakSplineBoundary() {
}

void OperationHierarchisationFundamentalNakSplineBoundary::doHierarchisation(
  DataVector& node_values) {
  HierarchisationFundamentalNakSplineBoundary func(grid);
  BreadthFirstSearch<HierarchisationFundamentalNakSplineBoundary>
  bfs(func, grid->getStorage());
  bfs.execute(node_values, node_values);
}

void OperationHierarchisationFundamentalNakSplineBoundary::doDehierarchisation(
  DataVector& alpha) {
  DehierarchisationFundamentalNakSplineBoundary func(grid);
  BreadthFirstSearch<DehierarchisationFundamentalNakSplineBoundary>
  bfs(func, grid->getStorage());
  DataVector alphaCopy(alpha);
  bfs.execute(alphaCopy, alpha);
}

void OperationHierarchisationFundamentalNakSplineBoundary::doHierarchisation(
  DataMatrix& node_values) {
  HierarchisationFundamentalNakSplineBoundary func(grid);
  BreadthFirstSearch<HierarchisationFundamentalNakSplineBoundary>
  bfs(func, grid->getStorage());
  bfs.execute(node_values, node_values);
}

void OperationHierarchisationFundamentalNakSplineBoundary::doDehierarchisation(
  DataMatrix& alpha) {
  DehierarchisationFundamentalNakSplineBoundary func(grid);
  BreadthFirstSearch<DehierarchisationFundamentalNakSplineBoundary>
  bfs(func, grid->getStorage());
  DataMatrix alphaCopy(alpha);
  bfs.execute(alphaCopy, alpha);
}

}  // namespace base
}  // namespace sgpp
