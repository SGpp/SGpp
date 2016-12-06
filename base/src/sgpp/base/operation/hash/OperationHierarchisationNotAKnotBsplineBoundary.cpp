// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/globaldef.hpp>
#include <sgpp/base/operation/hash/OperationHierarchisationNotAKnotBsplineBoundary.hpp>
#include <sgpp/base/operation/hash/common/algorithm_bfs/HierarchisationNotAKnotBsplineBoundary.hpp>
#include <sgpp/base/operation/hash/common/algorithm_bfs/DehierarchisationNotAKnotBsplineBoundary.hpp>
#include <sgpp/base/algorithm/BreadthFirstSearch.hpp>

namespace sgpp {
namespace base {

OperationHierarchisationNotAKnotBsplineBoundary::
OperationHierarchisationNotAKnotBsplineBoundary(NotAKnotBsplineBoundaryGrid* grid) :
  grid(grid) {
}

OperationHierarchisationNotAKnotBsplineBoundary::
~OperationHierarchisationNotAKnotBsplineBoundary() {
}

void OperationHierarchisationNotAKnotBsplineBoundary::doHierarchisation(
  DataVector& node_values) {
  HierarchisationNotAKnotBsplineBoundary func(grid);
  BreadthFirstSearch<HierarchisationNotAKnotBsplineBoundary>
  bfs(func, grid->getStorage());
  bfs.execute(node_values, node_values);
}

void OperationHierarchisationNotAKnotBsplineBoundary::doDehierarchisation(
  DataVector& alpha) {
  DehierarchisationNotAKnotBsplineBoundary func(grid);
  BreadthFirstSearch<DehierarchisationNotAKnotBsplineBoundary>
  bfs(func, grid->getStorage());
  DataVector alphaCopy(alpha);
  bfs.execute(alphaCopy, alpha);
}

void OperationHierarchisationNotAKnotBsplineBoundary::doHierarchisation(
  DataMatrix& node_values) {
  HierarchisationNotAKnotBsplineBoundary func(grid);
  BreadthFirstSearch<HierarchisationNotAKnotBsplineBoundary>
  bfs(func, grid->getStorage());
  bfs.execute(node_values, node_values);
}

void OperationHierarchisationNotAKnotBsplineBoundary::doDehierarchisation(
  DataMatrix& alpha) {
  DehierarchisationNotAKnotBsplineBoundary func(grid);
  BreadthFirstSearch<DehierarchisationNotAKnotBsplineBoundary>
  bfs(func, grid->getStorage());
  DataMatrix alphaCopy(alpha);
  bfs.execute(alphaCopy, alpha);
}

}  // namespace base
}  // namespace sgpp
