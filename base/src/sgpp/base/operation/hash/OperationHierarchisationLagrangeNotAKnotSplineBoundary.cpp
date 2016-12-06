// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/globaldef.hpp>
#include <sgpp/base/operation/hash/OperationHierarchisationLagrangeNotAKnotSplineBoundary.hpp>
#include <sgpp/base/operation/hash/common/algorithm_bfs/HierarchisationLagrangeNotAKnotSplineBoundary.hpp>
#include <sgpp/base/operation/hash/common/algorithm_bfs/DehierarchisationLagrangeNotAKnotSplineBoundary.hpp>
#include <sgpp/base/algorithm/BreadthFirstSearch.hpp>

namespace sgpp {
namespace base {

OperationHierarchisationLagrangeNotAKnotSplineBoundary::
OperationHierarchisationLagrangeNotAKnotSplineBoundary(LagrangeNotAKnotSplineBoundaryGrid* grid) :
  grid(grid) {
}

OperationHierarchisationLagrangeNotAKnotSplineBoundary::
~OperationHierarchisationLagrangeNotAKnotSplineBoundary() {
}

void OperationHierarchisationLagrangeNotAKnotSplineBoundary::doHierarchisation(
  DataVector& node_values) {
  HierarchisationLagrangeNotAKnotSplineBoundary func(grid);
  BreadthFirstSearch<HierarchisationLagrangeNotAKnotSplineBoundary>
  bfs(func, grid->getStorage());
  bfs.execute(node_values, node_values);
}

void OperationHierarchisationLagrangeNotAKnotSplineBoundary::doDehierarchisation(
  DataVector& alpha) {
  DehierarchisationLagrangeNotAKnotSplineBoundary func(grid);
  BreadthFirstSearch<DehierarchisationLagrangeNotAKnotSplineBoundary>
  bfs(func, grid->getStorage());
  DataVector alphaCopy(alpha);
  bfs.execute(alphaCopy, alpha);
}

void OperationHierarchisationLagrangeNotAKnotSplineBoundary::doHierarchisation(
  DataMatrix& node_values) {
  HierarchisationLagrangeNotAKnotSplineBoundary func(grid);
  BreadthFirstSearch<HierarchisationLagrangeNotAKnotSplineBoundary>
  bfs(func, grid->getStorage());
  bfs.execute(node_values, node_values);
}

void OperationHierarchisationLagrangeNotAKnotSplineBoundary::doDehierarchisation(
  DataMatrix& alpha) {
  DehierarchisationLagrangeNotAKnotSplineBoundary func(grid);
  BreadthFirstSearch<DehierarchisationLagrangeNotAKnotSplineBoundary>
  bfs(func, grid->getStorage());
  DataMatrix alphaCopy(alpha);
  bfs.execute(alphaCopy, alpha);
}

}  // namespace base
}  // namespace sgpp
