// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/globaldef.hpp>
#include <sgpp/base/operation/hash/OperationHierarchisationLagrangeSplineBoundary.hpp>
#include <sgpp/base/operation/hash/common/algorithm_bfs/HierarchisationLagrangeSplineBoundary.hpp>
#include <sgpp/base/operation/hash/common/algorithm_bfs/DehierarchisationLagrangeSplineBoundary.hpp>
#include <sgpp/base/algorithm/BreadthFirstSearch.hpp>

namespace sgpp {
namespace base {

OperationHierarchisationLagrangeSplineBoundary::
OperationHierarchisationLagrangeSplineBoundary(LagrangeSplineBoundaryGrid* grid) :
  grid(grid) {
}

OperationHierarchisationLagrangeSplineBoundary::
~OperationHierarchisationLagrangeSplineBoundary() {
}

void OperationHierarchisationLagrangeSplineBoundary::doHierarchisation(
  DataVector& node_values) {
  HierarchisationLagrangeSplineBoundary func(grid);
  BreadthFirstSearch<HierarchisationLagrangeSplineBoundary>
  bfs(func, grid->getStorage());
  bfs.execute(node_values, node_values);
}

void OperationHierarchisationLagrangeSplineBoundary::doDehierarchisation(
  DataVector& alpha) {
  DehierarchisationLagrangeSplineBoundary func(grid);
  BreadthFirstSearch<DehierarchisationLagrangeSplineBoundary>
  bfs(func, grid->getStorage());
  DataVector alphaCopy(alpha);
  bfs.execute(alphaCopy, alpha);
}

void OperationHierarchisationLagrangeSplineBoundary::doHierarchisation(
  DataMatrix& node_values) {
  HierarchisationLagrangeSplineBoundary func(grid);
  BreadthFirstSearch<HierarchisationLagrangeSplineBoundary>
  bfs(func, grid->getStorage());
  bfs.execute(node_values, node_values);
}

void OperationHierarchisationLagrangeSplineBoundary::doDehierarchisation(
  DataMatrix& alpha) {
  DehierarchisationLagrangeSplineBoundary func(grid);
  BreadthFirstSearch<DehierarchisationLagrangeSplineBoundary>
  bfs(func, grid->getStorage());
  DataMatrix alphaCopy(alpha);
  bfs.execute(alphaCopy, alpha);
}

}  // namespace base
}  // namespace sgpp
