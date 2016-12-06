// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/globaldef.hpp>
#include <sgpp/base/operation/hash/OperationHierarchisationNaturalBsplineBoundary.hpp>
#include <sgpp/base/operation/hash/common/algorithm_bfs/HierarchisationNaturalBsplineBoundary.hpp>
#include <sgpp/base/operation/hash/common/algorithm_bfs/DehierarchisationNaturalBsplineBoundary.hpp>
#include <sgpp/base/algorithm/BreadthFirstSearch.hpp>

namespace sgpp {
namespace base {

OperationHierarchisationNaturalBsplineBoundary::
OperationHierarchisationNaturalBsplineBoundary(NaturalBsplineBoundaryGrid* grid) :
  grid(grid) {
}

OperationHierarchisationNaturalBsplineBoundary::
~OperationHierarchisationNaturalBsplineBoundary() {
}

void OperationHierarchisationNaturalBsplineBoundary::doHierarchisation(
  DataVector& node_values) {
  HierarchisationNaturalBsplineBoundary func(grid);
  BreadthFirstSearch<HierarchisationNaturalBsplineBoundary>
  bfs(func, grid->getStorage());
  bfs.execute(node_values, node_values);
}

void OperationHierarchisationNaturalBsplineBoundary::doDehierarchisation(
  DataVector& alpha) {
  DehierarchisationNaturalBsplineBoundary func(grid);
  BreadthFirstSearch<DehierarchisationNaturalBsplineBoundary>
  bfs(func, grid->getStorage());
  DataVector alphaCopy(alpha);
  bfs.execute(alphaCopy, alpha);
}

void OperationHierarchisationNaturalBsplineBoundary::doHierarchisation(
  DataMatrix& node_values) {
  HierarchisationNaturalBsplineBoundary func(grid);
  BreadthFirstSearch<HierarchisationNaturalBsplineBoundary>
  bfs(func, grid->getStorage());
  bfs.execute(node_values, node_values);
}

void OperationHierarchisationNaturalBsplineBoundary::doDehierarchisation(
  DataMatrix& alpha) {
  DehierarchisationNaturalBsplineBoundary func(grid);
  BreadthFirstSearch<DehierarchisationNaturalBsplineBoundary>
  bfs(func, grid->getStorage());
  DataMatrix alphaCopy(alpha);
  bfs.execute(alphaCopy, alpha);
}

}  // namespace base
}  // namespace sgpp
