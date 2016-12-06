// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/globaldef.hpp>
#include <sgpp/base/operation/hash/OperationHierarchisationModNotAKnotBspline.hpp>
#include <sgpp/base/operation/hash/common/algorithm_bfs/HierarchisationModNotAKnotBspline.hpp>
#include <sgpp/base/operation/hash/common/algorithm_bfs/DehierarchisationModNotAKnotBspline.hpp>
#include <sgpp/base/algorithm/BreadthFirstSearch.hpp>

namespace sgpp {
namespace base {

OperationHierarchisationModNotAKnotBspline::
OperationHierarchisationModNotAKnotBspline(ModNotAKnotBsplineGrid* grid) :
  grid(grid) {
}

OperationHierarchisationModNotAKnotBspline::
~OperationHierarchisationModNotAKnotBspline() {
}

void OperationHierarchisationModNotAKnotBspline::doHierarchisation(
  DataVector& node_values) {
  HierarchisationModNotAKnotBspline func(grid);
  BreadthFirstSearch<HierarchisationModNotAKnotBspline>
  bfs(func, grid->getStorage());
  bfs.execute(node_values, node_values);
}

void OperationHierarchisationModNotAKnotBspline::doDehierarchisation(
  DataVector& alpha) {
  DehierarchisationModNotAKnotBspline func(grid);
  BreadthFirstSearch<DehierarchisationModNotAKnotBspline>
  bfs(func, grid->getStorage());
  DataVector alphaCopy(alpha);
  bfs.execute(alphaCopy, alpha);
}

void OperationHierarchisationModNotAKnotBspline::doHierarchisation(
  DataMatrix& node_values) {
  HierarchisationModNotAKnotBspline func(grid);
  BreadthFirstSearch<HierarchisationModNotAKnotBspline>
  bfs(func, grid->getStorage());
  bfs.execute(node_values, node_values);
}

void OperationHierarchisationModNotAKnotBspline::doDehierarchisation(
  DataMatrix& alpha) {
  DehierarchisationModNotAKnotBspline func(grid);
  BreadthFirstSearch<DehierarchisationModNotAKnotBspline>
  bfs(func, grid->getStorage());
  DataMatrix alphaCopy(alpha);
  bfs.execute(alphaCopy, alpha);
}

}  // namespace base
}  // namespace sgpp
