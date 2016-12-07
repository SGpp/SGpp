// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/globaldef.hpp>
#include <sgpp/base/operation/hash/OperationHierarchisationModLagrangeNotAKnotSpline.hpp>
#include <sgpp/base/operation/hash/common/algorithm_bfs/HierarchisationModLagrangeNotAKnotSpline.hpp>
#include <sgpp/base/operation/hash/common/algorithm_bfs/DehierarchisationModLagrangeNotAKnotSpline.hpp>
#include <sgpp/base/algorithm/BreadthFirstSearch.hpp>

namespace sgpp {
namespace base {

OperationHierarchisationModLagrangeNotAKnotSpline::
OperationHierarchisationModLagrangeNotAKnotSpline(ModLagrangeNotAKnotSplineGrid* grid) :
  grid(grid) {
}

OperationHierarchisationModLagrangeNotAKnotSpline::
~OperationHierarchisationModLagrangeNotAKnotSpline() {
}

void OperationHierarchisationModLagrangeNotAKnotSpline::doHierarchisation(
  DataVector& node_values) {
  HierarchisationModLagrangeNotAKnotSpline func(grid);
  BreadthFirstSearch<HierarchisationModLagrangeNotAKnotSpline>
  bfs(func, grid->getStorage());
  bfs.execute(node_values, node_values);
}

void OperationHierarchisationModLagrangeNotAKnotSpline::doDehierarchisation(
  DataVector& alpha) {
  DehierarchisationModLagrangeNotAKnotSpline func(grid);
  BreadthFirstSearch<DehierarchisationModLagrangeNotAKnotSpline>
  bfs(func, grid->getStorage());
  DataVector alphaCopy(alpha);
  bfs.execute(alphaCopy, alpha);
}

void OperationHierarchisationModLagrangeNotAKnotSpline::doHierarchisation(
  DataMatrix& node_values) {
  HierarchisationModLagrangeNotAKnotSpline func(grid);
  BreadthFirstSearch<HierarchisationModLagrangeNotAKnotSpline>
  bfs(func, grid->getStorage());
  bfs.execute(node_values, node_values);
}

void OperationHierarchisationModLagrangeNotAKnotSpline::doDehierarchisation(
  DataMatrix& alpha) {
  DehierarchisationModLagrangeNotAKnotSpline func(grid);
  BreadthFirstSearch<DehierarchisationModLagrangeNotAKnotSpline>
  bfs(func, grid->getStorage());
  DataMatrix alphaCopy(alpha);
  bfs.execute(alphaCopy, alpha);
}

}  // namespace base
}  // namespace sgpp
