// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/globaldef.hpp>
#include <sgpp/base/operation/hash/OperationHierarchisationModFundamentalSpline.hpp>
#include <sgpp/base/operation/hash/common/algorithm_bfs/HierarchisationModFundamentalSpline.hpp>
#include <sgpp/base/operation/hash/common/algorithm_bfs/DehierarchisationModFundamentalSpline.hpp>
#include <sgpp/base/algorithm/BreadthFirstSearch.hpp>

namespace SGPP {
  namespace base {

    OperationHierarchisationModFundamentalSpline::
    OperationHierarchisationModFundamentalSpline(
      ModFundamentalSplineGrid* grid) :
      grid(grid) {
    }

    OperationHierarchisationModFundamentalSpline::
    ~OperationHierarchisationModFundamentalSpline() {
    }

    void OperationHierarchisationModFundamentalSpline::doHierarchisation(
      DataVector& node_values) {
      HierarchisationModFundamentalSpline func(grid);
      BreadthFirstSearch<HierarchisationModFundamentalSpline>
      bfs(func, grid->getStorage());
      bfs.execute(node_values, node_values);
    }

    void OperationHierarchisationModFundamentalSpline::doDehierarchisation(
      DataVector& alpha) {
      DehierarchisationModFundamentalSpline func(grid);
      BreadthFirstSearch<DehierarchisationModFundamentalSpline>
      bfs(func, grid->getStorage());
      DataVector alphaCopy(alpha);
      bfs.execute(alphaCopy, alpha);
    }
  }
}
