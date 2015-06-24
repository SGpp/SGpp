// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/globaldef.hpp>
#include <sgpp/base/operation/hash/OperationHierarchisationFundamentalSpline.hpp>
#include <sgpp/base/operation/hash/common/algorithm_bfs/HierarchisationFundamentalSpline.hpp>
#include <sgpp/base/operation/hash/common/algorithm_bfs/DehierarchisationFundamentalSpline.hpp>
#include <sgpp/base/algorithm/BreadthFirstSearch.hpp>

namespace SGPP {
  namespace base {

    OperationHierarchisationFundamentalSpline::
    OperationHierarchisationFundamentalSpline(FundamentalSplineGrid* grid) :
      grid(grid) {
    }

    OperationHierarchisationFundamentalSpline::
    ~OperationHierarchisationFundamentalSpline() {
    }

    void OperationHierarchisationFundamentalSpline::doHierarchisation(
      DataVector& node_values) {
      HierarchisationFundamentalSpline func(grid);
      BreadthFirstSearch<HierarchisationFundamentalSpline>
      bfs(func, grid->getStorage());
      bfs.execute(node_values, node_values);
    }

    void OperationHierarchisationFundamentalSpline::doDehierarchisation(
      DataVector& alpha) {
      DehierarchisationFundamentalSpline func(grid);
      BreadthFirstSearch<DehierarchisationFundamentalSpline>
      bfs(func, grid->getStorage());
      DataVector alphaCopy(alpha);
      bfs.execute(alphaCopy, alpha);
    }
  }
}
