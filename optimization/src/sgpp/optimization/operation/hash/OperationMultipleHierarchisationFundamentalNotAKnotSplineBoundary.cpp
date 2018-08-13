// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/globaldef.hpp>

#include <sgpp/optimization/operation/hash/OperationMultipleHierarchisationFundamentalNotAKnotSplineBoundary.hpp>
#include <sgpp/optimization/tools/Printer.hpp>

namespace sgpp {
namespace optimization {

OperationMultipleHierarchisationFundamentalNotAKnotSplineBoundary::
    OperationMultipleHierarchisationFundamentalNotAKnotSplineBoundary(
        base::FundamentalNotAKnotSplineBoundaryGrid& grid)
    : grid(grid), op(&grid) {}

OperationMultipleHierarchisationFundamentalNotAKnotSplineBoundary::
    ~OperationMultipleHierarchisationFundamentalNotAKnotSplineBoundary() {}

bool OperationMultipleHierarchisationFundamentalNotAKnotSplineBoundary::doHierarchisation(
    base::DataVector& nodeValues) {
  Printer::getInstance().printStatusBegin("Hierarchization (BFS)...");
  op.doHierarchisation(nodeValues);
  Printer::getInstance().printStatusEnd();
  return true;
}

void OperationMultipleHierarchisationFundamentalNotAKnotSplineBoundary::doDehierarchisation(
    base::DataVector& alpha) {
  Printer::getInstance().printStatusBegin("Dehierarchization (BFS)...");
  op.doDehierarchisation(alpha);
  Printer::getInstance().printStatusEnd();
}

bool OperationMultipleHierarchisationFundamentalNotAKnotSplineBoundary::doHierarchisation(
    base::DataMatrix& nodeValues) {
  Printer::getInstance().printStatusBegin("Hierarchization (BFS)...");
  op.doHierarchisation(nodeValues);
  Printer::getInstance().printStatusEnd();
  return true;
}

void OperationMultipleHierarchisationFundamentalNotAKnotSplineBoundary::doDehierarchisation(
    base::DataMatrix& alpha) {
  Printer::getInstance().printStatusBegin("Dehierarchization (BFS)...");
  op.doDehierarchisation(alpha);
  Printer::getInstance().printStatusEnd();
}
}  // namespace optimization
}  // namespace sgpp
