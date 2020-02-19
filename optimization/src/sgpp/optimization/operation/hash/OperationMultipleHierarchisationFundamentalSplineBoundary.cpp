// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/globaldef.hpp>

#include <sgpp/base/tools/Printer.hpp>
#include <sgpp/optimization/operation/hash/OperationMultipleHierarchisationFundamentalSplineBoundary.hpp>

namespace sgpp {
namespace optimization {

OperationMultipleHierarchisationFundamentalSplineBoundary::
    OperationMultipleHierarchisationFundamentalSplineBoundary(
        base::FundamentalSplineBoundaryGrid& grid)
    : grid(grid), op(&grid) {}

OperationMultipleHierarchisationFundamentalSplineBoundary::
    ~OperationMultipleHierarchisationFundamentalSplineBoundary() {}

bool OperationMultipleHierarchisationFundamentalSplineBoundary::doHierarchisation(
    base::DataVector& nodeValues) {
  base::Printer::getInstance().printStatusBegin("Hierarchization (BFS)...");
  op.doHierarchisation(nodeValues);
  base::Printer::getInstance().printStatusEnd();
  return true;
}

void OperationMultipleHierarchisationFundamentalSplineBoundary::doDehierarchisation(
    base::DataVector& alpha) {
  base::Printer::getInstance().printStatusBegin("Dehierarchization (BFS)...");
  op.doDehierarchisation(alpha);
  base::Printer::getInstance().printStatusEnd();
}

bool OperationMultipleHierarchisationFundamentalSplineBoundary::doHierarchisation(
    base::DataMatrix& nodeValues) {
  base::Printer::getInstance().printStatusBegin("Hierarchization (BFS)...");
  op.doHierarchisation(nodeValues);
  base::Printer::getInstance().printStatusEnd();
  return true;
}

void OperationMultipleHierarchisationFundamentalSplineBoundary::doDehierarchisation(
    base::DataMatrix& alpha) {
  base::Printer::getInstance().printStatusBegin("Dehierarchization (BFS)...");
  op.doDehierarchisation(alpha);
  base::Printer::getInstance().printStatusEnd();
}
}  // namespace optimization
}  // namespace sgpp
