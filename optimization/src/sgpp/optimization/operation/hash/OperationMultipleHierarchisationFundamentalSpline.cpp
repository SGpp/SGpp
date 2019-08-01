// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/tools/Printer.hpp>
#include <sgpp/globaldef.hpp>
#include <sgpp/optimization/operation/hash/OperationMultipleHierarchisationFundamentalSpline.hpp>

namespace sgpp {
namespace optimization {

OperationMultipleHierarchisationFundamentalSpline::
    OperationMultipleHierarchisationFundamentalSpline(base::FundamentalSplineGrid& grid)
    : grid(grid), op(&grid) {}

OperationMultipleHierarchisationFundamentalSpline::
    ~OperationMultipleHierarchisationFundamentalSpline() {}

bool OperationMultipleHierarchisationFundamentalSpline::doHierarchisation(
    base::DataVector& nodeValues) {
  base::Printer::getInstance().printStatusBegin("Hierarchization (BFS)...");
  op.doHierarchisation(nodeValues);
  base::Printer::getInstance().printStatusEnd();
  return true;
}

void OperationMultipleHierarchisationFundamentalSpline::doDehierarchisation(
    base::DataVector& alpha) {
  base::Printer::getInstance().printStatusBegin("Dehierarchization (BFS)...");
  op.doDehierarchisation(alpha);
  base::Printer::getInstance().printStatusEnd();
}

bool OperationMultipleHierarchisationFundamentalSpline::doHierarchisation(
    base::DataMatrix& nodeValues) {
  base::Printer::getInstance().printStatusBegin("Hierarchization (BFS)...");
  op.doHierarchisation(nodeValues);
  base::Printer::getInstance().printStatusEnd();
  return true;
}

void OperationMultipleHierarchisationFundamentalSpline::doDehierarchisation(
    base::DataMatrix& alpha) {
  base::Printer::getInstance().printStatusBegin("Dehierarchization (BFS)...");
  op.doDehierarchisation(alpha);
  base::Printer::getInstance().printStatusEnd();
}
}  // namespace optimization
}  // namespace sgpp
