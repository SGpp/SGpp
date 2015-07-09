// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/globaldef.hpp>

#include <sgpp/optimization/operation/hash/OperationMultipleHierarchisationModFundamentalSpline.hpp>
#include <sgpp/optimization/tools/Printer.hpp>

namespace SGPP {
  namespace optimization {

    OperationMultipleHierarchisationModFundamentalSpline::OperationMultipleHierarchisationModFundamentalSpline(
      base::ModFundamentalSplineGrid& grid) :
      grid(grid), op(&grid) {
    }

    OperationMultipleHierarchisationModFundamentalSpline::~OperationMultipleHierarchisationModFundamentalSpline() {
    }

    bool OperationMultipleHierarchisationModFundamentalSpline::doHierarchisation(
      base::DataVector& nodeValues) {
      printer.printStatusBegin("Hierarchization (BFS)...");
      op.doHierarchisation(nodeValues);
      printer.printStatusEnd();
      return true;
    }

    void OperationMultipleHierarchisationModFundamentalSpline::doDehierarchisation(
      base::DataVector& alpha) {
      printer.printStatusBegin("Dehierarchization (BFS)...");
      op.doDehierarchisation(alpha);
      printer.printStatusEnd();
    }

    bool OperationMultipleHierarchisationModFundamentalSpline::doHierarchisation(
      base::DataMatrix& nodeValues) {
      printer.printStatusBegin("Hierarchization (BFS)...");
      op.doHierarchisation(nodeValues);
      printer.printStatusEnd();
      return true;
    }

    void OperationMultipleHierarchisationModFundamentalSpline::doDehierarchisation(
      base::DataMatrix& alpha) {
      printer.printStatusBegin("Dehierarchization (BFS)...");
      op.doDehierarchisation(alpha);
      printer.printStatusEnd();
    }

  }
}
