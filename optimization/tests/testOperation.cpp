// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <sgpp/base/tools/Printer.hpp>
#include <sgpp/base/tools/RandomNumberGenerator.hpp>
#include <sgpp/optimization/operation/OptimizationOpFactory.hpp>
#include <sgpp/optimization/test_problems/unconstrained/Sphere.hpp>

#include <vector>

#include "GridCreator.hpp"

using sgpp::base::Printer;
using sgpp::base::RandomNumberGenerator;
using sgpp::base::ScalarFunction;
using sgpp::optimization::OperationMultipleHierarchisation;
using sgpp::optimization::test_problems::Sphere;

BOOST_AUTO_TEST_CASE(TestOperationMultipleHierarchisation) {
  Printer::getInstance().setVerbosity(-1);
  RandomNumberGenerator::getInstance().setSeed(42);

  const size_t d = 2;
  const size_t p = 5;
  const size_t l = 4;
  const size_t m = 4;
  const double tol = 1e-4;

  Sphere testProblem(d);
  ScalarFunction& f = testProblem.getObjectiveFunction();

  // Test All The Grids!
  std::vector<std::unique_ptr<sgpp::base::Grid>> grids;
  createSupportedGrids(d, p, grids);

  for (auto& grid : grids) {
    sgpp::base::DataVector functionValues(0);
    testProblem.generateDisplacement();
    createSampleGrid(*grid, l, f, functionValues);

    sgpp::base::DataMatrix functionValuesMatrix(grid->getSize(), m);

    for (size_t j = 0; j < m; j++) {
      sgpp::base::DataVector column(0);
      testProblem.generateDisplacement();
      createSampleGrid(*grid, l, f, column);
      functionValuesMatrix.setColumn(j, column);
    }

    std::unique_ptr<OperationMultipleHierarchisation> op(
        sgpp::op_factory::createOperationMultipleHierarchisation(*grid));

    sgpp::base::DataVector functionValues2(functionValues);
    op->doHierarchisation(functionValues2);
    op->doDehierarchisation(functionValues2);

    for (size_t i = 0; i < grid->getSize(); i++) {
      BOOST_CHECK_CLOSE(functionValues[i], functionValues2[i], tol);
    }

    sgpp::base::DataMatrix functionValuesMatrix2(functionValuesMatrix);
    op->doHierarchisation(functionValuesMatrix2);
    op->doDehierarchisation(functionValuesMatrix2);

    for (size_t i = 0; i < grid->getSize(); i++) {
      for (size_t j = 0; j < m; j++) {
        BOOST_CHECK_CLOSE(functionValuesMatrix(i, j),
                          functionValuesMatrix2(i, j), tol);
      }
    }
  }
}
