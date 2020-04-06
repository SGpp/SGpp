// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/operation/BaseOpFactory.hpp>

#include <vector>

using sgpp::base::BoundingBox1D;
using sgpp::base::DataMatrix;
using sgpp::base::DataVector;
using sgpp::base::Grid;
using sgpp::base::OperationMultipleEval;
using sgpp::base::OperationEvalPartialDerivative;

BOOST_AUTO_TEST_SUITE(TestOperationMultipleEvalPartialDerivativeNaive)

BOOST_AUTO_TEST_CASE(testOperationMultipleEvalPartialDerivativeNaive) {
  const size_t d = 2;  // dims
  const size_t l = 4;  // level
  const size_t p = 3;  // degree

  std::vector<std::unique_ptr<Grid>> grids;
  grids.push_back(std::unique_ptr<Grid>(Grid::createBsplineGrid(d, p)));
  grids.push_back(std::unique_ptr<Grid>(Grid::createBsplineBoundaryGrid(d, p)));
  grids.push_back(std::unique_ptr<Grid>(Grid::createModBsplineGrid(d, p)));

  for (size_t k = 0; k < grids.size(); k++) {
    Grid& grid = *grids[k];

    // create regular sparse grid
    grid.getGenerator().regular(l);
    const size_t n = grid.getSize();

    // set random bounding box
    grid.getBoundingBox().setBoundary(0, BoundingBox1D(3.0, 5.0));
    grid.getBoundingBox().setBoundary(1, BoundingBox1D(-2.0, 2.0));

    DataVector alpha(n);

    for (int i = 0; i < static_cast<int>(n); ++i) {
      alpha[i] = static_cast<double>(i + 1);
    }

    // transformed to unit cube: {{0.5, 1.0}, {0.3, 0.4}, {0.9, 0.7}};
    const double points[3][d] = {{4.0, 2.0}, {3.6, -0.4}, {4.8, 0.8}};
    const size_t numberDataPoints = 3;

    DataVector result(numberDataPoints);
    DataMatrix dataset(numberDataPoints, d);

    for (unsigned int i = 0; i < (numberDataPoints); ++i) {
      DataVector temp(d);

      for (int j = 0; j < static_cast<int>(d); ++j) {
        temp[j] = points[i][j];
      }

      dataset.setRow(i, temp);
    }

    for (size_t dd = 0; dd < d; ++dd) {
      std::unique_ptr<OperationMultipleEval>(
          sgpp::op_factory::createOperationMultipleEvalPartialDerivativeNaive(grid, dataset, dd))
          ->mult(alpha, result);

      for (size_t i = 0; i < numberDataPoints; ++i) {
        DataVector pct(d);
        dataset.getRow(i, pct);
        double res;
        std::unique_ptr<OperationEvalPartialDerivative>(
            sgpp::op_factory::createOperationEvalPartialDerivativeNaive(grid))
            ->evalPartialDerivative(alpha, pct, dd, res);
        BOOST_CHECK_EQUAL(result[i], res);
      }
    }
  }
}

BOOST_AUTO_TEST_SUITE_END()
