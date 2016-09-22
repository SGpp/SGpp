// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/datatypes/DataMatrix.hpp>
// #include <sgpp/datadriven/DatadrivenOpFactory.hpp>
#include <sgpp/base/operation/BaseOpFactory.hpp>

using sgpp::base::BoundingBox1D;
using sgpp::base::DataMatrix;
using sgpp::base::DataVector;
using sgpp::base::Grid;
using sgpp::base::GridStorage;
using sgpp::base::OperationMultipleEval;

BOOST_AUTO_TEST_SUITE(TestOperationMultipleEval)

BOOST_AUTO_TEST_CASE(testOperationMultipleEval) {
  const size_t dim = 2;
  std::unique_ptr<Grid> grid(Grid::createLinearGrid(dim));
  grid->getGenerator().regular(2);

  grid->getBoundingBox().setBoundary(0, BoundingBox1D(3.0, 5.0));
  grid->getBoundingBox().setBoundary(1, BoundingBox1D(-2.0, 2.0));

  GridStorage& gS = grid->getStorage();

  size_t N = gS.getSize();

  DataVector alpha(N);

  for (int i = 0; i < static_cast<int>(N); ++i) {
    alpha[i] = static_cast<double>(i + 1);
  }

  // transformed to unit cube: {{0.5, 1.0}, {0.3, 0.4}, {0.9, 0.7}};
  const double points[3][2] = {{4.0, 2.0}, {3.6, -0.4}, {4.8, 0.8}};
  const size_t numberDataPoints = 3;

  DataVector result(numberDataPoints);

  for (unsigned int i = 0; i < (numberDataPoints); ++i) {
    result[i] = static_cast<double>(i);
  }

  DataMatrix dataset(numberDataPoints, dim);

  for (unsigned int i = 0; i < (numberDataPoints); ++i) {
    DataVector temp(dim);

    for (int j = 0; j < static_cast<int>(dim); ++j) {
      temp[j] = points[i][j];
    }

    dataset.setRow(i, temp);
  }

  sgpp::op_factory::createOperationMultipleEval(*grid, dataset)->mult(alpha, result);

  BOOST_TEST_MESSAGE(alpha.toString() + "\n");
  BOOST_TEST_MESSAGE(result.toString() + "\n");

  DataVector result_ref(numberDataPoints);
  result_ref[0] = 0.0;
  result_ref[1] = 2.72;
  result_ref[2] = 1.64;

  BOOST_CHECK_CLOSE(result[0], result_ref[0], 1e-7);
  BOOST_CHECK_CLOSE(result[1], result_ref[1], 1e-7);
  BOOST_CHECK_CLOSE(result[2], result_ref[2], 1e-7);
}

BOOST_AUTO_TEST_SUITE_END()
