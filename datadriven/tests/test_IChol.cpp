/* Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * test_IChol.cpp
 *
 *  Created on: Nov 26, 2016
 *      Author: Michael Lettrich
 */

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/grid/generation/GridGenerator.hpp>
#include <sgpp/base/operation/hash/OperationMatrix.hpp>
#include <sgpp/datadriven/algorithm/IChol.hpp>
#include <sgpp/pde/operation/PdeOpFactory.hpp>

#include <omp.h>
#include <cmath>
#include <vector>

using sgpp::base::DataMatrix;
using sgpp::base::DataVector;
using sgpp::base::Grid;
using sgpp::base::GridGenerator;
using sgpp::base::OperationMatrix;
using sgpp::datadriven::SparseDataMatrix;

BOOST_AUTO_TEST_SUITE(test_IChol)

BOOST_AUTO_TEST_CASE(decomp_identity) {
  const std::vector<double> data{1, 1, 1, 1, 1};
  const std::vector<size_t> colIdx{0, 1, 2, 3, 4};

  auto size = 5u;
  SparseDataMatrix A{size, size, data, colIdx, colIdx};

  // decomp:
  sgpp::datadriven::IChol::decompose(A, 1);

  // test
  const auto& aData = A.getDataVector();
  for (auto i = 0u; i < data.size(); i++) {
    BOOST_CHECK_CLOSE(aData[i], data[i], 10e-5);
  }
}

BOOST_AUTO_TEST_CASE(decomp_diag) {
  const std::vector<double> data{1, 4, 9, 16, 25};
  const std::vector<double> results{1, 2, 3, 4, 5};
  const std::vector<size_t> colIdx{0, 1, 2, 3, 4};

  auto size = 5u;
  SparseDataMatrix A{size, size, data, colIdx, colIdx};

  // decomp:
  sgpp::datadriven::IChol::decompose(A, 1);

  // test
  const auto& aData = A.getDataVector();
  for (auto i = 0u; i < data.size(); i++) {
    BOOST_CHECK_CLOSE(aData[i], results[i], 10e-5);
  }
}

BOOST_AUTO_TEST_CASE(decomp_arbitrary) {
  // we only get reproducable results if we run on 1 omp thread
  auto numThreads = 0;

#pragma omp parallel
  {
#pragma omp single
    { numThreads = omp_get_num_threads(); }
  }
  omp_set_num_threads(1);

  auto size = 5u;

  // initialize
  const std::vector<double> aSparse{2, 2, 6, 1, 6, 1, 6, 1, 4, 2, 16};
  const std::vector<size_t> colIdx{0, 0, 1, 0, 2, 2, 3, 0, 2, 3, 4};
  const std::vector<size_t> rowPtr{0, 1, 3, 5, 7};
  const std::vector<double> results{1.414213562373095, 1.414213562373095, 2.000000000000000,
                                    0.707106781186547, 2.345207879911715, 0.426401432711221,
                                    2.412090756622109, 0.707106781186547, 1.492405014489273,
                                    0.565333771083307, 3.599045012221992};

  SparseDataMatrix A{size, size, aSparse, colIdx, rowPtr};
  DataVector aNorm(size);
  // decomp:
  sgpp::datadriven::IChol::decompose(A, 1);

  // test
  const auto& aData = A.getDataVector();
  for (auto i = 0u; i < aData.size(); i++) {
    BOOST_CHECK_CLOSE(aData[i], results[i], 10e-5);
  }
  omp_set_num_threads(numThreads);
}

BOOST_AUTO_TEST_CASE(norm) {
  const std::vector<double> aSparse{1, 1, 4, 1, 9, 1, 16, 1, 25};
  const std::vector<double> bSparse{1.0, 1.0 / 1.0 * 1.0 / 2.0, 1.0, 1.0 / 2.0 * 1.0 / 3.0,
                                    1.0, 1.0 / 3.0 * 1.0 / 4.0, 1.0, 1.0 / 4.0 * 1.0 / 5.0,
                                    1.0};
  const double bVec[]{1.0, 1.0 / 2.0, 1.0 / 3.0, 1.0 / 4.0, 1.0 / 5.0};
  const std::vector<size_t> colIdx{0, 0, 1, 1, 2, 2, 3, 3, 4};
  const std::vector<size_t> rowPtr{0, 1, 3, 5, 7};

  auto size = 5u;

  // initialize
  SparseDataMatrix A{size, size, aSparse, colIdx, rowPtr};
  DataVector aNorm{size};
  // norm:
  sgpp::datadriven::IChol::normToUnitDiagonal(A, aNorm);

  // test
  const auto& aData = A.getDataVector();

  for (auto i = 0u; i < aData.size(); i++) {
    BOOST_CHECK_CLOSE(aData[i], bSparse[i], 10e-5);
  }

  for (auto i = 0u; i < aNorm.getSize(); i++) {
    BOOST_CHECK_CLOSE(aNorm[i], bVec[i], 10e3);
  }
}

BOOST_AUTO_TEST_CASE(reaplyNorm) {
  const std::vector<double> aSparse{1, 1, 4, 1, 9, 1, 16, 1, 25};
  const std::vector<double> bSparse{1.0, 1.0 / 1.0 * 1.0 / 2.0, 1.0, 1.0 / 2.0 * 1.0 / 3.0,
                                    1.0, 1.0 / 3.0 * 1.0 / 4.0, 1.0, 1.0 / 4.0 * 1.0 / 5.0,
                                    1.0};
  double bVec[]{1.0, 1.0 / 2.0, 1.0 / 3.0, 1.0 / 4.0, 1.0 / 5.0};
  const std::vector<size_t> colIdx{0, 0, 1, 1, 2, 2, 3, 3, 4};
  const std::vector<size_t> rowPtr{0, 1, 3, 5, 7};

  auto size = 5u;

  // initialize
  SparseDataMatrix A{size, size, bSparse, colIdx, rowPtr};
  DataVector aNorm{bVec, size};
  // norm:
  sgpp::datadriven::IChol::reaplyDiagonal(A, aNorm);

  // test
  const auto& aData = A.getDataVector();

  for (auto i = 0u; i < aData.size(); i++) {
    BOOST_CHECK_CLOSE(aData[i], aSparse[i], 10e-5);
  }
}

BOOST_AUTO_TEST_CASE(normSGMatrix) {
  auto dim = 5u;
  auto lvl = 3u;
  auto grid = std::unique_ptr<Grid>{Grid::createLinearGrid(dim)};
  grid->getGenerator().regular(lvl);
  auto size = grid->getSize();
  DataMatrix lhsMatrix{size, size};
  std::unique_ptr<OperationMatrix> op(
      sgpp::op_factory::createOperationLTwoDotExplicit(&lhsMatrix, *grid));

  // extract lower triangular matrix.
  for (size_t i = 0; i < lhsMatrix.getNrows() - 1; i++) {
    for (size_t j = i + 1; j < lhsMatrix.getNcols(); j++) {
      lhsMatrix.set(i, j, 0);
    }
  }
  SparseDataMatrix sparseLHS(lhsMatrix);
  DataVector norm{size};
  sgpp::datadriven::IChol::normToUnitDiagonal(sparseLHS, norm);
  sgpp::datadriven::IChol::reaplyDiagonal(sparseLHS, norm);

  DataMatrix normed{size, size};
  SparseDataMatrix::toDataMatrix(sparseLHS, normed);

  // test
  for (auto i = 0u; i < normed.getSize(); i++) {
    BOOST_CHECK_CLOSE(normed[i], lhsMatrix[i], 10e-5);
  }
}

BOOST_AUTO_TEST_SUITE_END()
