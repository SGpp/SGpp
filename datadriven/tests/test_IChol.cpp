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
#include <sgpp/datadriven/algorithm/IChol.hpp>

#include <cmath>
#include <vector>

using sgpp::base::DataMatrix;
using sgpp::base::DataVector;
using sgpp::datadriven::IChol;
using sgpp::datadriven::SparseDataMatrix;

BOOST_AUTO_TEST_SUITE(test_IChol)

BOOST_AUTO_TEST_CASE(decomp_identity) {
  const std::vector<double> data{1, 1, 1, 1, 1};
  const std::vector<size_t> colIdx{0, 1, 2, 3, 4};

  auto size = 5u;
  SparseDataMatrix A{size, size, data, colIdx, colIdx};

  // decomp:
  IChol::decompose(A, 1);

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
  IChol::decompose(A, 1);

  // test
  const auto& aData = A.getDataVector();
  for (auto i = 0u; i < data.size(); i++) {
    BOOST_CHECK_CLOSE(aData[i], results[i], 10e-5);
  }
}

BOOST_AUTO_TEST_CASE(decomp_arbitrary) {
  auto size = 5u;

  // initialize
  const std::vector<double> aSparse{2, 2, 6, 1, 6, 1, 6, 1, 4, 2, 16};
  const std::vector<size_t> colIdx{0, 0, 1, 0, 2, 2, 3, 0, 2, 3, 4};
  const std::vector<size_t> rowPtr{0, 1, 3, 5, 7};
  const std::vector<double> results{1.4142, 1.4142, 2.0000, 0.7071, 2.3452, 0.4264,
                                    2.4121, 0.7071, 1.4924, 0.5653, 3.5990};

  SparseDataMatrix A{size, size, aSparse, colIdx, rowPtr};
  DataVector aNorm(size);
  // decomp:
  IChol::normToUnitDiagonal(A, aNorm);
  IChol::decompose(A, 1);
  // IChol::reaplyDiagonal(A, aNorm);

  // test
  const auto& aData = A.getDataVector();
  for (auto i = 0u; i < aData.size(); i++) {
    BOOST_CHECK_CLOSE(aData[i], results[i], 10e-5);
  }
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
  IChol::normToUnitDiagonal(A, aNorm);

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
  IChol::reaplyDiagonal(A, aNorm);

  // test
  const auto& aData = A.getDataVector();

  for (auto i = 0u; i < aData.size(); i++) {
    BOOST_CHECK_CLOSE(aData[i], aSparse[i], 10e-5);
  }
}

BOOST_AUTO_TEST_SUITE_END()
