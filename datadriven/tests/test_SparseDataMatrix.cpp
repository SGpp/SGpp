/* Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * test_SparseDataMatrix.cpp
 *
 *  Created on: Feb 5, 2017
 *      Author: Michael Lettrich
 */

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/datadriven/algorithm/SparseDataMatrix.hpp>

#include <algorithm>
#include <iostream>
#include <vector>

using sgpp::datadriven::SparseDataMatrix;
using sgpp::base::DataMatrix;

BOOST_AUTO_TEST_SUITE(test_SparseDataMatrix)

BOOST_AUTO_TEST_CASE(testToDataMatrix) {
  const double a_vec[]{-2, 1, 0, 0, 0, 1,  -2, 1, 0, 0, 0, 1, -2,
                       1,  0, 0, 0, 1, -2, 1,  0, 0, 0, 1, -2};
  const std::vector<double> a_sparse{-2, 1, 1, -2, 1, 1, -2, 1, 1, -2, 1, 1, -2};
  const std::vector<size_t> a_colIdx{0, 1, 0, 1, 2, 1, 2, 3, 2, 3, 4, 3, 4};
  const std::vector<size_t> a_rowPtr{0, 2, 5, 8, 11};

  SparseDataMatrix A{5, 5};
  auto& aData = A.getDataVector();
  auto& aColIdx = A.getColIndexVector();
  auto& aRowPtr = A.getRowPtrVector();

  aData.insert(std::begin(aData), std::begin(a_sparse), std::end(a_sparse));
  aColIdx.insert(std::begin(aColIdx), std::begin(a_colIdx), std::end(a_colIdx));
  std::copy(std::begin(a_rowPtr), std::end(a_rowPtr), std::begin(aRowPtr));
  //  aData = a_sparse;
  //  aColIdx = a_colIdx;
  //  aRowPtr = a_rowPtr;

  DataMatrix B{};

  SparseDataMatrix::toDataMatrix(A, B);

  BOOST_CHECK_EQUAL(B.getNrows(), 5);
  BOOST_CHECK_EQUAL(B.getNcols(), 5);
  BOOST_CHECK_EQUAL(B.getSize(), 25);
  for (auto i = 0u; i < B.getSize(); i++) {
    BOOST_CHECK_CLOSE(a_vec[i], B[i], 10 - 5);
  }
}

BOOST_AUTO_TEST_CASE(testResize_shrink) {
  const double a_vec[]{-2, 1, 0, 0, 1, -2, 1, 0, 0, 1, -2, 1, 0, 0, 1, -2};
  const std::vector<double> a_sparse{-2, 1, 1, -2, 1, 1, -2, 1, 1, -2, 1, 1, -2};
  const std::vector<size_t> a_colIdx{0, 1, 0, 1, 2, 1, 2, 3, 2, 3, 4, 3, 4};
  const std::vector<size_t> a_rowPtr{0, 2, 5, 8, 11};

  SparseDataMatrix A{5, 5};
  auto& aData = A.getDataVector();
  auto& aColIdx = A.getColIndexVector();
  auto& aRowPtr = A.getRowPtrVector();

  aData.insert(std::begin(aData), std::begin(a_sparse), std::end(a_sparse));
  aColIdx.insert(std::begin(aColIdx), std::begin(a_colIdx), std::end(a_colIdx));
  std::copy(std::begin(a_rowPtr), std::end(a_rowPtr), std::begin(aRowPtr));

  A.resize(4, 4);

  DataMatrix B{};
  SparseDataMatrix::toDataMatrix(A, B);

  std::cout << "data = [";
  for (auto item : aData) {
    std::cout << item << ",";
  }
  std::cout << "]\n";

  std::cout << "colIdx = [";
  for (auto item : aColIdx) {
    std::cout << item << ",";
  }
  std::cout << "]\n";

  std::cout << "rowPtr = [";
  for (auto item : aRowPtr) {
    std::cout << item << ",";
  }
  std::cout << "]\n";

  printf("%s\n", B.toString().c_str());

  BOOST_CHECK_EQUAL(B.getNrows(), 4);
  BOOST_CHECK_EQUAL(B.getNcols(), 4);
  BOOST_CHECK_EQUAL(B.getSize(), 16);
  for (auto i = 0u; i < B.getSize(); i++) {
    BOOST_CHECK_CLOSE(a_vec[i], B[i], 10 - 5);
  }
}

BOOST_AUTO_TEST_SUITE_END()
