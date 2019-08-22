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
#include <sgpp/base/operation/hash/OperationDiagonal.hpp>

#include <array>

using sgpp::base::DataMatrix;
using sgpp::base::DataVector;
using sgpp::base::OperationDiagonal;
using sgpp::base::OperationMatrix;

BOOST_AUTO_TEST_SUITE(TestOperationDiagonal)

BOOST_AUTO_TEST_CASE(testOperationDiagonalIdentity) {
  size_t dim = 3;
  std::unique_ptr<sgpp::base::Grid> grid(sgpp::base::Grid::createLinearGrid(dim));
  auto& gen = grid->getGenerator();
  auto& gridStorage = grid->getStorage();
  gen.regular(2);

  double multFactor = 1.0;
  // OperationDiagonal with multFactor of 1.0 should be identical to an identity matrix!
  std::unique_ptr<OperationMatrix> opDiag(
      sgpp::op_factory::createOperationDiagonal(*grid, multFactor));

  const auto size = gridStorage.getSize();
  auto testVec = sgpp::base::DataVector(size, 42.0);
  auto resultVec = sgpp::base::DataVector(size);

  opDiag->mult(testVec, resultVec);

  for (size_t i = 0; i < size; ++i) {
    BOOST_CHECK_CLOSE(testVec[i], resultVec[i], 1e-7);
  }
}

BOOST_AUTO_TEST_CASE(testOperationDiagonal) {
  size_t dim = 2;
  std::unique_ptr<sgpp::base::Grid> grid(sgpp::base::Grid::createLinearGrid(dim));
  auto& gen = grid->getGenerator();
  auto& gridStorage = grid->getStorage();
  gen.regular(2);

  std::unique_ptr<OperationMatrix> opDiag(sgpp::op_factory::createOperationDiagonal(*grid));

  const auto size = gridStorage.getSize();
  auto testVec = sgpp::base::DataVector(size, 1.0);
  auto resultVec = sgpp::base::DataVector(size);
  const auto shouldVec = std::array<double, 5>{{1, 4.0, 4.0, 4.0, 4.0}};

  opDiag->mult(testVec, resultVec);

  for (size_t i = 0; i < size; ++i) {
    BOOST_CHECK_CLOSE(resultVec[i], shouldVec[i], 1e-7);
  }
}

BOOST_AUTO_TEST_SUITE_END()
