// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/tools/Printer.hpp>
#include <sgpp/base/tools/RandomNumberGenerator.hpp>
#include <sgpp/base/tools/sle/system/FullSLE.hpp>
#include <sgpp/optimization/gridgen/IterativeGridGeneratorRitterNovak.hpp>
#include <sgpp/optimization/test_problems/unconstrained/Sphere.hpp>
#include <sgpp/optimization/tools/FileIO.hpp>
#include <sgpp/optimization/tools/Math.hpp>

#include <string>
#include <vector>

#include "ObjectiveFunctions.hpp"

using sgpp::base::Printer;
using sgpp::base::RandomNumberGenerator;

void gridEqualityTest(sgpp::base::Grid& grid1, sgpp::base::Grid& grid2) {
  sgpp::base::GridStorage& storage1 = grid1.getStorage();
  sgpp::base::GridStorage& storage2 = grid2.getStorage();
  const size_t d = storage1.getDimension();
  BOOST_CHECK_EQUAL(d, storage2.getDimension());
  const size_t n = storage1.getSize();
  BOOST_CHECK_EQUAL(n, storage2.getSize());

  for (size_t k = 0; k < n; k++) {
    for (size_t t = 0; t < d; t++) {
      BOOST_CHECK_EQUAL(storage1[k].getLevel(t), storage2[k].getLevel(t));
      BOOST_CHECK_EQUAL(storage1[k].getIndex(t), storage2[k].getIndex(t));
    }
  }
}

void orthogonalityTest(sgpp::base::DataMatrix& A) {
  const size_t n = A.getNrows();
  BOOST_CHECK_EQUAL(n, A.getNcols());

  for (size_t i = 0; i < n; i++) {
    for (size_t j = 0; j <= i; j++) {
      double entry = 0.0;

      for (size_t l = 0; l < n; l++) {
        entry += A(i, l) * A(j, l);
      }

      BOOST_CHECK_SMALL(entry - ((i == j) ? 1.0 : 0.0), 1e-10);
    }
  }
}

void symmetryTest(sgpp::base::DataMatrix& A) {
  const size_t n = A.getNrows();
  BOOST_CHECK_EQUAL(n, A.getNcols());

  for (size_t i = 0; i < n; i++) {
    for (size_t j = 0; j < i; j++) {
      BOOST_CHECK_CLOSE(A(i, j), A(j, i), 1e-10);
    }
  }
}

void similiarityTest(sgpp::base::DataMatrix& A, sgpp::base::DataMatrix& V,
                     sgpp::base::DataMatrix& B) {
  const size_t n = A.getNrows();
  BOOST_CHECK_EQUAL(n, A.getNcols());
  BOOST_CHECK_EQUAL(n, V.getNrows());
  BOOST_CHECK_EQUAL(n, V.getNcols());
  BOOST_CHECK_EQUAL(n, B.getNrows());
  BOOST_CHECK_EQUAL(n, B.getNcols());

  for (size_t i = 0; i < n; i++) {
    for (size_t j = 0; j < n; j++) {
      double entry1 = 0.0, entry2 = 0.0;

      for (size_t l = 0; l < n; l++) {
        entry1 += A(i, l) * V(l, j);
        entry2 += V(i, l) * B(l, j);
      }

      BOOST_CHECK_CLOSE(entry1, entry2, 5e-5);
    }
  }
}

void generateRandomMatrix(sgpp::base::DataMatrix& A) {
  for (size_t i = 0; i < A.getNrows(); i++) {
    for (size_t j = 0; j < A.getNcols(); j++) {
      A(i, j) = RandomNumberGenerator::getInstance().getGaussianRN();
    }
  }
}

template <class T>
void randomMatrixEntry(T& x) {
  x = static_cast<T>(RandomNumberGenerator::getInstance().getUniformIndexRN(100));
}

template <>
void randomMatrixEntry(float& x) {
  x = static_cast<float>(RandomNumberGenerator::getInstance().getUniformRN(-100.0, 100.0));
}

template <>
void randomMatrixEntry(double& x) {
  x = static_cast<double>(RandomNumberGenerator::getInstance().getUniformRN(-100.0, 100.0));
}

template <>
void randomMatrixEntry(std::string& x) {
  const size_t length = RandomNumberGenerator::getInstance().getUniformIndexRN(100);
  x.clear();

  for (size_t i = 0; i < length; i++) {
    x += static_cast<char>(32 + RandomNumberGenerator::getInstance().getUniformIndexRN(96));
  }
}

template <class T>
void testReadWriteMatrix(std::vector<T>& A1, size_t m1, size_t n1) {
  std::vector<T> A2;
  A1.clear();

  for (size_t i = 0; i < m1; i++) {
    for (size_t j = 0; j < n1; j++) {
      T entry;
      randomMatrixEntry(entry);
      A1.push_back(entry);
    }
  }

  size_t m2, n2;

  {
    const std::string fileName = "testTools_matrix.tmp";
    sgpp::optimization::file_io::writeMatrix(fileName, A1, m1, n1);
    sgpp::optimization::file_io::readMatrix(fileName, A2, m2, n2);
    std::remove(fileName.c_str());
  }

  BOOST_CHECK_EQUAL(m1, m2);
  BOOST_CHECK_EQUAL(n1, n2);
  BOOST_CHECK_EQUAL(A1.size(), A2.size());

  for (size_t k = 0; k < A1.size(); k++) {
    BOOST_CHECK_EQUAL(A1[k], A2[k]);
  }
}

BOOST_AUTO_TEST_CASE(TestFileIOReadWriteGrid) {
  // Test sgpp::optimization::sgpp::optimization::file_io::readGrid/writeGrid.
  Printer::getInstance().setVerbosity(-1);
  RandomNumberGenerator::getInstance().setSeed(42);

  for (size_t d = 1; d < 6; d++) {
    std::unique_ptr<sgpp::base::Grid> grid1(sgpp::base::Grid::createLinearGrid(d)),
        grid2(sgpp::base::Grid::createLinearGrid(d));
    grid1->getGenerator().regular(3);

    {
      const std::string fileName = "testTools_grid.tmp";
      sgpp::optimization::file_io::writeGrid(fileName, grid1->getStorage());
      sgpp::optimization::file_io::readGrid(fileName, grid2->getStorage());
      std::remove(fileName.c_str());
    }

    gridEqualityTest(*grid1, *grid2);

    sgpp::base::DataVector functionValues1(grid1->getSize());
    sgpp::base::DataVector functionValues2(0);

    for (size_t k = 0; k < functionValues1.getSize(); k++) {
      functionValues1[k] = RandomNumberGenerator::getInstance().getUniformRN();
    }

    {
      const std::string fileName = "testTools_grid.tmp";
      sgpp::optimization::file_io::writeGrid(fileName, grid1->getStorage(), functionValues1);
      std::unique_ptr<sgpp::base::Grid> grid2(sgpp::base::Grid::createLinearGrid(d));
      sgpp::optimization::file_io::readGrid(fileName, grid2->getStorage(), functionValues2);
      std::remove(fileName.c_str());
    }

    gridEqualityTest(*grid1, *grid2);
    BOOST_CHECK_EQUAL(functionValues1.getSize(), functionValues2.getSize());

    for (size_t k = 0; k < functionValues1.getSize(); k++) {
      BOOST_CHECK_EQUAL(functionValues1[k], functionValues2[k]);
    }
  }
}

BOOST_AUTO_TEST_CASE(TestFileIOReadWriteMatrix) {
  // Test sgpp::optimization::sgpp::optimization::file_io::readMatrix/writeMatrix.
  Printer::getInstance().setVerbosity(-1);
  RandomNumberGenerator::getInstance().setSeed(42);
  const size_t m1 = 100;
  const size_t n1 = 200;

  // test read/write with std::vector<T>
  {
    std::vector<double> A1;
    testReadWriteMatrix(A1, m1, n1);
  }

  {
    std::vector<float> A1;
    testReadWriteMatrix(A1, m1, n1);
  }

  {
    std::vector<double> A1;
    testReadWriteMatrix(A1, m1, n1);
  }

  {
    std::vector<uint8_t> A1;
    testReadWriteMatrix(A1, m1, n1);
  }

  {
    std::vector<uint16_t> A1;
    testReadWriteMatrix(A1, m1, n1);
  }

  {
    std::vector<uint32_t> A1;
    testReadWriteMatrix(A1, m1, n1);
  }

  {
    std::vector<uint64_t> A1;
    testReadWriteMatrix(A1, m1, n1);
  }

  {
    std::vector<std::string> A1;
    testReadWriteMatrix(A1, m1, n1);
  }

  // test read/write with DataMatrix
  {
    sgpp::base::DataMatrix A1(m1, n1), A2(0, 0);

    for (size_t i = 0; i < m1; i++) {
      for (size_t j = 0; j < n1; j++) {
        A1(i, j) = RandomNumberGenerator::getInstance().getUniformRN();
      }
    }

    {
      const std::string fileName = "testTools_matrix.tmp";
      sgpp::optimization::file_io::writeMatrix(fileName, A1);
      sgpp::optimization::file_io::readMatrix(fileName, A2);
      std::remove(fileName.c_str());
    }

    BOOST_CHECK_EQUAL(A1.getNrows(), A2.getNrows());
    BOOST_CHECK_EQUAL(A1.getNcols(), A2.getNcols());

    for (size_t i = 0; i < m1; i++) {
      for (size_t j = 0; j < n1; j++) {
        BOOST_CHECK_EQUAL(A1(i, j), A2(i, j));
      }
    }
  }
}

BOOST_AUTO_TEST_CASE(TestFileIOReadWriteVector) {
  // Test sgpp::optimization::sgpp::optimization::file_io::readVector/writeVector.
  Printer::getInstance().setVerbosity(-1);
  RandomNumberGenerator::getInstance().setSeed(42);
  const size_t n = 100;

  // test read/write with std::vector<double>
  {
    std::vector<double> v1, v2;

    for (size_t i = 0; i < n; i++) {
      v1.push_back(RandomNumberGenerator::getInstance().getUniformRN());
    }

    {
      const std::string fileName = "testTools_matrix.tmp";
      sgpp::optimization::file_io::writeVector(fileName, v1);
      sgpp::optimization::file_io::readVector(fileName, v2);
      std::remove(fileName.c_str());
    }

    BOOST_CHECK_EQUAL(v1.size(), v2.size());

    for (size_t k = 0; k < v1.size(); k++) {
      BOOST_CHECK_EQUAL(v1[k], v2[k]);
    }
  }

  // test read/write with DataVector
  {
    sgpp::base::DataVector v1(n), v2(0);

    for (size_t i = 0; i < n; i++) {
      v1[i] = RandomNumberGenerator::getInstance().getUniformRN();
    }

    {
      const std::string fileName = "testTools_vector.tmp";
      sgpp::optimization::file_io::writeVector(fileName, v1);
      sgpp::optimization::file_io::readVector(fileName, v2);
      std::remove(fileName.c_str());
    }

    BOOST_CHECK_EQUAL(v1.getSize(), v2.getSize());

    for (size_t i = 0; i < v1.getSize(); i++) {
      BOOST_CHECK_EQUAL(v1[i], v2[i]);
    }
  }
}

BOOST_AUTO_TEST_CASE(TestHouseholderTransformation) {
  // Test sgpp::optimization::sgpp::optimization::math::householderTransformation.
  RandomNumberGenerator::getInstance().setSeed(42);
  const size_t n = 20;
  sgpp::base::DataMatrix A(n, n);
  generateRandomMatrix(A);

  const std::vector<size_t> ps = {0, 2, 5};
  const std::vector<size_t> qs = {0, 1, 9};

  for (size_t k = 0; k < ps.size(); k++) {
    const size_t p = ps[k];
    const size_t q = qs[k];
    const size_t m = n - p;
    sgpp::base::DataMatrix Q(m, m);
    sgpp::optimization::math::householderTransformation(A, p, q, Q);
    symmetryTest(Q);
    orthogonalityTest(Q);

    for (size_t i = p + 1; i < n; i++) {
      double entry = 0.0;

      for (size_t l = p; l < n; l++) {
        entry += Q(i - p, l - p) * A(l, q);
      }

      BOOST_CHECK_SMALL(entry, 1e-10);
    }
  }
}

BOOST_AUTO_TEST_CASE(TestHessenbergForm) {
  // Test sgpp::optimization::sgpp::optimization::math::hessenbergForm.
  RandomNumberGenerator::getInstance().setSeed(42);
  const size_t n = 20;
  sgpp::base::DataMatrix A(n, n);
  generateRandomMatrix(A);
  sgpp::base::DataMatrix H(A);
  sgpp::base::DataMatrix V(n, n);
  sgpp::optimization::math::hessenbergForm(H, V);

  orthogonalityTest(V);
  similiarityTest(A, V, H);

  for (size_t i = 2; i < n; i++) {
    for (size_t j = 0; j < i - 1; j++) {
      BOOST_CHECK_SMALL(H(i, j), 1e-10);
    }
  }
}

BOOST_AUTO_TEST_CASE(TestQRDecomposition) {
  // Test sgpp::optimization::sgpp::optimization::math::QRDecomposition.
  RandomNumberGenerator::getInstance().setSeed(42);
  const size_t n = 20;
  sgpp::base::DataMatrix A(n, n);
  generateRandomMatrix(A);
  sgpp::base::DataMatrix R(A);
  sgpp::base::DataMatrix Q(n, n);
  sgpp::optimization::math::QRDecomposition(R, Q);

  orthogonalityTest(Q);

  for (size_t i = 1; i < n; i++) {
    for (size_t j = 0; j < i - 1; j++) {
      BOOST_CHECK_SMALL(R(i, j), 1e-10);
    }
  }

  for (size_t i = 0; i < n; i++) {
    for (size_t j = 0; j < n; j++) {
      double entry = 0.0;

      for (size_t l = 0; l < n; l++) {
        entry += Q(i, l) * R(l, j);
      }

      BOOST_CHECK_CLOSE(A(i, j), entry, 1e-10);
    }
  }
}

BOOST_AUTO_TEST_CASE(TestSchurDecomposition) {
  // Test sgpp::optimization::sgpp::optimization::math::schurDecomposition.
  RandomNumberGenerator::getInstance().setSeed(42);
  const size_t n = 20;
  sgpp::base::DataMatrix A(n, n);
  generateRandomMatrix(A);

  for (size_t i = 0; i < n; i++) {
    for (size_t j = 0; j < i; j++) {
      A(i, j) = A(j, i);
    }
  }

  sgpp::base::DataMatrix S(A);
  sgpp::base::DataMatrix V(n, n);
  sgpp::optimization::math::schurDecomposition(S, V);

  orthogonalityTest(V);
  similiarityTest(A, V, S);

  for (size_t i = 0; i < n; i++) {
    for (size_t j = 0; j < n; j++) {
      if (i != j) {
        BOOST_CHECK_SMALL(S(i, j), 1e-8);
      }
    }
  }
}
