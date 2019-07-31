// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifdef ZLIB

#define BOOST_TEST_DYN_LINK
#include <boost/lexical_cast.hpp>
#include <boost/test/unit_test.hpp>

#include <zlib.h>

#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/base/grid/generation/functors/SurplusRefinementFunctor.hpp>
#include <sgpp/base/operation/BaseOpFactory.hpp>
#include <sgpp/globaldef.hpp>
#include <sgpp/pde/operation/PdeOpFactory.hpp>

#include <algorithm>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

using sgpp::base::DataMatrix;
using sgpp::base::DataVector;
using sgpp::base::Grid;
using sgpp::base::GridGenerator;
using sgpp::base::GridStorage;
using sgpp::base::HashGridPoint;
using sgpp::base::OperationMatrix;
using sgpp::base::SurplusRefinementFunctor;

DataMatrix* generateLaplaceMatrix(Grid& grid, size_t level) {
  GridStorage& storage = grid.getStorage();

  grid.getGenerator().regular(level);

  std::unique_ptr<OperationMatrix> laplace(sgpp::op_factory::createOperationLaplace(grid));

  // create vector
  DataVector alpha(storage.getSize());
  DataVector erg(storage.getSize());

  // create stiffness matrix
  DataMatrix* m = new DataMatrix(storage.getSize(), storage.getSize());
  m->setAll(0);

  for (size_t i = 0; i < storage.getSize(); i++) {
    alpha.setAll(0.0);
    alpha[i] = 1.0;
    laplace->mult(alpha, erg);
    m->setColumn(i, erg);
  }

  return m;
}

DataMatrix* generateLaplaceEnhancedMatrix(Grid& grid, size_t level) {
  GridStorage& storage = grid.getStorage();

  grid.getGenerator().regular(level);

  std::unique_ptr<OperationMatrix> laplace(sgpp::op_factory::createOperationLaplaceEnhanced(grid));

  // create vector
  DataVector alpha(storage.getSize());
  DataVector erg(storage.getSize());

  // create stiffness matrix
  DataMatrix* m = new DataMatrix(storage.getSize(), storage.getSize());
  m->setAll(0);

  for (size_t i = 0; i < storage.getSize(); i++) {
    alpha.setAll(0.0);
    alpha[i] = 1.0;
    laplace->mult(alpha, erg);
    m->setColumn(i, erg);
  }

  return m;
}

std::string uncompressFile(std::string fileName) {
  gzFile inFileZ = gzopen(fileName.c_str(), "rb");

  if (inFileZ == NULL) {
    std::cout << "Error: Failed to gzopen file " << fileName << std::endl;
    exit(0);
  }

  unsigned char unzipBuffer[8192];
  unsigned int unzippedBytes;
  std::vector<unsigned char> unzippedData;

  while (true) {
    unzippedBytes = gzread(inFileZ, unzipBuffer, 8192);

    if (unzippedBytes > 0) {
      for (size_t i = 0; i < unzippedBytes; i++) {
        unzippedData.push_back(unzipBuffer[i]);
      }
    } else {
      break;
    }
  }

  gzclose(inFileZ);

  std::stringstream convert;

  for (size_t i = 0; i < unzippedData.size(); i++) {
    convert << unzippedData[i];
  }

  return convert.str();
}

DataMatrix* readReferenceMatrix(GridStorage& storage, std::string fileName) {
  std::string content = uncompressFile(fileName);

  std::stringstream contentStream;
  contentStream << content;
  std::string line;

  DataMatrix* m = new DataMatrix(0, storage.getSize());

  size_t currentRow = 0;

  while (!contentStream.eof()) {
    std::getline(contentStream, line);

    // for lines that only contain a newline
    if (line.size() == 0) {
      break;
    }

    m->appendRow();

    size_t curPos = 0;
    size_t curFind = 0;
    std::string curValue;
    double floatValue;

    for (size_t i = 0; i < storage.getSize(); i++) {
      curFind = std::min(line.find(" ", curPos), line.find("\t", curPos));
      curValue = line.substr(curPos, curFind - curPos);
      floatValue = boost::lexical_cast<double>(curValue);
      m->set(currentRow, i, floatValue);
      curPos = curFind + 1;
    }

    currentRow += 1;
  }

  return m;
}

void compareStiffnessMatrices(DataMatrix* m1, DataMatrix* m2) {
  double tolerance = 1E-5;

  // check dimensions
  BOOST_CHECK_EQUAL(m1->getNrows(), m2->getNrows());
  BOOST_CHECK_EQUAL(m1->getNcols(), m2->getNcols());

  size_t rows = m1->getNrows();  // was n

  size_t cols = m1->getNcols();  // was m

  // check diagonal
  std::vector<double> values;

  for (size_t i = 0; i < rows; i++) {
    values.push_back(m1->get(i, i));
  }

  std::sort(values.begin(), values.end());
  std::vector<double> valuesRef;

  for (size_t i = 0; i < rows; i++) {
    valuesRef.push_back(m2->get(i, i));
  }

  std::sort(valuesRef.begin(), valuesRef.end());

  for (size_t i = 0; i < rows; i++) {
    BOOST_CHECK_CLOSE(values[i], valuesRef[i], tolerance);
  }

  // check row sum
  DataVector v(cols);

  values.clear();

  for (size_t i = 0; i < rows; i++) {
    m1->getRow(i, v);
    values.push_back(v.sum());
  }

  std::sort(values.begin(), values.end());

  std::vector<double> valuesReference;

  for (size_t i = 0; i < rows; i++) {
    m2->getRow(i, v);
    valuesReference.push_back(v.sum());
  }

  std::sort(valuesReference.begin(), valuesReference.end());

  for (size_t i = 0; i < rows; i++) {
    BOOST_CHECK_SMALL(values[i] - valuesReference[i], tolerance);
  }

  // check col sum
  values.clear();

  DataVector vRows(rows);

  for (size_t i = 0; i < cols; i++) {
    m1->getColumn(i, vRows);
    values.push_back(vRows.sum());
  }

  std::sort(values.begin(), values.end());

  valuesReference.clear();

  for (size_t i = 0; i < cols; i++) {
    m2->getColumn(i, vRows);
    valuesReference.push_back(vRows.sum());
  }

  std::sort(valuesReference.begin(), valuesReference.end());

  for (size_t i = 0; i < rows; i++) {
    BOOST_CHECK_SMALL(values[i] - valuesReference[i], tolerance);
  }
}

BOOST_AUTO_TEST_SUITE(TestOperationLaplaceLinear)

BOOST_AUTO_TEST_CASE(testHatRegular1D) {
  size_t dim = 1;
  size_t level = 7;

  std::unique_ptr<Grid> grid(sgpp::base::Grid::createLinearGrid(dim));
  GridStorage& storage = grid->getStorage();

  grid->getGenerator().regular(level);

  std::unique_ptr<OperationMatrix> laplace(sgpp::op_factory::createOperationLaplace(*grid));

  DataVector alpha(storage.getSize());
  DataVector result(storage.getSize());

  alpha.setAll(1.0);

  laplace->mult(alpha, result);
  HashGridPoint::index_type idx;
  HashGridPoint::level_type lvl;

  for (size_t seq = 0; seq < storage.getSize(); seq++) {
    storage.getPoint(seq).get(0, lvl, idx);
    BOOST_CHECK_CLOSE(result[seq], pow(2.0, static_cast<double>(lvl + 1)), 0.0);
  }
}

BOOST_AUTO_TEST_CASE(testHatRegulardD) {
  size_t dim = 3;
  size_t level = 3;
  std::unique_ptr<Grid> grid(sgpp::base::Grid::createLinearGrid(dim));

  std::unique_ptr<DataMatrix> m(generateLaplaceMatrix(*grid, level));

  std::string fileName("misc/tests/data/C_laplace_phi_li_hut_dim_3_nopsgrid_31_float.dat.gz");
  std::unique_ptr<DataMatrix> mRef(readReferenceMatrix(grid->getStorage(), fileName));

  compareStiffnessMatrices(m.get(), mRef.get());
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(TestOperationLaplaceEnhancedLinear)

BOOST_AUTO_TEST_CASE(testHatRegular1D) {
  size_t dim = 1;
  size_t level = 7;

  std::unique_ptr<Grid> grid(sgpp::base::Grid::createLinearGrid(dim));
  GridStorage& storage = grid->getStorage();

  grid->getGenerator().regular(level);

  std::unique_ptr<OperationMatrix> laplace(sgpp::op_factory::createOperationLaplaceEnhanced(*grid));

  DataVector alpha(storage.getSize());
  DataVector result(storage.getSize());

  alpha.setAll(1.0);

  laplace->mult(alpha, result);
  HashGridPoint::index_type idx;
  HashGridPoint::level_type lvl;

  for (size_t seq = 0; seq < storage.getSize(); seq++) {
    storage.getPoint(seq).get(0, lvl, idx);
    BOOST_CHECK_CLOSE(result[seq], pow(2.0, static_cast<double>(lvl + 1)), 0.0);
  }
}

BOOST_AUTO_TEST_CASE(testHatRegulardD) {
  size_t dim = 3;
  size_t level = 3;
  std::unique_ptr<Grid> grid(sgpp::base::Grid::createLinearGrid(dim));

  std::unique_ptr<DataMatrix> m(generateLaplaceEnhancedMatrix(*grid, level));

  std::string fileName("misc/tests/data/C_laplace_phi_li_hut_dim_3_nopsgrid_31_float.dat.gz");
  std::unique_ptr<DataMatrix> mRef(readReferenceMatrix(grid->getStorage(), fileName));

  compareStiffnessMatrices(m.get(), mRef.get());
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(TestOperationLaplaceModLinear)

BOOST_AUTO_TEST_CASE(testHatRegular1D) {
  size_t dim = 1;
  size_t level = 5;

  std::unique_ptr<Grid> grid(sgpp::base::Grid::createModLinearGrid(dim));

  std::unique_ptr<DataMatrix> m(generateLaplaceMatrix(*grid, level));

  std::string fileName(
      "misc/tests/data/C_laplace_phi_li_ausgeklappt_dim_1_nopsgrid_31_float.dat.gz");
  std::unique_ptr<DataMatrix> mRef(readReferenceMatrix(grid->getStorage(), fileName));

  compareStiffnessMatrices(m.get(), mRef.get());
}

BOOST_AUTO_TEST_CASE(testHatRegulardD) {
  size_t dim = 3;
  size_t level = 3;
  std::unique_ptr<Grid> grid(sgpp::base::Grid::createModLinearGrid(dim));

  std::unique_ptr<DataMatrix> m(generateLaplaceMatrix(*grid, level));

  std::string fileName(
      "misc/tests/data/C_laplace_phi_li_ausgeklappt_dim_3_nopsgrid_31_float.dat.gz");
  std::unique_ptr<DataMatrix> mRef(readReferenceMatrix(grid->getStorage(), fileName));

  compareStiffnessMatrices(m.get(), mRef.get());
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(TestOperationLaplaceLinearTruncatedBoundary)

BOOST_AUTO_TEST_CASE(testHatRegular1D_one) {
  size_t dim = 1;
  size_t level = 4;
  std::unique_ptr<Grid> grid(sgpp::base::Grid::createLinearBoundaryGrid(dim));

  std::unique_ptr<DataMatrix> m(generateLaplaceMatrix(*grid, level));

  std::string fileName(
      "misc/tests/data/C_laplace_phi_li_hut_trapezrand_dim_1_nopsgrid_17_float.dat.gz");
  std::unique_ptr<DataMatrix> mRef(readReferenceMatrix(grid->getStorage(), fileName));

  compareStiffnessMatrices(m.get(), mRef.get());
}

BOOST_AUTO_TEST_CASE(testHatRegular1D_two) {
  size_t dim = 1;
  size_t level = 5;
  std::unique_ptr<Grid> grid(sgpp::base::Grid::createLinearBoundaryGrid(dim));

  std::unique_ptr<DataMatrix> m(generateLaplaceMatrix(*grid, level));

  std::string fileName(
      "misc/tests/data/C_laplace_phi_li_hut_trapezrand_dim_1_nopsgrid_33_float.dat.gz");
  std::unique_ptr<DataMatrix> mRef(readReferenceMatrix(grid->getStorage(), fileName));

  compareStiffnessMatrices(m.get(), mRef.get());
}

BOOST_AUTO_TEST_CASE(testHatRegulardD_one) {
  size_t dim = 3;
  size_t level = 3;
  std::unique_ptr<Grid> grid(sgpp::base::Grid::createLinearBoundaryGrid(dim));

  std::unique_ptr<DataMatrix> m(generateLaplaceMatrix(*grid, level));

  std::string fileName(
      "misc/tests/data/C_laplace_phi_li_hut_trapezrand_dim_3_nopsgrid_225_float.dat.gz");
  std::unique_ptr<DataMatrix> mRef(readReferenceMatrix(grid->getStorage(), fileName));

  compareStiffnessMatrices(m.get(), mRef.get());
}

BOOST_AUTO_TEST_CASE(testHatRegulardD_two) {
  size_t dim = 3;
  size_t level = 2;
  std::unique_ptr<Grid> grid(sgpp::base::Grid::createLinearBoundaryGrid(dim));

  std::unique_ptr<DataMatrix> m(generateLaplaceMatrix(*grid, level));

  std::string fileName(
      "misc/tests/data/C_laplace_phi_li_hut_trapezrand_dim_3_nopsgrid_81_float.dat.gz");
  std::unique_ptr<DataMatrix> mRef(readReferenceMatrix(grid->getStorage(), fileName));

  compareStiffnessMatrices(m.get(), mRef.get());
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(TestOperationLaplaceEnhancedLinearTruncatedBoundary)

BOOST_AUTO_TEST_CASE(testHatRegular1D_one) {
  size_t dim = 1;
  size_t level = 4;
  std::unique_ptr<Grid> grid(sgpp::base::Grid::createLinearBoundaryGrid(dim));

  std::unique_ptr<DataMatrix> m(generateLaplaceEnhancedMatrix(*grid, level));

  std::string fileName(
      "misc/tests/data/C_laplace_phi_li_hut_trapezrand_dim_1_nopsgrid_17_float.dat.gz");
  std::unique_ptr<DataMatrix> mRef(readReferenceMatrix(grid->getStorage(), fileName));

  compareStiffnessMatrices(m.get(), mRef.get());
}

BOOST_AUTO_TEST_CASE(testHatRegular1D_two) {
  size_t dim = 1;
  size_t level = 5;
  std::unique_ptr<Grid> grid(sgpp::base::Grid::createLinearBoundaryGrid(dim));

  std::unique_ptr<DataMatrix> m(generateLaplaceEnhancedMatrix(*grid, level));

  std::string fileName(
      "misc/tests/data/C_laplace_phi_li_hut_trapezrand_dim_1_nopsgrid_33_float.dat.gz");
  std::unique_ptr<DataMatrix> mRef(readReferenceMatrix(grid->getStorage(), fileName));

  compareStiffnessMatrices(m.get(), mRef.get());
}

BOOST_AUTO_TEST_CASE(testHatRegulardD_one) {
  size_t dim = 3;
  size_t level = 3;
  std::unique_ptr<Grid> grid(sgpp::base::Grid::createLinearBoundaryGrid(dim));

  std::unique_ptr<DataMatrix> m(generateLaplaceEnhancedMatrix(*grid, level));

  std::string fileName(
      "misc/tests/data/C_laplace_phi_li_hut_trapezrand_dim_3_nopsgrid_225_float.dat.gz");
  std::unique_ptr<DataMatrix> mRef(readReferenceMatrix(grid->getStorage(), fileName));

  compareStiffnessMatrices(m.get(), mRef.get());
}

BOOST_AUTO_TEST_CASE(testHatRegulardD_two) {
  size_t dim = 3;
  size_t level = 2;
  std::unique_ptr<Grid> grid(sgpp::base::Grid::createLinearBoundaryGrid(dim));

  std::unique_ptr<DataMatrix> m(generateLaplaceEnhancedMatrix(*grid, level));

  std::string fileName(
      "misc/tests/data/C_laplace_phi_li_hut_trapezrand_dim_3_nopsgrid_81_float.dat.gz");
  std::unique_ptr<DataMatrix> mRef(readReferenceMatrix(grid->getStorage(), fileName));

  compareStiffnessMatrices(m.get(), mRef.get());
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(TestOperationLaplaceLinearBoundary)

BOOST_AUTO_TEST_CASE(testHatRegular1D_one) {
  size_t dim = 1;
  size_t level = 4;
  std::unique_ptr<Grid> grid(sgpp::base::Grid::createLinearBoundaryGrid(dim, 0));

  std::unique_ptr<DataMatrix> m(generateLaplaceMatrix(*grid, level));

  std::string fileName(
      "misc/tests/data/C_laplace_phi_li_hut_l0_rand_dim_1_nopsgrid_17_float.dat.gz");
  std::unique_ptr<DataMatrix> mRef(readReferenceMatrix(grid->getStorage(), fileName));

  compareStiffnessMatrices(m.get(), mRef.get());
}

BOOST_AUTO_TEST_CASE(testHatRegular1D_two) {
  size_t dim = 1;
  size_t level = 5;
  std::unique_ptr<Grid> grid(sgpp::base::Grid::createLinearBoundaryGrid(dim, 0));

  std::unique_ptr<DataMatrix> m(generateLaplaceMatrix(*grid, level));

  std::string fileName(
      "misc/tests/data/C_laplace_phi_li_hut_l0_rand_dim_1_nopsgrid_33_float.dat.gz");
  std::unique_ptr<DataMatrix> mRef(readReferenceMatrix(grid->getStorage(), fileName));

  compareStiffnessMatrices(m.get(), mRef.get());
}

BOOST_AUTO_TEST_CASE(testHatRegulardD_one) {
  size_t dim = 3;
  size_t level = 3;
  std::unique_ptr<Grid> grid(sgpp::base::Grid::createLinearBoundaryGrid(dim, 0));

  std::unique_ptr<DataMatrix> m(generateLaplaceMatrix(*grid, level));

  std::string fileName(
      "misc/tests/data/C_laplace_phi_li_hut_l0_rand_dim_3_nopsgrid_123_float.dat.gz");
  std::unique_ptr<DataMatrix> mRef(readReferenceMatrix(grid->getStorage(), fileName));

  compareStiffnessMatrices(m.get(), mRef.get());
}

BOOST_AUTO_TEST_CASE(testHatRegulardD_two) {
  size_t dim = 3;
  size_t level = 4;
  std::unique_ptr<Grid> grid(sgpp::base::Grid::createLinearBoundaryGrid(dim, 0));

  std::unique_ptr<DataMatrix> m(generateLaplaceMatrix(*grid, level));

  std::string fileName(
      "misc/tests/data/C_laplace_phi_li_hut_l0_rand_dim_3_nopsgrid_297_float.dat.gz");
  std::unique_ptr<DataMatrix> mRef(readReferenceMatrix(grid->getStorage(), fileName));

  compareStiffnessMatrices(m.get(), mRef.get());
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(TestOperationLaplaceEnhancedLinearBoundary)

BOOST_AUTO_TEST_CASE(testHatRegular1D_one) {
  size_t dim = 1;
  size_t level = 4;
  std::unique_ptr<Grid> grid(sgpp::base::Grid::createLinearBoundaryGrid(dim, 0));

  std::unique_ptr<DataMatrix> m(generateLaplaceMatrix(*grid, level));

  std::string fileName(
      "misc/tests/data/C_laplace_phi_li_hut_l0_rand_dim_1_nopsgrid_17_float.dat.gz");
  std::unique_ptr<DataMatrix> mRef(readReferenceMatrix(grid->getStorage(), fileName));

  compareStiffnessMatrices(m.get(), mRef.get());
}

BOOST_AUTO_TEST_CASE(testHatRegular1D_two) {
  size_t dim = 1;
  size_t level = 5;
  std::unique_ptr<Grid> grid(sgpp::base::Grid::createLinearBoundaryGrid(dim, 0));

  std::unique_ptr<DataMatrix> m(generateLaplaceMatrix(*grid, level));

  std::string fileName(
      "misc/tests/data/C_laplace_phi_li_hut_l0_rand_dim_1_nopsgrid_33_float.dat.gz");
  std::unique_ptr<DataMatrix> mRef(readReferenceMatrix(grid->getStorage(), fileName));

  compareStiffnessMatrices(m.get(), mRef.get());
}

BOOST_AUTO_TEST_CASE(testHatRegulardD_one) {
  size_t dim = 3;
  size_t level = 3;
  std::unique_ptr<Grid> grid(sgpp::base::Grid::createLinearBoundaryGrid(dim, 0));

  std::unique_ptr<DataMatrix> m(generateLaplaceMatrix(*grid, level));

  std::string fileName(
      "misc/tests/data/C_laplace_phi_li_hut_l0_rand_dim_3_nopsgrid_123_float.dat.gz");
  std::unique_ptr<DataMatrix> mRef(readReferenceMatrix(grid->getStorage(), fileName));

  compareStiffnessMatrices(m.get(), mRef.get());
}

BOOST_AUTO_TEST_CASE(testHatRegulardD_two) {
  size_t dim = 3;
  size_t level = 4;
  std::unique_ptr<Grid> grid(sgpp::base::Grid::createLinearBoundaryGrid(dim, 0));

  std::unique_ptr<DataMatrix> m(generateLaplaceMatrix(*grid, level));

  std::string fileName(
      "misc/tests/data/C_laplace_phi_li_hut_l0_rand_dim_3_nopsgrid_297_float.dat.gz");
  std::unique_ptr<DataMatrix> mRef(readReferenceMatrix(grid->getStorage(), fileName));

  compareStiffnessMatrices(m.get(), mRef.get());
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(TestOperationLaplacePrewavelet)

BOOST_AUTO_TEST_CASE(testPrewavelet1D_one) {
  size_t dim = 1;
  size_t level = 4;
  std::unique_ptr<Grid> grid(sgpp::base::Grid::createPrewaveletGrid(dim));

  std::unique_ptr<DataMatrix> m(generateLaplaceMatrix(*grid, level));

  std::string fileName("misc/tests/data/C_laplace_prewavelet_dim_1_nopsgrid_15_float.dat.gz");
  std::unique_ptr<DataMatrix> mRef(readReferenceMatrix(grid->getStorage(), fileName));

  compareStiffnessMatrices(m.get(), mRef.get());
}

BOOST_AUTO_TEST_CASE(testPrewavelet1D_two) {
  size_t dim = 1;
  size_t level = 5;
  std::unique_ptr<Grid> grid(sgpp::base::Grid::createPrewaveletGrid(dim));

  std::unique_ptr<DataMatrix> m(generateLaplaceMatrix(*grid, level));

  std::string fileName("misc/tests/data/C_laplace_prewavelet_dim_1_nopsgrid_31_float.dat.gz");
  std::unique_ptr<DataMatrix> mRef(readReferenceMatrix(grid->getStorage(), fileName));

  compareStiffnessMatrices(m.get(), mRef.get());
}

BOOST_AUTO_TEST_CASE(testPrewaveletD_one) {
  size_t dim = 3;
  size_t level = 3;
  std::unique_ptr<Grid> grid(sgpp::base::Grid::createPrewaveletGrid(dim));

  std::unique_ptr<DataMatrix> m(generateLaplaceMatrix(*grid, level));

  std::string fileName("misc/tests/data/C_laplace_prewavelet_dim_3_nopsgrid_31_float.dat.gz");
  std::unique_ptr<DataMatrix> mRef(readReferenceMatrix(grid->getStorage(), fileName));

  compareStiffnessMatrices(m.get(), mRef.get());
}

BOOST_AUTO_TEST_CASE(testPrewaveletD_two) {
  size_t dim = 3;
  size_t level = 4;
  std::unique_ptr<Grid> grid(sgpp::base::Grid::createPrewaveletGrid(dim));

  std::unique_ptr<DataMatrix> m(generateLaplaceMatrix(*grid, level));

  std::string fileName("misc/tests/data/C_laplace_prewavelet_dim_3_nopsgrid_111_float.dat.gz");
  std::unique_ptr<DataMatrix> mRef(readReferenceMatrix(grid->getStorage(), fileName));

  compareStiffnessMatrices(m.get(), mRef.get());
}

BOOST_AUTO_TEST_CASE(testPrewaveletAdaptiveD_two) {
  size_t dim = 4;
  size_t level = 2;
  std::unique_ptr<Grid> grid(sgpp::base::Grid::createPrewaveletGrid(dim));
  GridGenerator& generator = grid->getGenerator();
  generator.regular(level);

  GridStorage& gridStorage = grid->getStorage();
  DataVector alpha(gridStorage.getSize());

  for (size_t i = 0; i < gridStorage.getSize(); i++) {
    alpha[i] = static_cast<double>(i + 1);
  }

  size_t refinements_num = 1;
  double threshold = 0.0;
  SurplusRefinementFunctor srf = SurplusRefinementFunctor(alpha, refinements_num, threshold);
  generator.refine(srf);

  std::unique_ptr<OperationMatrix> laplace(sgpp::op_factory::createOperationLaplace(*grid));

  DataVector alpha2(gridStorage.getSize());
  DataVector erg(gridStorage.getSize());
  DataMatrix m(gridStorage.getSize(), gridStorage.getSize());

  m.setAll(0.0);

  for (size_t i = 0; i < gridStorage.getSize(); i++) {
    alpha2.setAll(0.0);
    alpha2[i] = 1.0;
    laplace->mult(alpha2, erg);
    m.setColumn(i, erg);
  }

  std::string fileName("misc/tests/data/C_laplace_prewavelet_dim_4_nopsgrid_17_adapt_float.dat.gz");
  std::unique_ptr<DataMatrix> mRef(readReferenceMatrix(gridStorage, fileName));

  compareStiffnessMatrices(&m, mRef.get());
}
BOOST_AUTO_TEST_SUITE_END()
#endif
