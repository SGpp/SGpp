// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <zlib.h>
#include <string>
#include <algorithm>
#include <vector>

#include "sgpp/base/datatypes/DataVector.hpp"
#include "sgpp/base/datatypes/DataMatrix.hpp"
#include "sgpp/datadriven/tools/ARFFTools.hpp"
#include "sgpp/base/operation/BaseOpFactory.hpp"
#include "sgpp/datadriven/DatadrivenOpFactory.hpp"
#include "sgpp/base/grid/generation/functors/SurplusRefinementFunctor.hpp"
#include "sgpp/globaldef.hpp"
#include "test_datadrivenCommon.hpp"

using SGPP::base::DataMatrix;
using SGPP::base::DataVector;
using SGPP::base::Grid;
using SGPP::base::GridStorage;
using SGPP::base::GridGenerator;
using SGPP::base::OperationMultipleEval;

DataMatrix* generateBBTMatrix(Grid& grid, DataMatrix& training) {
  GridStorage* storage = grid.getStorage();

  OperationMultipleEval* b = SGPP::op_factory::createOperationMultipleEval(grid, training);

  DataVector alpha(storage->size());
  DataVector erg(storage->size());
  DataVector temp(training.getNrows());

  // create BT matrix
  DataMatrix* m = new DataMatrix(storage->size(), storage->size());

  for (size_t i = 0; i < storage->size(); i++) {
    temp.setAll(0.0);
    erg.setAll(0.0);
    alpha.setAll(0.0);
    alpha[i] = 1.0;
    b->mult(alpha, temp);
    b->multTranspose(temp, erg);
    m->setColumn(i, erg);
  }

  //  for (size_t i = 0; i < storage->size(); i++) {
  //    for (size_t j = 0; j< storage->size(); j++) {
  //      std::cout << "BBT[" << i << ", " << j << "] = " << m->get(i, j) << std::endl;
  //    }
  //  }

  return m;
}

void compareBBTMatrices(DataMatrix* m1, DataMatrix* m2) {
#if USE_DOUBLE_PRECISION
  double tolerance = 1E-2;
#else
  double tolerance = 1E-1;
#endif

  // check dimensions
  BOOST_CHECK_EQUAL(m1->getNrows(), m2->getNrows());
  BOOST_CHECK_EQUAL(m1->getNcols(), m2->getNcols());

  size_t rows = m1->getNrows();  // was n

  size_t cols = m1->getNcols();  // was m

  // check diagonal
  std::vector<SGPP::float_t> valuesDiag;

  for (size_t i = 0; i < rows; i++) {
    valuesDiag.push_back(m1->get(i, i));
  }

  std::sort(valuesDiag.begin(), valuesDiag.end());

  std::vector<SGPP::float_t> valuesDiagReference;

  for (size_t i = 0; i < rows; i++) {
    valuesDiagReference.push_back(m2->get(i, i));
  }

  std::sort(valuesDiagReference.begin(), valuesDiagReference.end());

  for (size_t i = 0; i < rows; i++) {
    BOOST_CHECK_CLOSE(valuesDiag[i], valuesDiagReference[i], tolerance);
  }

  // check row sum
  DataVector v(cols);

  std::vector<SGPP::float_t> values;

  for (size_t i = 0; i < rows; i++) {
    m1->getRow(i, v);
    values.push_back(v.sum());
  }

  std::sort(values.begin(), values.end());

  std::vector<SGPP::float_t> valuesReference;

  for (size_t i = 0; i < rows; i++) {
    m2->getRow(i, v);
    valuesReference.push_back(v.sum());
  }

  std::sort(valuesReference.begin(), valuesReference.end());

  for (size_t i = 0; i < rows; i++) {
    BOOST_CHECK_CLOSE(values[i], valuesReference[i], tolerance);
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
    BOOST_CHECK_CLOSE(values[i], valuesReference[i], tolerance);
  }
}

BOOST_AUTO_TEST_SUITE(TestOperationBBTModLinear)

BOOST_AUTO_TEST_CASE(testHatRegular1D_one) {
  size_t level = 3;
  std::string fileName("datadriven/tests/data/data_dim_1_nops_8_float.arff.gz");
  std::string referenceMatrixFileName(
      "datadriven/tests/data/BBT_phi_li_ausgeklappt_dim_1_nopsgrid_7_float.dat.gz");
  std::string content = uncompressFile(fileName);
  SGPP::datadriven::ARFFTools arffTools;
  SGPP::datadriven::Dataset dataset = arffTools.readARFFFromString(content);
  DataMatrix& trainingData = dataset.getData();

  //  for (size_t i = 0; i < trainingData->getNrows(); i++) {
  //    for (size_t j = 0; j < trainingData->getNcols(); j++) {
  //      std::cout << "training[" << i << "," << j << "] = " << trainingData->get(i, j) <<
  //      std::endl;
  //    }
  //  }

  size_t dim = dataset.getDimension();

  std::unique_ptr<Grid> grid = SGPP::base::Grid::createModLinearGrid(dim);
  GridGenerator* generator = grid->createGridGenerator();
  generator->regular(level);
  GridStorage* gridStorage = grid->getStorage();

  DataMatrix* m = generateBBTMatrix(*grid, trainingData);

  DataMatrix* mRef = readReferenceMatrix(gridStorage, referenceMatrixFileName);

  compareBBTMatrices(m, mRef);
}

BOOST_AUTO_TEST_CASE(testHatRegular1D_two) {
  size_t level = 5;
  std::string fileName("datadriven/tests/data/data_dim_1_nops_8_float.arff.gz");
  std::string referenceMatrixFileName(
      "datadriven/tests/data/BBT_phi_li_ausgeklappt_dim_1_nopsgrid_31_float.dat.gz");
  std::string content = uncompressFile(fileName);
  SGPP::datadriven::ARFFTools arffTools;
  SGPP::datadriven::Dataset dataset = arffTools.readARFFFromString(content);
  DataMatrix& trainingData = dataset.getData();

  size_t dim = dataset.getDimension();

  std::unique_ptr<Grid> grid = SGPP::base::Grid::createModLinearGrid(dim);
  GridGenerator* generator = grid->createGridGenerator();
  generator->regular(level);
  GridStorage* gridStorage = grid->getStorage();

  DataMatrix* m = generateBBTMatrix(*grid, trainingData);

  DataMatrix* mRef = readReferenceMatrix(gridStorage, referenceMatrixFileName);

  compareBBTMatrices(m, mRef);
}

BOOST_AUTO_TEST_CASE(testHatRegulardD_one) {
  size_t level = 3;
  std::string fileName("datadriven/tests/data/data_dim_3_nops_512_float.arff.gz");
  std::string referenceMatrixFileName(
      "datadriven/tests/data/BBT_phi_li_ausgeklappt_dim_3_nopsgrid_31_float.dat.gz");
  std::string content = uncompressFile(fileName);
  SGPP::datadriven::ARFFTools arffTools;
  SGPP::datadriven::Dataset dataset = arffTools.readARFFFromString(content);
  DataMatrix& trainingData = dataset.getData();

  size_t dim = dataset.getDimension();

  std::unique_ptr<Grid> grid = SGPP::base::Grid::createModLinearGrid(dim);
  GridGenerator* generator = grid->createGridGenerator();
  generator->regular(level);
  GridStorage* gridStorage = grid->getStorage();

  DataMatrix* m = generateBBTMatrix(*grid, trainingData);

  DataMatrix* mRef = readReferenceMatrix(gridStorage, referenceMatrixFileName);

  compareBBTMatrices(m, mRef);
}

BOOST_AUTO_TEST_CASE(testHatRegulardD_two) {
  size_t level = 4;
  std::string fileName("datadriven/tests/data/data_dim_3_nops_512_float.arff.gz");
  std::string referenceMatrixFileName(
      "datadriven/tests/data/BBT_phi_li_ausgeklappt_dim_3_nopsgrid_111_float.dat.gz");
  std::string content = uncompressFile(fileName);
  SGPP::datadriven::ARFFTools arffTools;
  SGPP::datadriven::Dataset dataset = arffTools.readARFFFromString(content);
  DataMatrix& trainingData = dataset.getData();

  size_t dim = dataset.getDimension();

  std::unique_ptr<Grid> grid = SGPP::base::Grid::createModLinearGrid(dim);
  GridGenerator* generator = grid->createGridGenerator();
  generator->regular(level);
  GridStorage* gridStorage = grid->getStorage();

  DataMatrix* m = generateBBTMatrix(*grid, trainingData);

  DataMatrix* mRef = readReferenceMatrix(gridStorage, referenceMatrixFileName);

  compareBBTMatrices(m, mRef);
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(TestOperationBBTLinear)

BOOST_AUTO_TEST_CASE(testHatRegular1D_one) {
  size_t level = 3;
  std::string fileName("datadriven/tests/data/data_dim_1_nops_8_float.arff.gz");
  std::string referenceMatrixFileName(
      "datadriven/tests/data/BBT_phi_li_hut_dim_1_nopsgrid_7_float.dat.gz");
  std::string content = uncompressFile(fileName);
  SGPP::datadriven::ARFFTools arffTools;
  SGPP::datadriven::Dataset dataset = arffTools.readARFFFromString(content);

  DataMatrix& trainingData = dataset.getData();

  size_t dim = dataset.getDimension();

  std::unique_ptr<Grid> grid = SGPP::base::Grid::createLinearGrid(dim);
  GridGenerator* generator = grid->createGridGenerator();
  generator->regular(level);
  GridStorage* gridStorage = grid->getStorage();

  DataMatrix* m = generateBBTMatrix(*grid, trainingData);

  DataMatrix* mRef = readReferenceMatrix(gridStorage, referenceMatrixFileName);

  compareBBTMatrices(m, mRef);
}

BOOST_AUTO_TEST_CASE(testHatRegular1D_two) {
  size_t level = 5;
  std::string fileName("datadriven/tests/data/data_dim_1_nops_8_float.arff.gz");
  std::string referenceMatrixFileName(
      "datadriven/tests/data/BBT_phi_li_hut_dim_1_nopsgrid_31_float.dat.gz");
  std::string content = uncompressFile(fileName);
  SGPP::datadriven::ARFFTools arffTools;
  SGPP::datadriven::Dataset dataset = arffTools.readARFFFromString(content);
  DataMatrix& trainingData = dataset.getData();

  size_t dim = dataset.getDimension();

  std::unique_ptr<Grid> grid = SGPP::base::Grid::createLinearGrid(dim);
  GridGenerator* generator = grid->createGridGenerator();
  generator->regular(level);
  GridStorage* gridStorage = grid->getStorage();

  DataMatrix* m = generateBBTMatrix(*grid, trainingData);

  DataMatrix* mRef = readReferenceMatrix(gridStorage, referenceMatrixFileName);

  compareBBTMatrices(m, mRef);
}

BOOST_AUTO_TEST_CASE(testHatRegulardD_one) {
  size_t level = 3;
  std::string fileName("datadriven/tests/data/data_dim_3_nops_512_float.arff.gz");
  std::string referenceMatrixFileName(
      "datadriven/tests/data/BBT_phi_li_hut_dim_3_nopsgrid_31_float.dat.gz");
  std::string content = uncompressFile(fileName);
  SGPP::datadriven::ARFFTools arffTools;
  SGPP::datadriven::Dataset dataset = arffTools.readARFFFromString(content);
  DataMatrix& trainingData = dataset.getData();

  size_t dim = dataset.getDimension();

  std::unique_ptr<Grid> grid = SGPP::base::Grid::createLinearGrid(dim);
  GridGenerator* generator = grid->createGridGenerator();
  generator->regular(level);
  GridStorage* gridStorage = grid->getStorage();

  DataMatrix* m = generateBBTMatrix(*grid, trainingData);

  DataMatrix* mRef = readReferenceMatrix(gridStorage, referenceMatrixFileName);

  compareBBTMatrices(m, mRef);
}

BOOST_AUTO_TEST_CASE(testHatRegulardD_two) {
  size_t level = 4;
  std::string fileName("datadriven/tests/data/data_dim_3_nops_512_float.arff.gz");
  std::string referenceMatrixFileName(
      "datadriven/tests/data/BBT_phi_li_hut_dim_3_nopsgrid_111_float.dat.gz");
  std::string content = uncompressFile(fileName);
  SGPP::datadriven::ARFFTools arffTools;
  SGPP::datadriven::Dataset dataset = arffTools.readARFFFromString(content);
  DataMatrix& trainingData = dataset.getData();

  size_t dim = dataset.getDimension();

  std::unique_ptr<Grid> grid = SGPP::base::Grid::createLinearGrid(dim);
  GridGenerator* generator = grid->createGridGenerator();
  generator->regular(level);
  GridStorage* gridStorage = grid->getStorage();

  DataMatrix* m = generateBBTMatrix(*grid, trainingData);

  DataMatrix* mRef = readReferenceMatrix(gridStorage, referenceMatrixFileName);

  compareBBTMatrices(m, mRef);
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(TestOperationBBTPrewavelet)

BOOST_AUTO_TEST_CASE(testPrewavelet1D_one) {
  size_t level = 3;
  std::string fileName("datadriven/tests/data/data_dim_1_nops_8_float.arff.gz");
  std::string referenceMatrixFileName(
      "datadriven/tests/data/BBT_prewavelet_dim_1_nopsgrid_7_float.dat.gz");
  std::string content = uncompressFile(fileName);
  SGPP::datadriven::ARFFTools arffTools;
  SGPP::datadriven::Dataset dataset = arffTools.readARFFFromString(content);
  DataMatrix& trainingData = dataset.getData();

  size_t dim = dataset.getDimension();

  std::unique_ptr<Grid> grid = SGPP::base::Grid::createPrewaveletGrid(dim);
  GridGenerator* generator = grid->createGridGenerator();
  generator->regular(level);
  GridStorage* gridStorage = grid->getStorage();

  DataMatrix* m = generateBBTMatrix(*grid, trainingData);

  DataMatrix* mRef = readReferenceMatrix(gridStorage, referenceMatrixFileName);

  compareBBTMatrices(m, mRef);
}

BOOST_AUTO_TEST_CASE(testPrewavelet1D_two) {
  size_t level = 5;
  std::string fileName("datadriven/tests/data/data_dim_1_nops_8_float.arff.gz");
  std::string referenceMatrixFileName(
      "datadriven/tests/data/BBT_prewavelet_dim_1_nopsgrid_31_float.dat.gz");
  std::string content = uncompressFile(fileName);
  SGPP::datadriven::ARFFTools arffTools;
  SGPP::datadriven::Dataset dataset = arffTools.readARFFFromString(content);
  DataMatrix& trainingData = dataset.getData();

  size_t dim = dataset.getDimension();

  std::unique_ptr<Grid> grid = SGPP::base::Grid::createPrewaveletGrid(dim);
  GridGenerator* generator = grid->createGridGenerator();
  generator->regular(level);
  GridStorage* gridStorage = grid->getStorage();

  DataMatrix* m = generateBBTMatrix(*grid, trainingData);

  DataMatrix* mRef = readReferenceMatrix(gridStorage, referenceMatrixFileName);

  compareBBTMatrices(m, mRef);
}

BOOST_AUTO_TEST_CASE(testPrewaveletdD_one) {
  size_t level = 3;
  std::string fileName("datadriven/tests/data/data_dim_3_nops_512_float.arff.gz");
  std::string referenceMatrixFileName(
      "datadriven/tests/data/BBT_prewavelet_dim_3_nopsgrid_31_float.dat.gz");
  std::string content = uncompressFile(fileName);
  SGPP::datadriven::ARFFTools arffTools;
  SGPP::datadriven::Dataset dataset = arffTools.readARFFFromString(content);
  DataMatrix& trainingData = dataset.getData();

  size_t dim = dataset.getDimension();

  std::unique_ptr<Grid> grid = SGPP::base::Grid::createPrewaveletGrid(dim);
  GridGenerator* generator = grid->createGridGenerator();
  generator->regular(level);
  GridStorage* gridStorage = grid->getStorage();

  DataMatrix* m = generateBBTMatrix(*grid, trainingData);

  DataMatrix* mRef = readReferenceMatrix(gridStorage, referenceMatrixFileName);

  compareBBTMatrices(m, mRef);
}

BOOST_AUTO_TEST_CASE(testPrewaveletdD_two) {
  size_t level = 4;
  std::string fileName("datadriven/tests/data/data_dim_3_nops_512_float.arff.gz");
  std::string referenceMatrixFileName(
      "datadriven/tests/data/BBT_prewavelet_dim_3_nopsgrid_111_float.dat.gz");
  std::string content = uncompressFile(fileName);
  SGPP::datadriven::ARFFTools arffTools;
  SGPP::datadriven::Dataset dataset = arffTools.readARFFFromString(content);
  DataMatrix& trainingData = dataset.getData();

  size_t dim = dataset.getDimension();

  std::unique_ptr<Grid> grid = SGPP::base::Grid::createPrewaveletGrid(dim);
  GridGenerator* generator = grid->createGridGenerator();
  generator->regular(level);
  GridStorage* gridStorage = grid->getStorage();

  DataMatrix* m = generateBBTMatrix(*grid, trainingData);

  DataMatrix* mRef = readReferenceMatrix(gridStorage, referenceMatrixFileName);

  compareBBTMatrices(m, mRef);
}

BOOST_AUTO_TEST_CASE(testPrewaveletAdaptivedD_two) {
  size_t level = 2;
  std::string fileName("datadriven/tests/data/data_dim_4_nops_4096_float.arff.gz");
  std::string referenceMatrixFileName(
      "datadriven/tests/data/BBT_prewavelet_dim_4_nopsgrid_17_adapt_float.dat.gz");
  std::string content = uncompressFile(fileName);
  SGPP::datadriven::ARFFTools arffTools;
  SGPP::datadriven::Dataset dataset = arffTools.readARFFFromString(content);
  DataMatrix& trainingData = dataset.getData();

  size_t dim = dataset.getDimension();

  std::unique_ptr<Grid> grid = SGPP::base::Grid::createPrewaveletGrid(dim);
  GridGenerator* generator = grid->createGridGenerator();
  generator->regular(level);
  GridStorage* gridStorage = grid->getStorage();

  DataVector alpha(gridStorage->size());

  for (size_t i = 0; i < gridStorage->size(); i++) {
    alpha[i] = static_cast<double>(i + 1);
  }

  SGPP::base::SurplusRefinementFunctor functor(&alpha, 1);
  generator->refine(&functor);

  DataMatrix* m = generateBBTMatrix(*grid, trainingData);

  DataMatrix* mRef = readReferenceMatrix(gridStorage, referenceMatrixFileName);

  compareBBTMatrices(m, mRef);
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(TestOperationBBTLinearBoundary)

BOOST_AUTO_TEST_CASE(testHatRegular1D_one) {
  size_t level = 4;
  std::string fileName("datadriven/tests/data/data_dim_1_nops_8_float.arff.gz");
  std::string referenceMatrixFileName(
      "datadriven/tests/data/BBT_phi_li_hut_l0_rand_dim_1_nopsgrid_17_float.dat.gz");
  std::string content = uncompressFile(fileName);
  SGPP::datadriven::ARFFTools arffTools;
  SGPP::datadriven::Dataset dataset = arffTools.readARFFFromString(content);
  DataMatrix& trainingData = dataset.getData();

  size_t dim = dataset.getDimension();

  std::unique_ptr<Grid> grid = SGPP::base::Grid::createLinearBoundaryGrid(dim, 0);
  GridGenerator* generator = grid->createGridGenerator();
  generator->regular(level);
  GridStorage* gridStorage = grid->getStorage();

  DataMatrix* m = generateBBTMatrix(*grid, trainingData);

  DataMatrix* mRef = readReferenceMatrix(gridStorage, referenceMatrixFileName);

  compareBBTMatrices(m, mRef);
}

BOOST_AUTO_TEST_CASE(testHatRegulardD_one) {
  size_t level = 3;
  std::string fileName("datadriven/tests/data/data_dim_3_nops_512_float.arff.gz");
  std::string referenceMatrixFileName(
      "datadriven/tests/data/BBT_phi_li_hut_l0_rand_dim_3_nopsgrid_123_float.dat.gz");
  std::string content = uncompressFile(fileName);
  SGPP::datadriven::ARFFTools arffTools;
  SGPP::datadriven::Dataset dataset = arffTools.readARFFFromString(content);
  DataMatrix& trainingData = dataset.getData();

  size_t dim = dataset.getDimension();

  std::unique_ptr<Grid> grid = SGPP::base::Grid::createLinearBoundaryGrid(dim, 0);
  GridGenerator* generator = grid->createGridGenerator();
  generator->regular(level);
  GridStorage* gridStorage = grid->getStorage();

  DataMatrix* m = generateBBTMatrix(*grid, trainingData);

  DataMatrix* mRef = readReferenceMatrix(gridStorage, referenceMatrixFileName);

  compareBBTMatrices(m, mRef);
}

BOOST_AUTO_TEST_CASE(testHatRegulardD_two) {
  size_t level = 4;
  std::string fileName("datadriven/tests/data/data_dim_3_nops_512_float.arff.gz");
  std::string referenceMatrixFileName(
      "datadriven/tests/data/BBT_phi_li_hut_l0_rand_dim_3_nopsgrid_297_float.dat.gz");
  std::string content = uncompressFile(fileName);
  SGPP::datadriven::ARFFTools arffTools;
  SGPP::datadriven::Dataset dataset = arffTools.readARFFFromString(content);
  DataMatrix& trainingData = dataset.getData();

  size_t dim = dataset.getDimension();

  std::unique_ptr<Grid> grid = SGPP::base::Grid::createLinearBoundaryGrid(dim, 0);
  GridGenerator* generator = grid->createGridGenerator();
  generator->regular(level);
  GridStorage* gridStorage = grid->getStorage();

  DataMatrix* m = generateBBTMatrix(*grid, trainingData);

  DataMatrix* mRef = readReferenceMatrix(gridStorage, referenceMatrixFileName);

  compareBBTMatrices(m, mRef);
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(TestOperationBBTLinearTruncatedBoundary)

BOOST_AUTO_TEST_CASE(testHatRegular1D_one) {
  size_t level = 4;
  std::string fileName("datadriven/tests/data/data_dim_1_nops_8_float.arff.gz");
  std::string referenceMatrixFileName(
      "datadriven/tests/data/BBT_phi_li_hut_trapezrand_dim_1_nopsgrid_17_float.dat.gz");
  std::string content = uncompressFile(fileName);
  SGPP::datadriven::ARFFTools arffTools;
  SGPP::datadriven::Dataset dataset = arffTools.readARFFFromString(content);
  DataMatrix& trainingData = dataset.getData();

  size_t dim = dataset.getDimension();

  std::unique_ptr<Grid> grid = SGPP::base::Grid::createLinearBoundaryGrid(dim);
  GridGenerator* generator = grid->createGridGenerator();
  generator->regular(level);
  GridStorage* gridStorage = grid->getStorage();

  DataMatrix* m = generateBBTMatrix(*grid, trainingData);

  DataMatrix* mRef = readReferenceMatrix(gridStorage, referenceMatrixFileName);

  compareBBTMatrices(m, mRef);
}

BOOST_AUTO_TEST_CASE(testHatRegular1D_two) {
  size_t level = 5;
  std::string fileName("datadriven/tests/data/data_dim_1_nops_8_float.arff.gz");
  std::string referenceMatrixFileName(
      "datadriven/tests/data/BBT_phi_li_hut_trapezrand_dim_1_nopsgrid_33_float.dat.gz");
  std::string content = uncompressFile(fileName);
  SGPP::datadriven::ARFFTools arffTools;
  SGPP::datadriven::Dataset dataset = arffTools.readARFFFromString(content);
  DataMatrix& trainingData = dataset.getData();

  size_t dim = dataset.getDimension();

  std::unique_ptr<Grid> grid = SGPP::base::Grid::createLinearBoundaryGrid(dim);
  GridGenerator* generator = grid->createGridGenerator();
  generator->regular(level);
  GridStorage* gridStorage = grid->getStorage();

  DataMatrix* m = generateBBTMatrix(*grid, trainingData);

  DataMatrix* mRef = readReferenceMatrix(gridStorage, referenceMatrixFileName);

  compareBBTMatrices(m, mRef);
}

BOOST_AUTO_TEST_CASE(testHatRegulardD_one) {
  size_t level = 2;
  std::string fileName("datadriven/tests/data/data_dim_3_nops_512_float.arff.gz");
  std::string referenceMatrixFileName(
      "datadriven/tests/data/BBT_phi_li_hut_trapezrand_dim_3_nopsgrid_81_float.dat.gz");
  std::string content = uncompressFile(fileName);
  SGPP::datadriven::ARFFTools arffTools;
  SGPP::datadriven::Dataset dataset = arffTools.readARFFFromString(content);
  DataMatrix& trainingData = dataset.getData();

  size_t dim = dataset.getDimension();

  std::unique_ptr<Grid> grid = SGPP::base::Grid::createLinearBoundaryGrid(dim);
  GridGenerator* generator = grid->createGridGenerator();
  generator->regular(level);
  GridStorage* gridStorage = grid->getStorage();

  DataMatrix* m = generateBBTMatrix(*grid, trainingData);

  DataMatrix* mRef = readReferenceMatrix(gridStorage, referenceMatrixFileName);

  compareBBTMatrices(m, mRef);
}

BOOST_AUTO_TEST_CASE(testHatRegulardD_two) {
  size_t level = 3;
  std::string fileName("datadriven/tests/data/data_dim_3_nops_512_float.arff.gz");
  std::string referenceMatrixFileName(
      "datadriven/tests/data/BBT_phi_li_hut_trapezrand_dim_3_nopsgrid_225_float.dat.gz");
  std::string content = uncompressFile(fileName);
  SGPP::datadriven::ARFFTools arffTools;
  SGPP::datadriven::Dataset dataset = arffTools.readARFFFromString(content);
  DataMatrix& trainingData = dataset.getData();

  size_t dim = dataset.getDimension();

  std::unique_ptr<Grid> grid = SGPP::base::Grid::createLinearBoundaryGrid(dim);
  GridGenerator* generator = grid->createGridGenerator();
  generator->regular(level);
  GridStorage* gridStorage = grid->getStorage();

  DataMatrix* m = generateBBTMatrix(*grid, trainingData);

  DataMatrix* mRef = readReferenceMatrix(gridStorage, referenceMatrixFileName);

  compareBBTMatrices(m, mRef);
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(TestLinearGrid)

BOOST_AUTO_TEST_CASE(testOperationTest_test) {
  size_t level = 1;
  std::unique_ptr<Grid> grid = SGPP::base::Grid::createLinearGrid(1);
  GridGenerator* generator = grid->createGridGenerator();
  generator->regular(level);
  GridStorage* gridStorage = grid->getStorage();

  DataVector alpha(gridStorage->size());
  DataMatrix data(1, 1);
  data.setAll(0.25);
  DataVector classes(1);
  classes.setAll(1.0);

  SGPP::datadriven::OperationTest* testOP = SGPP::op_factory::createOperationTest(*grid);

  alpha.setAll(1.0);
  SGPP::float_t c = testOP->test(alpha, data, classes);
  BOOST_CHECK(c > 0.0);

  alpha.setAll(-1.0);
  c = testOP->test(alpha, data, classes);
  BOOST_CHECK(c == 0.0);
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(TestLinearBoundaryGrid)

BOOST_AUTO_TEST_CASE(testOperationTest_test) {
  size_t level = 1;
  std::unique_ptr<Grid> grid = SGPP::base::Grid::createLinearBoundaryGrid(1);
  GridGenerator* generator = grid->createGridGenerator();
  generator->regular(level);
  GridStorage* gridStorage = grid->getStorage();

  DataVector alpha(gridStorage->size());
  DataMatrix data(1, 1);
  data.setAll(0.25);
  DataVector classes(1);
  classes.setAll(1.0);

  SGPP::datadriven::OperationTest* testOP = SGPP::op_factory::createOperationTest(*grid);

  alpha[0] = 0.0;
  alpha[1] = 0.0;
  alpha[2] = 1.0;
  SGPP::float_t c = testOP->test(alpha, data, classes);
  BOOST_CHECK(c > 0.0);

  alpha[0] = 0.0;
  alpha[1] = 0.0;
  alpha[2] = -1.0;
  c = testOP->test(alpha, data, classes);
  BOOST_CHECK(c == 0.0);
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(TestLinearL0BoundaryGrid)

BOOST_AUTO_TEST_CASE(testOperationTest_test) {
  size_t level = 1;
  std::unique_ptr<Grid> grid = SGPP::base::Grid::createLinearBoundaryGrid(1, 0);
  GridGenerator* generator = grid->createGridGenerator();
  generator->regular(level);
  GridStorage* gridStorage = grid->getStorage();

  DataVector alpha(gridStorage->size());
  DataMatrix data(1, 1);
  data.setAll(0.25);
  DataVector classes(1);
  classes.setAll(1.0);

  SGPP::datadriven::OperationTest* testOP = SGPP::op_factory::createOperationTest(*grid);

  alpha[0] = 0.0;
  alpha[1] = 0.0;
  alpha[2] = 1.0;
  SGPP::float_t c = testOP->test(alpha, data, classes);
  BOOST_CHECK(c > 0.0);

  alpha[0] = 0.0;
  alpha[1] = 0.0;
  alpha[2] = -1.0;
  c = testOP->test(alpha, data, classes);
  BOOST_CHECK(c == 0.0);
}

BOOST_AUTO_TEST_SUITE_END()
