// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifdef ZLIB

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <boost/lexical_cast.hpp>

#include <algorithm>
#include <vector>
#include <string>

#include "sgpp/base/datatypes/DataVector.hpp"
#include "sgpp/base/datatypes/DataMatrix.hpp"
#include "sgpp/datadriven/tools/ARFFTools.hpp"
#include "sgpp/base/operation/BaseOpFactory.hpp"
#include "sgpp/globaldef.hpp"
#include "test_datadrivenCommon.hpp"

using sgpp::base::DataMatrix;
using sgpp::base::DataVector;
using sgpp::base::Grid;
using sgpp::base::GridStorage;
using sgpp::base::GridGenerator;
using sgpp::base::OperationMultipleEval;

DataMatrix* generateBTMatrix(Grid& grid, DataMatrix& training) {
  GridStorage& storage = grid.getStorage();

  std::unique_ptr<OperationMultipleEval> b(
      sgpp::op_factory::createOperationMultipleEval(grid, training));

  DataVector alpha(storage.getSize());
  DataVector temp(training.getNrows());

  // create BT matrix
  DataMatrix* m = new DataMatrix(training.getNrows(), storage.getSize());

  for (size_t i = 0; i < storage.getSize(); i++) {
    temp.setAll(0.0);
    alpha.setAll(0.0);
    alpha[i] = 1.0;
    b->mult(alpha, temp);

    m->setColumn(i, temp);
  }

  return m;
}

void compareBTMatrices(DataMatrix* m1, DataMatrix* m2) {
  double tolerance = 1E-5;

  // check dimensions
  BOOST_CHECK_EQUAL(m1->getNrows(), m2->getNrows());
  BOOST_CHECK_EQUAL(m1->getNcols(), m2->getNcols());

  size_t rows = m1->getNrows();  // was n

  size_t cols = m1->getNcols();  // was m

  // check row sum
  DataVector v(cols);

  std::vector<double> values;

  for (size_t i = 0; i < rows; i++) {
    m1->getRow(i, v);
    values.push_back(v.sum());

    //    std::cout << "row[" << i << "]: ";
    //    for (size_t j = 0; j < v.getSize(); j++) {
    //      if (j > 0) {
    //        std::cout << ", ";
    //      }
    //      std::cout << "(" << v[j] << ")";
    //    }
    //    std::cout << std::endl;
  }

  //  std::cout << "------------------------------" << std::endl;

  // std::sort(values.begin(), values.end());

  std::vector<double> valuesReference;

  for (size_t i = 0; i < rows; i++) {
    m2->getRow(i, v);
    valuesReference.push_back(v.sum());

    //    std::cout << "row[" << i << "]: ";
    //    for (size_t j = 0; j < v.getSize(); j++) {
    //      if (j > 0) {
    //        std::cout << ", ";
    //      }
    //      std::cout << "(" << v[j] << ")";
    //    }
    //    std::cout << std::endl;
  }

  // std::sort(valuesReference.begin(), valuesReference.end());

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

BOOST_AUTO_TEST_SUITE(TestOperationBTModLinear)

BOOST_AUTO_TEST_CASE(testHatRegular1D_one) {
  size_t level = 3;
  std::string fileName("datadriven/tests/data/data_dim_1_nops_8_float.arff.gz");
  std::string referenceMatrixFileName(
      "datadriven/tests/data/BT_phi_li_ausgeklappt_dim_1_nopsgrid_7_float.dat.gz");
  std::string content = uncompressFile(fileName);
  sgpp::datadriven::ARFFTools arffTools;
  sgpp::datadriven::Dataset dataset = arffTools.readARFFFromString(content, false);
  DataMatrix& trainingData = dataset.getData();

  size_t dim = dataset.getDimension();

  std::unique_ptr<Grid> grid(sgpp::base::Grid::createModLinearGrid(dim));
  grid->getGenerator().regular(level);
  GridStorage& gridStorage = grid->getStorage();

  DataMatrix* m = generateBTMatrix(*grid, trainingData);

  DataMatrix* mRef = readReferenceMatrix(gridStorage, referenceMatrixFileName);

  compareBTMatrices(m, mRef);
}

BOOST_AUTO_TEST_CASE(testHatRegular1D_two) {
  size_t level = 5;
  std::string fileName("datadriven/tests/data/data_dim_1_nops_8_float.arff.gz");
  std::string referenceMatrixFileName(
      "datadriven/tests/data/BT_phi_li_ausgeklappt_dim_1_nopsgrid_31_float.dat.gz");
  std::string content = uncompressFile(fileName);
  sgpp::datadriven::ARFFTools arffTools;
  sgpp::datadriven::Dataset dataset = arffTools.readARFFFromString(content, false);
  DataMatrix& trainingData = dataset.getData();

  size_t dim = dataset.getDimension();

  std::unique_ptr<Grid> grid(sgpp::base::Grid::createModLinearGrid(dim));
  grid->getGenerator().regular(level);
  GridStorage& gridStorage = grid->getStorage();

  DataMatrix* m = generateBTMatrix(*grid, trainingData);

  DataMatrix* mRef = readReferenceMatrix(gridStorage, referenceMatrixFileName);

  compareBTMatrices(m, mRef);
}

BOOST_AUTO_TEST_CASE(testHatRegulardD_one) {
  size_t level = 3;
  std::string fileName("datadriven/tests/data/data_dim_3_nops_512_float.arff.gz");
  std::string referenceMatrixFileName(
      "datadriven/tests/data/BT_phi_li_ausgeklappt_dim_3_nopsgrid_31_float.dat.gz");
  std::string content = uncompressFile(fileName);
  sgpp::datadriven::ARFFTools arffTools;
  sgpp::datadriven::Dataset dataset = arffTools.readARFFFromString(content, false);
  DataMatrix& trainingData = dataset.getData();

  size_t dim = dataset.getDimension();

  std::unique_ptr<Grid> grid(sgpp::base::Grid::createModLinearGrid(dim));
  grid->getGenerator().regular(level);
  GridStorage& gridStorage = grid->getStorage();

  DataMatrix* m = generateBTMatrix(*grid, trainingData);

  DataMatrix* mRef = readReferenceMatrix(gridStorage, referenceMatrixFileName);

  compareBTMatrices(m, mRef);
}

BOOST_AUTO_TEST_CASE(testHatRegulardD_two) {
  size_t level = 4;
  std::string fileName("datadriven/tests/data/data_dim_3_nops_512_float.arff.gz");
  std::string referenceMatrixFileName(
      "datadriven/tests/data/BT_phi_li_ausgeklappt_dim_3_nopsgrid_111_float.dat.gz");
  std::string content = uncompressFile(fileName);
  sgpp::datadriven::ARFFTools arffTools;
  sgpp::datadriven::Dataset dataset = arffTools.readARFFFromString(content, false);
  DataMatrix& trainingData = dataset.getData();

  size_t dim = dataset.getDimension();

  std::unique_ptr<Grid> grid(sgpp::base::Grid::createModLinearGrid(dim));
  grid->getGenerator().regular(level);
  GridStorage& gridStorage = grid->getStorage();

  DataMatrix* m = generateBTMatrix(*grid, trainingData);

  DataMatrix* mRef = readReferenceMatrix(gridStorage, referenceMatrixFileName);

  compareBTMatrices(m, mRef);
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(TestOperationBTLinear)

BOOST_AUTO_TEST_CASE(testHatRegular1D_one) {
  size_t level = 3;
  std::string fileName("datadriven/tests/data/data_dim_1_nops_8_float.arff.gz");
  std::string referenceMatrixFileName(
      "datadriven/tests/data/BT_phi_li_hut_dim_1_nopsgrid_7_float.dat.gz");
  std::string content = uncompressFile(fileName);
  sgpp::datadriven::ARFFTools arffTools;
  sgpp::datadriven::Dataset dataset = arffTools.readARFFFromString(content, false);
  DataMatrix& trainingData = dataset.getData();

  size_t dim = dataset.getDimension();

  std::unique_ptr<Grid> grid(sgpp::base::Grid::createLinearGrid(dim));
  grid->getGenerator().regular(level);
  GridStorage& gridStorage = grid->getStorage();

  DataMatrix* m = generateBTMatrix(*grid, trainingData);

  DataMatrix* mRef = readReferenceMatrix(gridStorage, referenceMatrixFileName);

  compareBTMatrices(m, mRef);
}

BOOST_AUTO_TEST_CASE(testHatRegular1D_two) {
  size_t level = 5;
  std::string fileName("datadriven/tests/data/data_dim_1_nops_8_float.arff.gz");
  std::string referenceMatrixFileName(
      "datadriven/tests/data/BT_phi_li_hut_dim_1_nopsgrid_31_float.dat.gz");
  std::string content = uncompressFile(fileName);
  sgpp::datadriven::ARFFTools arffTools;
  sgpp::datadriven::Dataset dataset = arffTools.readARFFFromString(content, false);
  DataMatrix& trainingData = dataset.getData();

  size_t dim = dataset.getDimension();

  std::unique_ptr<Grid> grid(sgpp::base::Grid::createLinearGrid(dim));
  grid->getGenerator().regular(level);
  GridStorage& gridStorage = grid->getStorage();

  DataMatrix* m = generateBTMatrix(*grid, trainingData);

  DataMatrix* mRef = readReferenceMatrix(gridStorage, referenceMatrixFileName);

  compareBTMatrices(m, mRef);
}

BOOST_AUTO_TEST_CASE(testHatRegulardD_one) {
  size_t level = 3;
  std::string fileName("datadriven/tests/data/data_dim_3_nops_512_float.arff.gz");
  std::string referenceMatrixFileName(
      "datadriven/tests/data/BT_phi_li_hut_dim_3_nopsgrid_31_float.dat.gz");
  std::string content = uncompressFile(fileName);
  sgpp::datadriven::ARFFTools arffTools;
  sgpp::datadriven::Dataset dataset = arffTools.readARFFFromString(content, false);
  DataMatrix& trainingData = dataset.getData();

  size_t dim = dataset.getDimension();

  std::unique_ptr<Grid> grid(sgpp::base::Grid::createLinearGrid(dim));
  grid->getGenerator().regular(level);
  GridStorage& gridStorage = grid->getStorage();

  DataMatrix* m = generateBTMatrix(*grid, trainingData);

  DataMatrix* mRef = readReferenceMatrix(gridStorage, referenceMatrixFileName);

  compareBTMatrices(m, mRef);
}

BOOST_AUTO_TEST_CASE(testHatRegulardD_two) {
  size_t level = 4;
  std::string fileName("datadriven/tests/data/data_dim_3_nops_512_float.arff.gz");
  std::string referenceMatrixFileName(
      "datadriven/tests/data/BT_phi_li_hut_dim_3_nopsgrid_111_float.dat.gz");
  std::string content = uncompressFile(fileName);
  sgpp::datadriven::ARFFTools arffTools;
  sgpp::datadriven::Dataset dataset = arffTools.readARFFFromString(content, false);
  DataMatrix& trainingData = dataset.getData();

  size_t dim = dataset.getDimension();

  std::unique_ptr<Grid> grid(sgpp::base::Grid::createLinearGrid(dim));
  grid->getGenerator().regular(level);
  GridStorage& gridStorage = grid->getStorage();

  DataMatrix* m = generateBTMatrix(*grid, trainingData);

  DataMatrix* mRef = readReferenceMatrix(gridStorage, referenceMatrixFileName);

  compareBTMatrices(m, mRef);
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(TestOperationBTLinearBoundary)

BOOST_AUTO_TEST_CASE(testHatRegular1D_one) {
  size_t level = 3;
  std::string fileName("datadriven/tests/data/data_dim_1_nops_8_float.arff.gz");
  std::string referenceMatrixFileName(
      "datadriven/tests/data/BT_phi_li_hut_l0_rand_dim_1_nopsgrid_17_float.dat.gz");
  std::string content = uncompressFile(fileName);
  sgpp::datadriven::ARFFTools arffTools;
  sgpp::datadriven::Dataset dataset = arffTools.readARFFFromString(content, false);
  DataMatrix& trainingData = dataset.getData();

  size_t dim = dataset.getDimension();

  std::unique_ptr<Grid> grid(sgpp::base::Grid::createLinearBoundaryGrid(dim, 0));
  grid->getGenerator().regular(level);
  GridStorage& gridStorage = grid->getStorage();

  DataMatrix* m = generateBTMatrix(*grid, trainingData);

  DataMatrix* mRef = readReferenceMatrix(gridStorage, referenceMatrixFileName);

  compareBTMatrices(m, mRef);
}

BOOST_AUTO_TEST_CASE(testHatRegular1D_two) {
  size_t level = 5;
  std::string fileName("datadriven/tests/data/data_dim_1_nops_8_float.arff.gz");
  std::string referenceMatrixFileName(
      "datadriven/tests/data/BT_phi_li_hut_l0_rand_dim_1_nopsgrid_33_float.dat.gz");
  std::string content = uncompressFile(fileName);
  sgpp::datadriven::ARFFTools arffTools;
  sgpp::datadriven::Dataset dataset = arffTools.readARFFFromString(content, false);
  DataMatrix& trainingData = dataset.getData();

  size_t dim = dataset.getDimension();

  std::unique_ptr<Grid> grid(sgpp::base::Grid::createLinearBoundaryGrid(dim, 0));
  grid->getGenerator().regular(level);
  GridStorage& gridStorage = grid->getStorage();

  DataMatrix* m = generateBTMatrix(*grid, trainingData);

  DataMatrix* mRef = readReferenceMatrix(gridStorage, referenceMatrixFileName);

  compareBTMatrices(m, mRef);
}

BOOST_AUTO_TEST_CASE(testHatRegulardD_one) {
  size_t level = 3;
  std::string fileName("datadriven/tests/data/data_dim_3_nops_512_float.arff.gz");
  std::string referenceMatrixFileName(
      "datadriven/tests/data/BT_phi_li_hut_l0_rand_dim_3_nopsgrid_123_float.dat.gz");
  std::string content = uncompressFile(fileName);
  sgpp::datadriven::ARFFTools arffTools;
  sgpp::datadriven::Dataset dataset = arffTools.readARFFFromString(content, false);
  DataMatrix& trainingData = dataset.getData();

  size_t dim = dataset.getDimension();

  std::unique_ptr<Grid> grid(sgpp::base::Grid::createLinearBoundaryGrid(dim, 0));
  grid->getGenerator().regular(level);
  GridStorage& gridStorage = grid->getStorage();

  DataMatrix* m = generateBTMatrix(*grid, trainingData);

  DataMatrix* mRef = readReferenceMatrix(gridStorage, referenceMatrixFileName);

  compareBTMatrices(m, mRef);
}

BOOST_AUTO_TEST_CASE(testHatRegulardD_two) {
  size_t level = 4;
  std::string fileName("datadriven/tests/data/data_dim_3_nops_512_float.arff.gz");
  std::string referenceMatrixFileName(
      "datadriven/tests/data/BT_phi_li_hut_l0_rand_dim_3_nopsgrid_297_float.dat.gz");
  std::string content = uncompressFile(fileName);
  sgpp::datadriven::ARFFTools arffTools;
  sgpp::datadriven::Dataset dataset = arffTools.readARFFFromString(content, false);
  DataMatrix& trainingData = dataset.getData();

  size_t dim = dataset.getDimension();

  std::unique_ptr<Grid> grid(sgpp::base::Grid::createLinearBoundaryGrid(dim, 0));
  grid->getGenerator().regular(level);
  GridStorage& gridStorage = grid->getStorage();

  DataMatrix* m = generateBTMatrix(*grid, trainingData);

  DataMatrix* mRef = readReferenceMatrix(gridStorage, referenceMatrixFileName);

  compareBTMatrices(m, mRef);
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(TestOperationBTLinearTruncatedBoundary)

BOOST_AUTO_TEST_CASE(testHatRegular1D_one) {
  size_t level = 4;
  std::string fileName("datadriven/tests/data/data_dim_1_nops_8_float.arff.gz");
  std::string referenceMatrixFileName(
      "datadriven/tests/data/BT_phi_li_hut_trapezrand_dim_1_nopsgrid_17_float.dat.gz");
  std::string content = uncompressFile(fileName);
  sgpp::datadriven::ARFFTools arffTools;
  sgpp::datadriven::Dataset dataset = arffTools.readARFFFromString(content, false);
  DataMatrix& trainingData = dataset.getData();

  size_t dim = dataset.getDimension();

  std::unique_ptr<Grid> grid(sgpp::base::Grid::createLinearBoundaryGrid(dim));
  grid->getGenerator().regular(level);
  GridStorage& gridStorage = grid->getStorage();

  DataMatrix* m = generateBTMatrix(*grid, trainingData);

  DataMatrix* mRef = readReferenceMatrix(gridStorage, referenceMatrixFileName);

  compareBTMatrices(m, mRef);
}

BOOST_AUTO_TEST_CASE(testHatRegular1D_two) {
  size_t level = 5;
  std::string fileName("datadriven/tests/data/data_dim_1_nops_8_float.arff.gz");
  std::string referenceMatrixFileName(
      "datadriven/tests/data/BT_phi_li_hut_trapezrand_dim_1_nopsgrid_33_float.dat.gz");
  std::string content = uncompressFile(fileName);
  sgpp::datadriven::ARFFTools arffTools;
  sgpp::datadriven::Dataset dataset = arffTools.readARFFFromString(content, false);
  DataMatrix& trainingData = dataset.getData();

  size_t dim = dataset.getDimension();

  std::unique_ptr<Grid> grid(sgpp::base::Grid::createLinearBoundaryGrid(dim));
  grid->getGenerator().regular(level);
  GridStorage& gridStorage = grid->getStorage();

  DataMatrix* m = generateBTMatrix(*grid, trainingData);

  DataMatrix* mRef = readReferenceMatrix(gridStorage, referenceMatrixFileName);

  compareBTMatrices(m, mRef);
}

BOOST_AUTO_TEST_CASE(testHatRegulardD_one) {
  size_t level = 2;
  std::string fileName("datadriven/tests/data/data_dim_3_nops_512_float.arff.gz");
  std::string referenceMatrixFileName(
      "datadriven/tests/data/BT_phi_li_hut_trapezrand_dim_3_nopsgrid_81_float.dat.gz");
  std::string content = uncompressFile(fileName);
  sgpp::datadriven::ARFFTools arffTools;
  sgpp::datadriven::Dataset dataset = arffTools.readARFFFromString(content, false);
  DataMatrix& trainingData = dataset.getData();

  size_t dim = dataset.getDimension();

  std::unique_ptr<Grid> grid(sgpp::base::Grid::createLinearBoundaryGrid(dim));
  grid->getGenerator().regular(level);
  GridStorage& gridStorage = grid->getStorage();

  DataMatrix* m = generateBTMatrix(*grid, trainingData);

  DataMatrix* mRef = readReferenceMatrix(gridStorage, referenceMatrixFileName);

  compareBTMatrices(m, mRef);
}

BOOST_AUTO_TEST_CASE(testHatRegulardD_two) {
  size_t level = 3;
  std::string fileName("datadriven/tests/data/data_dim_3_nops_512_float.arff.gz");
  std::string referenceMatrixFileName(
      "datadriven/tests/data/BT_phi_li_hut_trapezrand_dim_3_nopsgrid_225_float.dat.gz");
  std::string content = uncompressFile(fileName);
  sgpp::datadriven::ARFFTools arffTools;
  sgpp::datadriven::Dataset dataset = arffTools.readARFFFromString(content, false);
  DataMatrix& trainingData = dataset.getData();

  size_t dim = dataset.getDimension();

  std::unique_ptr<Grid> grid(sgpp::base::Grid::createLinearBoundaryGrid(dim));
  grid->getGenerator().regular(level);
  GridStorage& gridStorage = grid->getStorage();

  DataMatrix* m = generateBTMatrix(*grid, trainingData);

  DataMatrix* mRef = readReferenceMatrix(gridStorage, referenceMatrixFileName);

  compareBTMatrices(m, mRef);
}

BOOST_AUTO_TEST_SUITE_END()

#endif
