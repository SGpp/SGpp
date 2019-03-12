
/* Copyright (C) 2018-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * testCSVTools.cpp
 *
 *  Created on: 01.09.2018
 *      Author: Sebastian Kreisel
 */

#define BOOST_TEST_DYN_LINK

#include <boost/test/unit_test.hpp>
#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>

#include <sgpp/datadriven/tools/Dataset.hpp>
#include <sgpp/datadriven/tools/CSVTools.hpp>

#include <string>
#include <iostream>
#include <vector>

using sgpp::base::DataMatrix;
using sgpp::base::DataVector;
using sgpp::datadriven::Dataset;
using sgpp::datadriven::CSVTools;


BOOST_AUTO_TEST_SUITE(test_dataread_csv)

BOOST_AUTO_TEST_CASE(test_fullread_hastargets) {
  std::string fileName = "datadriven/datasets/dataread/simple.csv";
  double eps = 1e-06;
  double data[] = {
    0.0, 0.1, 0.2, 0.3, 0.4,
    -1, 0, 1, 2, 3,
    7e-01, 7e-02, -7e-01, 7, -7,
    8.0, 8.0, 8.0, 8.0, 8.0,
    9, 9, 9, 9, 9
  };
  double targets[] = {
    42.42, -5, 7.0, 7, 7e+00
  };
  DataMatrix controlMatrix(data, 5, 5);
  DataVector controlVector(targets, 5);
  Dataset d = CSVTools::readCSVFromFile(fileName, true, true);
  controlVector.mult(-1.0); controlVector.add(d.getTargets());
  controlMatrix.mult(-1.0); controlMatrix.add(d.getData()); controlMatrix.abs();
  BOOST_CHECK_SMALL(controlVector.l2Norm(), eps);
  BOOST_CHECK_SMALL(controlMatrix.max(), eps);
}

BOOST_AUTO_TEST_CASE(test_fullread_notargets) {
  std::string fileName = "datadriven/datasets/dataread/simple.csv";
  double eps = 1e-06;
  double data[] = {
    0.0, 0.1, 0.2, 0.3, 0.4,
    -1, 0, 1, 2, 3,
    7e-01, 7e-02, -7e-01, 7, -7,
    8.0, 8.0, 8.0, 8.0, 8.0,
    9, 9, 9, 9, 9
  };
  double targets[] = {
    42.42, -5, 7.0, 7, 7e+00
  };
  DataMatrix controlMatrix(data, 5, 5);
  DataVector controlVector(targets, 5);
  controlMatrix.appendCol(controlVector);
  Dataset d = CSVTools::readCSVFromFile(fileName, true, false);
  controlMatrix.mult(-1.0); controlMatrix.add(d.getData()); controlMatrix.abs();
  BOOST_CHECK_SMALL(controlMatrix.max(), eps);
}

BOOST_AUTO_TEST_CASE(test_partialread_cutoff) {
  std::string fileName = "datadriven/datasets/dataread/simple.csv";
  double eps = 1e-06;
  double data[] = {
    0.0, 0.1, 0.2, 0.3, 0.4,
    -1, 0, 1, 2, 3,
  };
  double targets[] = {
    42.42, -5
  };
  DataMatrix controlMatrix(data, 2, 5);
  DataVector controlVector(targets, 2);
  Dataset d = CSVTools::readCSVFromFile(fileName, true, true, 2);
  controlVector.mult(-1.0); controlVector.add(d.getTargets());
  controlMatrix.mult(-1.0); controlMatrix.add(d.getData()); controlMatrix.abs();
  BOOST_CHECK_SMALL(controlVector.l2Norm(), eps);
  BOOST_CHECK_SMALL(controlMatrix.max(), eps);
}

BOOST_AUTO_TEST_CASE(test_partialread_columns) {
  std::string fileName = "datadriven/datasets/dataread/simple.csv";
  double eps = 1e-06;
  // col 4 2 0
  double data[] = {
    0.4, 0.2, 0.0,
    3, 1, -1,
    -7, -7e-01, 7e-01,
    8.0, 8.0, 8.0,
    9, 9, 9
  };
  double targets[] = {
    42.42, -5, 7.0, 7, 7e+00
  };
  DataMatrix controlMatrix(data, 5, 3);
  DataVector controlVector(targets, 5);
  std::vector<size_t> cols;
  cols.push_back(4); cols.push_back(2); cols.push_back(0);
  Dataset d = CSVTools::readCSVFromFile(fileName, true, true, -1, cols);
  controlVector.mult(-1.0); controlVector.add(d.getTargets());
  controlMatrix.mult(-1.0); controlMatrix.add(d.getData()); controlMatrix.abs();
  BOOST_CHECK_SMALL(controlVector.l2Norm(), eps);
  BOOST_CHECK_SMALL(controlMatrix.max(), eps);
}

BOOST_AUTO_TEST_CASE(test_partialread_classes) {
  std::string fileName = "datadriven/datasets/dataread/simple.csv";
  double eps = 1e-06;
  double data[] = {
    0.0, 0.1, 0.2, 0.3, 0.4,
    7e-01, 7e-02, -7e-01, 7, -7,
    8.0, 8.0, 8.0, 8.0, 8.0,
    9, 9, 9, 9, 9
  };
  double targets[] = {
    42.42, 7.0, 7, 7e+00
  };
  DataMatrix controlMatrix(data, 4, 5);
  DataVector controlVector(targets, 4);
  std::vector<double> cls;
  cls.push_back(7.0); cls.push_back(42.42);
  Dataset d = CSVTools::readCSVFromFile(fileName, true, true, -1,
                                        std::vector<size_t>(), cls);
  controlVector.mult(-1.0); controlVector.add(d.getTargets());
  controlMatrix.mult(-1.0); controlMatrix.add(d.getData()); controlMatrix.abs();
  BOOST_CHECK_SMALL(controlVector.l2Norm(), eps);
  BOOST_CHECK_SMALL(controlMatrix.max(), eps);
}

// (sebastian) old test by Eric Koepke, Michael Lettrich
BOOST_AUTO_TEST_CASE(test_read_csv_old) {
  double testPoints[10][3] = {{0.307143, 0.130137, 0.050000}, {0.365584, 0.105479, 0.050000},
                              {0.178571, 0.201027, 0.050000}, {0.272078, 0.145548, 0.050000},
                              {0.318831, 0.065411, 0.050000}, {0.190260, 0.086986, 0.050000},
                              {0.190260, 0.062329, 0.072500}, {0.120130, 0.068493, 0.072500},
                              {0.225325, 0.056164, 0.072500}, {0.213636, 0.050000, 0.072500}};

  double testValues[10] = {-1., 1., 1., 1., 1., 1., -1., -1., -1., -1.};

  sgpp::datadriven::Dataset dataSet =
    sgpp::datadriven::CSVTools::readCSVFromFile(
      "datadriven/datasets/liver/liver-disorders_normalized_small.csv", true);

  sgpp::base::DataVector& classes = dataSet.getTargets();
  sgpp::base::DataMatrix& data = dataSet.getData();
  size_t nrows = data.getNrows();
  size_t ncols = data.getNcols();

  double tolerance = 1E-5;

  // Check if all dimensions agree
  BOOST_CHECK_EQUAL(10, classes.getSize());
  BOOST_CHECK_EQUAL(10, data.getNrows());
  BOOST_CHECK_EQUAL(3, data.getNcols());

  sgpp::base::DataVector testVector = sgpp::base::DataVector(ncols);
  for (size_t rowIdx = 0; rowIdx < nrows; rowIdx++) {
    data.getRow(rowIdx, testVector);
    for (size_t colIdx = 0; colIdx < ncols; colIdx++) {
      // this only works because dataset does not contain any zeros. if the dataset to test with
      // changed, we need a else if here and a boost check small.
      BOOST_CHECK_CLOSE(data.get(rowIdx, colIdx), testPoints[rowIdx][colIdx], tolerance);
    }
    // this only works because dataset does not contain any zeros. if the dataset to test with
    // changed, we need a else if here and a boost check small.
    BOOST_CHECK_CLOSE(classes.get(rowIdx), testValues[rowIdx], tolerance);
  }
}

BOOST_AUTO_TEST_SUITE_END()
