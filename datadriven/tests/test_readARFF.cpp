
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
#include <sgpp/datadriven/tools/ARFFTools.hpp>

#include <string>
#include <iostream>
#include <vector>

using sgpp::base::DataMatrix;
using sgpp::base::DataVector;
using sgpp::datadriven::Dataset;
using sgpp::datadriven::ARFFTools;


BOOST_AUTO_TEST_SUITE(test_dataread_arff)

BOOST_AUTO_TEST_CASE(test_fullread_hastargets) {
  std::string fileName = "datadriven/tests/datasets/dataread/simple.arff";
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
  Dataset d = ARFFTools::readARFFFromFile(fileName, true);
  controlVector.mult(-1.0); controlVector.add(d.getTargets());
  controlMatrix.mult(-1.0); controlMatrix.add(d.getData()); controlMatrix.abs();
  BOOST_CHECK_SMALL(controlVector.l2Norm(), eps);
  BOOST_CHECK_SMALL(controlMatrix.max(), eps);
}

BOOST_AUTO_TEST_CASE(test_fullread_notargets) {
  std::string fileName = "datadriven/tests/datasets/dataread/simple.arff";
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
  Dataset d = ARFFTools::readARFFFromFile(fileName, false);
  controlMatrix.mult(-1.0); controlMatrix.add(d.getData()); controlMatrix.abs();
  BOOST_CHECK_SMALL(controlMatrix.max(), eps);
}

BOOST_AUTO_TEST_CASE(test_partialread_cutoff) {
  std::string fileName = "datadriven/tests/datasets/dataread/simple.arff";
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
  Dataset d = ARFFTools::readARFFFromFile(fileName, true, 2);
  controlVector.mult(-1.0); controlVector.add(d.getTargets());
  controlMatrix.mult(-1.0); controlMatrix.add(d.getData()); controlMatrix.abs();
  BOOST_CHECK_SMALL(controlVector.l2Norm(), eps);
  BOOST_CHECK_SMALL(controlMatrix.max(), eps);
}

BOOST_AUTO_TEST_CASE(test_partialread_columns) {
  std::string fileName = "datadriven/tests/datasets/dataread/simple.arff";
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
  Dataset d = ARFFTools::readARFFFromFile(fileName, true, -1, cols);
  controlVector.mult(-1.0); controlVector.add(d.getTargets());
  controlMatrix.mult(-1.0); controlMatrix.add(d.getData()); controlMatrix.abs();
  BOOST_CHECK_SMALL(controlVector.l2Norm(), eps);
  BOOST_CHECK_SMALL(controlMatrix.max(), eps);
}

BOOST_AUTO_TEST_CASE(test_partialread_classes) {
  std::string fileName = "datadriven/tests/datasets/dataread/simple.arff";
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
  Dataset d = ARFFTools::readARFFFromFile(fileName, true, -1,
                                          std::vector<size_t>(), cls);
  controlVector.mult(-1.0); controlVector.add(d.getTargets());
  controlMatrix.mult(-1.0); controlMatrix.add(d.getData()); controlMatrix.abs();
  BOOST_CHECK_SMALL(controlVector.l2Norm(), eps);
  BOOST_CHECK_SMALL(controlMatrix.max(), eps);
}

BOOST_AUTO_TEST_SUITE_END()
