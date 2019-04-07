/* Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * datamingArffSampleProviderTest.cpp
 *
 *  Created on: 01.04.2016
 *      Author: Michael Lettrich
 */

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/datadriven/datamining/modules/dataSource/ArffFileSampleProvider.hpp>
#include <sgpp/datadriven/tools/Dataset.hpp>
#include <sgpp/globaldef.hpp>

#include <string>

using sgpp::datadriven::ArffFileSampleProvider;
using sgpp::base::DataVector;
using sgpp::base::DataMatrix;
using sgpp::datadriven::Dataset;

BOOST_AUTO_TEST_SUITE(datamingArffSampleProviderTest)

const std::string datasetPath = "datadriven/datasets/liver/liver-disorders_normalized_small.arff";

const double testPoints[10][3] = {{0.307143, 0.130137, 0.050000}, {0.365584, 0.105479, 0.050000},
                                  {0.178571, 0.201027, 0.050000}, {0.272078, 0.145548, 0.050000},
                                  {0.318831, 0.065411, 0.050000}, {0.190260, 0.086986, 0.050000},
                                  {0.190260, 0.062329, 0.072500}, {0.120130, 0.068493, 0.072500},
                                  {0.225325, 0.056164, 0.072500}, {0.213636, 0.050000, 0.072500}};

const double testValues[10] = {-1., 1., 1., 1., 1., 1., -1., -1., -1., -1.};
const size_t datasetDim = 3;
const size_t datasetSize = 10;
const double tolerance = 1E-5;

BOOST_AUTO_TEST_CASE(arffTestReadFile) {
  auto sampleProvider = ArffFileSampleProvider();
  sampleProvider.readFile(datasetPath, true);
  auto dataset = sampleProvider.getAllSamples();

  DataVector& classes = dataset->getTargets();
  DataMatrix& data = dataset->getData();
  size_t nrows = data.getNrows();
  size_t ncols = data.getNcols();

  // Check if all dimensions agree
  BOOST_CHECK_EQUAL(datasetSize, classes.getSize());
  BOOST_CHECK_EQUAL(datasetSize, data.getNrows());
  BOOST_CHECK_EQUAL(datasetDim, data.getNcols());

  DataVector testVector = DataVector(ncols);
  for (size_t rowIdx = 0; rowIdx < nrows; rowIdx++) {
    data.getRow(rowIdx, testVector);
    for (size_t colIdx = 0; colIdx < ncols; colIdx++) {
      // this only works because dataset does not contain any zeros. if the dataset to test with
      // changed, we need a else if here and a boost check small.
      BOOST_CHECK_CLOSE(data.get(rowIdx, colIdx), testPoints[rowIdx][colIdx], tolerance);
    }
    // this only works because dataset does not contain any zeros. if the dataset to test with is
    // changed, we need a else if here and a boost check small.
    BOOST_CHECK_CLOSE(classes.get(rowIdx), testValues[rowIdx], tolerance);
  }
}

BOOST_AUTO_TEST_CASE(arffTestgetSize) {
  auto sampleProvider = ArffFileSampleProvider();
  sampleProvider.readFile(datasetPath, true);
  BOOST_CHECK_EQUAL(datasetSize, sampleProvider.getNumSamples());
}

BOOST_AUTO_TEST_CASE(arffTestGetNextSamples) {
  size_t sampleSize1 = 5;
  size_t sampleSize2 = 3;

  auto sampleProvider = ArffFileSampleProvider();
  sampleProvider.readFile(datasetPath, true);

  auto dataset = std::unique_ptr<Dataset>(sampleProvider.getNextSamples(sampleSize1));
  DataVector& classesSample1 = dataset->getTargets();
  DataMatrix& dataSample1 = dataset->getData();
  size_t nrows = dataSample1.getNrows();
  size_t ncols = dataSample1.getNcols();

  // Check if all dimensions agree
  BOOST_CHECK_EQUAL(sampleSize1, classesSample1.getSize());
  BOOST_CHECK_EQUAL(sampleSize1, dataSample1.getNrows());
  BOOST_CHECK_EQUAL(datasetDim, dataSample1.getNcols());

  DataVector testVector = DataVector(ncols);
  for (size_t rowIdx = 0; rowIdx < nrows; rowIdx++) {
    dataSample1.getRow(rowIdx, testVector);
    for (size_t colIdx = 0; colIdx < ncols; colIdx++) {
      // this only works because dataset does not contain any zeros. if the dataset to test with
      // changed, we need a else if here and a boost check small.
      BOOST_CHECK_CLOSE(dataSample1.get(rowIdx, colIdx), testPoints[rowIdx][colIdx], tolerance);
    }
    // this only works because dataset does not contain any zeros. if the dataset to test with is
    // changed, we need a else if here and a boost check small.
    BOOST_CHECK_CLOSE(classesSample1.get(rowIdx), testValues[rowIdx], tolerance);
  }

  dataset.reset();

  // check if we get the correct samples in a second run
  dataset = std::unique_ptr<Dataset>(sampleProvider.getNextSamples(sampleSize2));
  DataVector& classesSample2 = dataset->getTargets();
  DataMatrix& dataSample2 = dataset->getData();
  nrows = dataSample2.getNrows();
  ncols = dataSample2.getNcols();

  // Check if all dimensions agree
  BOOST_CHECK_EQUAL(sampleSize2, classesSample2.getSize());
  BOOST_CHECK_EQUAL(sampleSize2, dataSample2.getNrows());
  BOOST_CHECK_EQUAL(datasetDim, dataSample2.getNcols());

  testVector = DataVector(ncols);
  for (size_t rowIdx = 0; rowIdx < nrows; rowIdx++) {
    dataSample2.getRow(rowIdx, testVector);
    for (size_t colIdx = 0; colIdx < ncols; colIdx++) {
      // this only works because dataset does not contain any zeros. if the dataset to test with
      // changed, we need a else if here and a boost check small.
      BOOST_CHECK_CLOSE(dataSample2.get(rowIdx, colIdx), testPoints[rowIdx + sampleSize1][colIdx],
                        tolerance);
    }
    // this only works because dataset does not contain any zeros. if the dataset to test with is
    // changed, we need a else if here and a boost check small.
    BOOST_CHECK_CLOSE(classesSample2.get(rowIdx), testValues[rowIdx + sampleSize1], tolerance);
  }
}

BOOST_AUTO_TEST_SUITE_END()
