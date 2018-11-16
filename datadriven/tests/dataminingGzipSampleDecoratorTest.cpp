/* This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * datamingArffSampleProviderTest.cpp
 *
 *  Created on: 01.04.2016
 *      Author: Michael Lettrich
 */
#ifdef ZLIB

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/datadriven/datamining/modules/dataSource/ArffFileSampleProvider.hpp>
#include <sgpp/datadriven/datamining/modules/dataSource/GzipFileSampleDecorator.hpp>
#include <sgpp/datadriven/tools/Dataset.hpp>
#include <sgpp/globaldef.hpp>

#include <string>

BOOST_AUTO_TEST_SUITE(dataminingGzipSampleDecoratorTest)

using sgpp::datadriven::GzipFileSampleDecorator;
using sgpp::datadriven::ArffFileSampleProvider;
using sgpp::base::DataMatrix;
using sgpp::base::DataVector;
using sgpp::datadriven::Dataset;

BOOST_AUTO_TEST_CASE(gzipTestReadFile) {
  double testPoints[10][3] = {{0.307143, 0.130137, 0.050000}, {0.365584, 0.105479, 0.050000},
                              {0.178571, 0.201027, 0.050000}, {0.272078, 0.145548, 0.050000},
                              {0.318831, 0.065411, 0.050000}, {0.190260, 0.086986, 0.050000},
                              {0.190260, 0.062329, 0.072500}, {0.120130, 0.068493, 0.072500},
                              {0.225325, 0.056164, 0.072500}, {0.213636, 0.050000, 0.072500}};

  double testValues[10] = {-1., 1., 1., 1., 1., 1., -1., -1., -1., -1.};

  GzipFileSampleDecorator sampleProvider = GzipFileSampleDecorator(new ArffFileSampleProvider());

  const auto datasetPath = "datadriven/tests/datasets/liver-disorders_normalized.arff.gz";

  sampleProvider.readFile(datasetPath, true);
  auto dataset = sampleProvider.getAllSamples();

  DataVector& classes = dataset->getTargets();
  DataMatrix& data = dataset->getData();
  size_t nrows = data.getNrows();
  size_t ncols = data.getNcols();

  double tolerance = 1E-5;

  // Check if all dimensions agree
  BOOST_CHECK_EQUAL(10, classes.getSize());
  BOOST_CHECK_EQUAL(10, data.getNrows());
  BOOST_CHECK_EQUAL(3, data.getNcols());

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

BOOST_AUTO_TEST_SUITE_END()
#endif
