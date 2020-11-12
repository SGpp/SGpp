// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifdef ZLIB

#define BOOST_TEST_DYN_LINK

#include <boost/test/unit_test.hpp>

#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/datadriven/datamining/modules/dataSource/ArffFileSampleProvider.hpp>
#include <sgpp/datadriven/datamining/modules/dataSource/DataSource.hpp>
#include <sgpp/datadriven/datamining/modules/dataSource/DataSourceSplitting.hpp>
#include <sgpp/datadriven/datamining/modules/dataSource/DataSourceConfig.hpp>
#include <sgpp/datadriven/datamining/modules/dataSource/GzipFileSampleDecorator.hpp>
#include <sgpp/datadriven/tools/Dataset.hpp>
#include <sgpp/globaldef.hpp>

#include <array>
#include <memory>
#include <string>

using sgpp::datadriven::DataSource;
using sgpp::datadriven::DataSourceSplitting;
using sgpp::datadriven::Dataset;
using sgpp::datadriven::SampleProvider;
using sgpp::datadriven::GzipFileSampleDecorator;
using sgpp::datadriven::ArffFileSampleProvider;
using sgpp::datadriven::DataSourceConfig;
using sgpp::base::DataMatrix;
using sgpp::base::DataVector;

struct State {
  State()
      : path("datadriven/datasets/liver/liver-disorders_normalized_small.arff.gz"),
        testPoints({{{{0.307143, 0.130137, 0.050000}},
                     {{0.365584, 0.105479, 0.050000}},
                     {{0.178571, 0.201027, 0.050000}},
                     {{0.272078, 0.145548, 0.050000}},
                     {{0.318831, 0.065411, 0.050000}},
                     {{0.190260, 0.086986, 0.050000}},
                     {{0.190260, 0.062329, 0.072500}},
                     {{0.120130, 0.068493, 0.072500}},
                     {{0.225325, 0.056164, 0.072500}},
                     {{0.213636, 0.050000, 0.072500}}}}),
        testValues({{-1., 1., 1., 1., 1., 1., -1., -1., -1., -1.}}) {}
  ~State() {}
  std::string path;
  std::array<std::array<double, 3>, 10> testPoints;
  std::array<double, 10> testValues;
};

// TODO(Michael Lettrich): Test with other configurations like different batch sizes
BOOST_FIXTURE_TEST_SUITE(dataSourceGetNextSamplesAllTest, State)

BOOST_AUTO_TEST_CASE(dataSourcegetNextSamplesAllSamplesTest) {
  SampleProvider* sampleProvider = new GzipFileSampleDecorator(new ArffFileSampleProvider());
  DataSourceConfig config;
  config.filePath_ = path;

  DataSource* dataSource = new DataSourceSplitting(config, sampleProvider);
  auto dataset = std::unique_ptr<Dataset>(dataSource->getNextSamples());

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
  delete dataSource;
}

BOOST_AUTO_TEST_CASE(dataSourceGetAllIteratorTest) {
  SampleProvider* sampleProvider = new GzipFileSampleDecorator(new ArffFileSampleProvider());
  DataSourceConfig config;
  config.filePath_ = path;

  DataSource* dataSource = new DataSourceSplitting(config, sampleProvider);

  for (auto dsPtr : *dataSource) {
    auto dataset = std::unique_ptr<Dataset>(dsPtr);
    DataVector& classes = dataset->getTargets();
    DataMatrix& data = dataset->getData();

    // Check if all dimensions agree
    BOOST_CHECK_EQUAL(10, classes.getSize());
    BOOST_CHECK_EQUAL(10, data.getNrows());
    BOOST_CHECK_EQUAL(3, data.getNcols());
    //    dataset.reset();
  }
  BOOST_CHECK_EQUAL(1, dataSource->getCurrentIteration());
  delete dataSource;
}

BOOST_AUTO_TEST_SUITE_END()
#endif
