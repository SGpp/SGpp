// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#define BOOST_TEST_DYN_LINK

#include <boost/test/unit_test.hpp>

#include <sgpp/datadriven/datamining/builder/DataSourceBuilder.hpp>
#include <sgpp/datadriven/datamining/configuration/DataMiningConfigParser.hpp>
#include <sgpp/datadriven/datamining/modules/dataSource/ArffFileSampleProvider.hpp>
#include <sgpp/datadriven/datamining/modules/dataSource/DataSource.hpp>
#include <sgpp/datadriven/datamining/modules/dataSource/DataTransformation.hpp>
#include <sgpp/datadriven/datamining/modules/dataSource/DataTransformationBuilder.hpp>
#include <sgpp/datadriven/datamining/modules/dataSource/DataTransformationTypeParser.hpp>
#include <sgpp/datadriven/datamining/modules/dataSource/RosenblattTransformation.hpp>

#include <vector>
#include <string>

using sgpp::datadriven::ArffFileSampleProvider;
using sgpp::datadriven::DataSourceBuilder;
using sgpp::datadriven::DataSourceConfig;
using sgpp::datadriven::Dataset;
using sgpp::datadriven::DataTransformation;
using sgpp::datadriven::DataTransformationBuilder;
using sgpp::datadriven::DataTransformationType;
using sgpp::datadriven::DataTransformationTypeParser;
using sgpp::datadriven::RosenblattTransformation;
using sgpp::datadriven::DataMiningConfigParser;
using sgpp::base::DataVector;

BOOST_AUTO_TEST_SUITE(testRosenblattTransformationInPipeline)

BOOST_AUTO_TEST_CASE(testRosenblattWrapper) {
  double tolerance = 1e-10;
  DataSourceConfig config;
  config.dataTransformation = DataTransformationType::ROSENBLATT;
  config.numSamplesForTranformation = 1000;

  // read arff file
  ArffFileSampleProvider arffsp = ArffFileSampleProvider();
  arffsp.readFile("datadriven/tests/data/chess_5d_2000.arff");
  Dataset* dataset = arffsp.getAllSamples();

  // do transformations
  DataTransformationBuilder dataTrBuilder;
  DataTransformation* dataTransformation = dataTrBuilder.buildTransformation(config, dataset);

  Dataset* datasetTr = dataTransformation->doTransformation(dataset);
  Dataset* datasetInvTr = dataTransformation->doInverseTransformation(datasetTr);

  // check error between original and transformed dataset
  DataVector sampleOrigin(dataset->getDimension());
  DataVector sampleTransformed(dataset->getDimension());

  for (size_t isample = 0; isample < dataset->getNumberInstances(); isample++) {
    dataset->getData().getRow(isample, sampleOrigin);
    datasetInvTr->getData().getRow(isample, sampleTransformed);
    for (size_t idim = 0; idim < dataset->getDimension(); idim++) {
      // assert that sampleOrigin and sampleTransformed contain the same samples
      double inversionError =
          std::abs(sampleOrigin[idim] - sampleTransformed[idim]) / sampleOrigin[idim];
      BOOST_CHECK_SMALL(inversionError, tolerance);
    }
  }
}

BOOST_AUTO_TEST_CASE(testDataTransformationParser) {
  double tolerance = 1e-10;
  DataSourceConfig config;
  DataSourceConfig defaults;

  // Config for "automatic" transformation
  std::string path = "datadriven/tests/data/test_Rosenblatt_in_Pipeline.json";
  DataSourceBuilder builder;
  DataMiningConfigParser parser(path);
  parser.getDataSourceConfig(config, defaults);

  // "automatic" transformation in DataSource->getNextSamples()
  auto dataSource = builder.fromConfig(config);
  Dataset* datasetAuto = dataSource->getNextSamples();

  // Read arff file manually
  ArffFileSampleProvider arffsp = ArffFileSampleProvider();
  arffsp.readFile("datadriven/tests/data/chess_5d_2000.arff");
  Dataset* dataset = arffsp.getAllSamples();

  // "manual" transformation
  DataTransformationBuilder dataTrBuilder;
  DataTransformation* dataTransformation = dataTrBuilder.buildTransformation(config, dataset);
  Dataset* datasetMan = dataTransformation->doTransformation(dataset);

  // check error between original and transformed dataset
  DataVector sampleAuto(dataset->getDimension());
  DataVector sampleMan(dataset->getDimension());

  for (size_t isample = 0; isample < dataset->getNumberInstances(); isample++) {
    datasetAuto->getData().getRow(isample, sampleAuto);
    datasetMan->getData().getRow(isample, sampleMan);

    for (size_t idim = 0; idim < dataset->getDimension(); idim++) {
      // assert that sampleAuto and sampleMan contain the same samples
      double inversionError =
          std::abs(sampleAuto[idim] - sampleMan[idim]) / sampleAuto[idim];
      BOOST_CHECK_SMALL(inversionError, tolerance);
    }
  }
}

BOOST_AUTO_TEST_SUITE_END()
