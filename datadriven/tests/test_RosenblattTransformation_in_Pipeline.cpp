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
  double tolerance = 1e-14;

  // read arff file
  ArffFileSampleProvider arffsp = ArffFileSampleProvider();
  arffsp.readFile("datadriven/tests/datasets/liver-disorders_normalized.arff");
  Dataset* dataset = arffsp.getAllSamples();
  DataSourceConfig config;
  config.dataTransformation = DataTransformationType::ROSENBLATT;
  config.numSamplesForTranformation = 2;

  // do transformations
  DataTransformationBuilder dataTrBuilder;
  DataTransformation* dataTransformation = dataTrBuilder.buildTransformation(config, dataset);
  std::cout << "Wrapper: dataTransformation gebuildet" << std::endl;

  Dataset* datasetTr = dataTransformation->doTransformation(dataset);
    std::cout << "Wrapper: erste Transformation" << std::endl;
  Dataset* datasetInvTr = dataTransformation->doInverseTransformation(datasetTr);
  std::cout << "Wrapper: dataTransformation durchgeführt" << std::endl;

  // check error between original and transformed dataset
  DataVector sampleOrigin(dataset->getDimension());
  DataVector sampleTransformed(dataset->getDimension());

  /*
  for (size_t i=0; i < datasetTr->getNumberInstances(); i++) {
    std::cout << "Sample " << i << ":" << std::endl;
    dataset->getData().getRow(i, sampleOrigin);
    std::cout << "Dataset: " << sampleOrigin[0] << " " <<
        sampleOrigin[1] << " " << sampleOrigin[2] << std::endl;
    datasetTr->getData().getRow(i, sampleOrigin);
    std::cout << "DatasetTr: " << sampleOrigin[0] << " " <<
        sampleOrigin[1] << " " << sampleOrigin[2] << std::endl;
    datasetInvTr->getData().getRow(i, sampleOrigin);
    std::cout << "DatasetInvTr: " << sampleOrigin[0] << " " <<
        sampleOrigin[1] << " " << sampleOrigin[2] << std::endl;
    }
    */

  std::cout << "Wrapper: dataset. Samples: " << dataset->getNumberInstances()
        << ", Dimension: " << dataset->getDimension() << std::endl;
  std::cout << "Wrapper: datasetInvTr. Samples: " << datasetInvTr->getNumberInstances()
        << ", Dimension: " << datasetInvTr->getDimension() << std::endl;

  for (size_t isample = 0; isample < dataset->getNumberInstances(); isample++) {
    dataset->getData().getRow(isample, sampleOrigin);
    datasetInvTr->getData().getRow(isample, sampleTransformed);
    for (size_t idim = 0; idim < dataset->getDimension(); idim++) {
      std::cout << "SampleOrigin: " << sampleOrigin[idim]
                      << ", SampleTransformed: " << sampleTransformed[idim] << std::endl;
      // assert that sampleOrigin and sampleTransformed contain the same samples
      double inversionError =
          std::abs(sampleOrigin[idim] - sampleTransformed[idim]) / sampleOrigin[idim];

      BOOST_CHECK_SMALL(inversionError, tolerance);
    }
  }
}


BOOST_AUTO_TEST_CASE(testDataTransformationParser) {
  double tolerance = 1e-14;

  // Config for "automatic" transformation
  std::string path = "datadriven/tests/data/test_Rosenblatt_in_Pipeline.json";
  DataSourceBuilder builder;
  DataSourceConfig config;
  DataSourceConfig defaults;

  std::cout << "Parser: Los gehts automatisch mit Transformationstyp: "
      << DataTransformationTypeParser::toString(config.dataTransformation) << std::endl;

  // "automatic" transformation in DataSource->getNextSamples
  DataMiningConfigParser parser(path);
  parser.getDataSourceConfig(config, defaults);

  std::cout << "Parser: getDataSourceConfig, type: "
      << DataTransformationTypeParser::toString(config.dataTransformation) << std::endl;

  auto dataSource = builder.fromConfig(config);
  Dataset* datasetAuto = (dataSource->getNextSamples());

  std::cout << "Parser: datasetAuto. Samples: " << datasetAuto->getNumberInstances()
      << ", Dimension: " << datasetAuto->getDimension() << std::endl;
  std::cout << "Parser: weiter manuell" << std::endl;

  // Read arff file manually
  ArffFileSampleProvider arffsp = ArffFileSampleProvider();
  arffsp.readFile("datadriven/tests/datasets/liver-disorders_normalized.arff");
  Dataset* dataset = arffsp.getAllSamples();

  std::cout << "Parser: dataset. Samples: " << dataset->getNumberInstances()
      << ", Dimension: " << dataset->getDimension() << std::endl;

  // "manual" transformation
  DataTransformationBuilder dataTrBuilder;
  DataTransformation* dataTransformation = dataTrBuilder.buildTransformation(config, dataset);
  Dataset* datasetMan = dataTransformation->doTransformation(dataset);

  std::cout << "Parser: datasetMan. Samples: " << datasetMan->getNumberInstances()
        << ", Dimension: " << datasetMan->getDimension() << std::endl;

  std::cout << "Parser: manuell transformed" << std::endl;

  // check error between original and transformed dataset
  DataVector sampleAuto(dataset->getDimension());
  DataVector sampleMan(dataset->getDimension());
  std::cout << "Parser: Jetzt Fehler checken" << std::endl;

  for (size_t isample = 0; isample < dataset->getNumberInstances(); isample++) {
    datasetAuto->getData().getRow(isample, sampleAuto);
    datasetMan->getData().getRow(isample, sampleMan);

    for (size_t idim = 0; idim < dataset->getDimension(); idim++) {
      // assert that sampleAuto and sampleMan contain the same samples
      std::cout << "SampleAuto: " << sampleAuto[idim]
                << ", SampleMan: " << sampleMan[idim] << std::endl;
      double inversionError =
          std::abs(sampleAuto[idim] - sampleMan[idim]) / sampleAuto[idim];
      BOOST_CHECK_SMALL(inversionError, tolerance);
    }
  }
}


BOOST_AUTO_TEST_SUITE_END()

