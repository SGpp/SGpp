/* Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 * DBMatOfflineDatabaseTest.cpp
 *
 * densityEstimationPipelineTest.cpp
 *
 *  Created on: May 27, 2018
 *      Author: dominik
 */
#ifdef USE_GSL

#include <boost/test/test_tools.hpp>
#include <boost/test/unit_test_suite.hpp>
#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/datadriven/datamining/base/SparseGridMiner.hpp>
#include <sgpp/datadriven/datamining/builder/ClassificationMinerFactory.hpp>
#include <sgpp/datadriven/datamining/modules/dataSource/CSVFileSampleProvider.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/FitterConfigurationDensityEstimation.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingBase.hpp>
#include <sgpp/datadriven/scalapack/BlacsProcessGrid.hpp>

#include <cmath>
#include <fstream>
#include <iostream>
#include <random>
#include <string>

using sgpp::datadriven::BlacsProcessGrid;
using sgpp::datadriven::ClassificationMinerFactory;
using sgpp::datadriven::CSVFileSampleProvider;
using sgpp::datadriven::DataVector;
using sgpp::datadriven::ModelFittingBase;
using sgpp::datadriven::SparseGridMiner;

double testModel(std::string configFile) {
  ClassificationMinerFactory factory;
  SparseGridMiner *miner = factory.buildMiner(configFile);
  miner->learn(false);
  ModelFittingBase *model = miner->getModel();

  // Test
  CSVFileSampleProvider csv;
  csv.readFile("datadriven/datasets/gmm/gmm_test.csv", true);
  auto testDataset = *(csv.getAllSamples());
  DataVector predictions(testDataset.getNumberInstances());
  model->evaluate(testDataset.getData(), predictions);
  size_t correct = 0;
  for (size_t idx = 0; idx < predictions.size(); idx++) {
    if (predictions.get(idx) == testDataset.getTargets().get(idx)) {
      correct++;
    }
  }

  if (predictions.size() == 0) {
    return 0.0;
  }

  double accuracy = static_cast<double>(correct) / static_cast<double>(predictions.size());
  return accuracy;
}

BOOST_AUTO_TEST_SUITE(testClassification)

BOOST_AUTO_TEST_CASE(testOnOff) {
  std::string configFile = "datadriven/tests/gmm_on_off.json";
  double accuracy = testModel(configFile);
  std::cout << "Accuracy " << accuracy << std::endl;
  BOOST_CHECK(accuracy > 0.7);
}
BOOST_AUTO_TEST_CASE(testCG) {
  std::string configFile = "datadriven/tests/gmm_cg.json";
  double accuracy = testModel(configFile);
  std::cout << "Accuracy " << accuracy << std::endl;
  BOOST_CHECK(accuracy > 0.7);
}

#ifdef USE_SCALAPACK
BOOST_AUTO_TEST_CASE(testOnOffParallelChol) {
  std::string configFileParallel = "datadriven/tests/gmm_on_off_parallel_chol.json";
  double accuracyParallel = testModel(configFileParallel);

  std::string configFile = "datadriven/tests/gmm_on_off_chol.json";
  double accuracy = testModel(configFile);
  if (BlacsProcessGrid::getCurrentProcess() == 0) {
    BOOST_CHECK(accuracyParallel > 0.7);
    BOOST_CHECK_CLOSE(accuracyParallel, accuracy, 1e-5);
  }
}

BOOST_AUTO_TEST_CASE(testOnOffParallelOrthoadapt) {
  std::string configFileParallel = "datadriven/tests/gmm_on_off_parallel_orthoadapt.json";
  double accuracyParallel = testModel(configFileParallel);

  std::string configFile = "datadriven/tests/gmm_on_off_orthoadapt.json";
  double accuracy = testModel(configFile);

  if (BlacsProcessGrid::getCurrentProcess() == 0) {
    BOOST_CHECK(accuracyParallel > 0.7);
    BOOST_CHECK_CLOSE(accuracyParallel, accuracy, 1e-5);
  }
}
#endif /* USE_SCALAPACK */

BOOST_AUTO_TEST_SUITE_END()

#endif /* USE_GSL */
