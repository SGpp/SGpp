// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <boost/test/test_tools.hpp>
#include <boost/test/unit_test_suite.hpp>
#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/datadriven/datamining/base/SparseGridMiner.hpp>
#include <sgpp/datadriven/datamining/builder/ClassificationMinerFactory.hpp>
#include <sgpp/datadriven/datamining/builder/DensityEstimationMinerFactory.hpp>
#include <sgpp/datadriven/datamining/modules/dataSource/CSVFileSampleProvider.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/FitterConfigurationDensityEstimation.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingBase.hpp>
#include <sgpp/datadriven/scalapack/BlacsProcessGrid.hpp>

#include <iostream>
#include <string>

using sgpp::datadriven::BlacsProcessGrid;
using sgpp::datadriven::ClassificationMinerFactory;
using sgpp::datadriven::DensityEstimationMinerFactory;
using sgpp::datadriven::CSVFileSampleProvider;
using sgpp::datadriven::DataVector;
using sgpp::datadriven::ModelFittingBase;
using sgpp::datadriven::SparseGridMiner;

double testModelClassification(std::string configFile) {
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

double testModelDE(std::string configFile) {
  DensityEstimationMinerFactory factory;
  auto miner = std::unique_ptr<SparseGridMiner>(factory.buildMiner(configFile));
  miner->learn(false);
  ModelFittingBase *model = miner->getModel();

  // Test
  CSVFileSampleProvider csv;
  csv.readFile("datadriven/datasets/densityEstimation/2D_StroSkewB2F.csv", true);
  auto testDataset = *(csv.getAllSamples());
  DataVector predictions(testDataset.getNumberInstances());
  model->evaluate(testDataset.getData(), predictions);

  // Calculate MSE
  predictions.sub(testDataset.getTargets());
  return predictions.l2Norm() / static_cast<double>(predictions.getSize());
}

BOOST_AUTO_TEST_SUITE(testLambdaOptimization)
#ifdef USE_SCALAPACK

BOOST_AUTO_TEST_CASE(testLambdaOptimizationClassificationParallelOrthoadapt) {
  std::string configFileParallel =
      "datadriven/tests/gmm_lambda_optimization_parallel_orthoadapt.json";
  double accuracyLambdaOptimization = testModelClassification(configFileParallel);

  std::string configFile = "datadriven/tests/gmm_on_off_orthoadapt.json";
  double accuracy = testModelClassification(configFile);
  if (BlacsProcessGrid::getCurrentProcess() == 0) {
    BOOST_CHECK(accuracyLambdaOptimization > 0.7);
    BOOST_CHECK(accuracyLambdaOptimization > accuracy);
  }
}

BOOST_AUTO_TEST_CASE(testLambdaOptimizationClassificationParallelOrthoadaptSMW) {
  std::string configFileParallel =
      "datadriven/tests/gmm_lambda_optimization_parallel_orthoadapt_smw.json";
  double accuracyLambdaOptimization = testModelClassification(configFileParallel);

  std::string configFile = "datadriven/tests/gmm_on_off_orthoadapt.json";
  double accuracy = testModelClassification(configFile);

  if (BlacsProcessGrid::getCurrentProcess() == 0) {
    BOOST_CHECK(accuracyLambdaOptimization > 0.7);
    BOOST_CHECK(accuracyLambdaOptimization > accuracy);
  }
}

BOOST_AUTO_TEST_CASE(testLambdaOptimizationParallelOrthoadapt) {
  std::string configFileLambdaOptimization =
      "datadriven/tests/de_lambda_optimization_parallel_orthoadapt.json";
  double mseLambdaOptimization = testModelDE(configFileLambdaOptimization);

  std::string configFile = "datadriven/tests/de_lambda_optimization_reference.json";
  double mse = testModelClassification(configFile);

  if (BlacsProcessGrid::getCurrentProcess() == 0) {
    BOOST_CHECK(mseLambdaOptimization < 5e-2);
    BOOST_CHECK(mseLambdaOptimization < mse);
  }
}

BOOST_AUTO_TEST_CASE(testLambdaOptimizationParallelOrthoadaptSMW) {
  std::string configFileLambdaOptimization =
      "datadriven/tests/de_lambda_optimization_parallel_orthoadapt_smw.json";
  double mseLambdaOptimization = testModelDE(configFileLambdaOptimization);

  std::string configFile = "datadriven/tests/de_lambda_optimization_reference.json";
  double mse = testModelDE(configFile);

  if (BlacsProcessGrid::getCurrentProcess() == 0) {
    BOOST_CHECK(mseLambdaOptimization < 5e-2);
    BOOST_CHECK(mseLambdaOptimization < mse);
  }
}
#endif /* USE_SCALAPACK */

#ifdef USE_GSL

BOOST_AUTO_TEST_CASE(testLambdaOptimizationClassificationOrthoadapt) {
  std::string configFileLambdaOptimization =
      "datadriven/tests/gmm_lambda_optimization_orthoadapt.json";
  double accuracyLambdaOptimization = testModelClassification(configFileLambdaOptimization);

  std::string configFile = "datadriven/tests/gmm_on_off_orthoadapt.json";
  double accuracy = testModelClassification(configFile);

  BOOST_CHECK(accuracyLambdaOptimization > 0.7);
  BOOST_CHECK(accuracyLambdaOptimization > accuracy);
}

BOOST_AUTO_TEST_CASE(testLambdaOptimizationClassificationOrthoadaptSMW) {
  std::string configFileLambdaOptimization =
      "datadriven/tests/gmm_lambda_optimization_orthoadapt_smw.json";
  double accuracyLambdaOptimization = testModelClassification(configFileLambdaOptimization);

  std::string configFile = "datadriven/tests/gmm_on_off_orthoadapt.json";
  double accuracy = testModelClassification(configFile);

  BOOST_CHECK(accuracyLambdaOptimization > 0.7);
  BOOST_CHECK(accuracyLambdaOptimization > accuracy);
}

BOOST_AUTO_TEST_CASE(testLambdaOptimizationOrthoadapt) {
  std::string configFileLambdaOptimization =
      "datadriven/tests/de_lambda_optimization_orthoadapt.json";
  double mseLambdaOptimization = testModelDE(configFileLambdaOptimization);

  std::string configFile = "datadriven/tests/de_lambda_optimization_reference.json";
  double mse = testModelClassification(configFile);

  BOOST_CHECK(mseLambdaOptimization < 5e-2);
  BOOST_CHECK(mseLambdaOptimization < mse);
}

BOOST_AUTO_TEST_CASE(testLambdaOptimizationOrthoadaptSMW) {
  std::string configFileLambdaOptimization =
      "datadriven/tests/de_lambda_optimization_orthoadapt_smw.json";
  double mseLambdaOptimization = testModelDE(configFileLambdaOptimization);

  std::string configFile = "datadriven/tests/de_lambda_optimization_reference.json";
  double mse = testModelDE(configFile);

  BOOST_CHECK(mseLambdaOptimization < 5e-2);
  BOOST_CHECK(mseLambdaOptimization < mse);
}

#endif /* USE_GSL */

BOOST_AUTO_TEST_SUITE_END()
