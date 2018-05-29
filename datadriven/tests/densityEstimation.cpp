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

#include <boost/test/unit_test_suite.hpp>
#include <boost/test/test_tools.hpp>
#include <sgpp/datadriven/datamining/modules/dataSource/CSVFileSampleProvider.hpp>
#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/datadriven/datamining/base/SparseGridMiner.hpp>
#include <sgpp/datadriven/datamining/builder/DensityEstimationMinerFactory.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingBase.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/FitterConfigurationDensityEstimation.hpp>

#include <string>
#include <iostream>
#include <fstream>
#include <random>
#include <cmath>

using sgpp::datadriven::DensityEstimationMinerFactory;
using sgpp::datadriven::SparseGridMiner;
using sgpp::datadriven::ModelFittingBase;
using sgpp::datadriven::CSVFileSampleProvider;
using sgpp::datadriven::DataVector;

double testDistribution(std::string testCSV, std::string config) {
  // Train model
  DensityEstimationMinerFactory factory;
  auto miner = std::unique_ptr<SparseGridMiner>(factory.buildMiner(config));
  miner->learn();
  ModelFittingBase *model = miner->getModel();
  // Test
  CSVFileSampleProvider csv;
  csv.readFile(testCSV);
  auto testDataset = *(csv.getAllSamples());
  DataVector predictions(testDataset.getNumberInstances());
  model->evaluate(testDataset.getData(), predictions);
  // Calculate MSE
  predictions.sub(testDataset.getTargets());
  return predictions.l2Norm() / static_cast<double>(predictions.getSize());
}

BOOST_AUTO_TEST_SUITE(testDensityEstimation)


BOOST_AUTO_TEST_CASE(Test_2D_StroSkewB2) {
  std::string samples = "datadriven/tests/datasets/densityEstimation/2D_StroSkewB2.csv";

  // Create config file
  std::string config = "tmpsgdeconfig.json";
  std::ofstream stream(config);
  stream << "{" << "\"dataSource\" : { \"filePath\" : \"" << samples << "\"},\"scorer\" : "
      << "{ \"testing\" : { \"metric\" : \"MSE\"}},\"fitter\" : " <<
      "{ \"type\" : \"densityEstimation\", \"gridConfig\" : { \"gridType\" : \"linear\","
      << "\"level\" : 5},\"adaptivityConfig\" : {\"numRefinements\" : 3, \"threshold\" : 0.001,"
      << "\"maxLevelType\" : false, \"noPoints\" : 3},\"regularizationConfig\" : {\"lambda\" : "
      << "1}}}" << std::endl;

  double mse = testDistribution(
      "datadriven/tests/datasets/densityEstimation/2D_StroSkewB2F.csv", config);
  std::cout << "MSE on test " << mse << std::endl;
  BOOST_CHECK(mse <= 5e-2);
  remove(config.c_str());
}


BOOST_AUTO_TEST_CASE(Test_3D_KurB4B1) {
  std::string samples = "datadriven/tests/datasets/densityEstimation/3D_KurB4B1.csv";

  // Create config file
  std::string config = "tmpsgdeconfig.json";
  std::ofstream stream(config);
  stream << "{" << "\"dataSource\" : { \"filePath\" : \"" << samples << "\"},\"scorer\" : "
      << "{ \"testing\" : { \"metric\" : \"MSE\"}},\"fitter\" : " <<
      "{ \"type\" : \"densityEstimation\", \"gridConfig\" : { \"gridType\" : \"linear\","
      << "\"level\" : 5},\"adaptivityConfig\" : {\"numRefinements\" : 3, \"threshold\" : 0.001,"
      << "\"maxLevelType\" : false, \"noPoints\" : 3},\"regularizationConfig\" : {\"lambda\" : "
      << "1}}}" << std::endl;

  double mse = testDistribution(
      "datadriven/tests/datasets/densityEstimation/3D_KurB4B1F.csv", config);
  std::cout << "MSE on test " << mse << std::endl;
  BOOST_CHECK(mse <= 5e-2);
  remove(config.c_str());
}


BOOST_AUTO_TEST_CASE(Test_4D_SepBiKurB2B4) {
  std::string samples = "datadriven/tests/datasets/densityEstimation/4D_SepBiKurB2B4.csv";

  // Create config file
  std::string config = "tmpsgdeconfig.json";
  std::ofstream stream(config);
  stream << "{" << "\"dataSource\" : { \"filePath\" : \"" << samples << "\"},\"scorer\" : "
      << "{ \"testing\" : { \"metric\" : \"MSE\"}},\"fitter\" : " <<
      "{ \"type\" : \"densityEstimation\", \"gridConfig\" : { \"gridType\" : \"linear\","
      << "\"level\" : 5},\"adaptivityConfig\" : {\"numRefinements\" : 3, \"threshold\" : 0.001,"
      << "\"maxLevelType\" : false, \"noPoints\" : 3},\"regularizationConfig\" : {\"lambda\" : "
      << "1}}}" << std::endl;

  double mse = testDistribution(
      "datadriven/tests/datasets/densityEstimation/4D_SepBiKurB2B4F.csv", config);
  std::cout << "MSE on test " << mse << std::endl;
  BOOST_CHECK(mse <= 5e-2);
  remove(config.c_str());
}


BOOST_AUTO_TEST_CASE(Test_5D_B1B4SepBiB3G) {
  std::string samples = "datadriven/tests/datasets/densityEstimation/5D_B1B4SepBiB3G.csv";

  // Create config file
  std::string config = "tmpsgdeconfig.json";
  std::ofstream stream(config);
  stream << "{" << "\"dataSource\" : { \"filePath\" : \"" << samples << "\"},\"scorer\" : "
      << "{ \"testing\" : { \"metric\" : \"MSE\"}},\"fitter\" : " <<
      "{ \"type\" : \"densityEstimation\", \"gridConfig\" : { \"gridType\" : \"linear\","
      << "\"level\" : 5},\"adaptivityConfig\" : {\"numRefinements\" : 3, \"threshold\" : 0.001,"
      << "\"maxLevelType\" : false, \"noPoints\" : 3},\"regularizationConfig\" : {\"lambda\" : "
      << "1}}}" << std::endl;

  double mse = testDistribution(
      "datadriven/tests/datasets/densityEstimation/5D_B1B4SepBiB3GF.csv", config);
  std::cout << "MSE on test " << mse << std::endl;
  BOOST_CHECK(mse <= 5e-2);
  remove(config.c_str());
}


BOOST_AUTO_TEST_CASE(Test_6D_KurStroSkewGClawB3SepBi) {
  std::string samples =
      "datadriven/tests/datasets/densityEstimation/6D_KurStroSkewGClawB3SepBi.csv";

  // Create config file
  std::string config = "tmpsgdeconfig.json";
  std::ofstream stream(config);
  stream << "{" << "\"dataSource\" : { \"filePath\" : \"" << samples << "\"},\"scorer\" : "
      << "{ \"testing\" : { \"metric\" : \"MSE\"}},\"fitter\" : " <<
      "{ \"type\" : \"densityEstimation\", \"gridConfig\" : { \"gridType\" : \"linear\","
      << "\"level\" : 5},\"adaptivityConfig\" : {\"numRefinements\" : 3, \"threshold\" : 0.001,"
      << "\"maxLevelType\" : false, \"noPoints\" : 3},\"regularizationConfig\" : {\"lambda\" : "
      << "1}}}" << std::endl;

  double mse = testDistribution(
      "datadriven/tests/datasets/densityEstimation/6D_KurStroSkewGClawB3SepBiF.csv", config);
  std::cout << "MSE on test " << mse << std::endl;
  BOOST_CHECK(mse <= 2e-1);
  remove(config.c_str());
}


BOOST_AUTO_TEST_CASE(Test_7D_StroSkewB2KurB4SepBiGB3) {
  std::string samples =
      "datadriven/tests/datasets/densityEstimation/7D_StroSkewB2KurB4SepBiGB3.csv";

  // Create config file
  std::string config = "tmpsgdeconfig.json";
  std::ofstream stream(config);
  stream << "{" << "\"dataSource\" : { \"filePath\" : \"" << samples << "\"},\"scorer\" : "
      << "{ \"testing\" : { \"metric\" : \"MSE\"}},\"fitter\" : " <<
      "{ \"type\" : \"densityEstimation\", \"gridConfig\" : { \"gridType\" : \"linear\","
      << "\"level\" : 4},\"adaptivityConfig\" : {\"numRefinements\" : 3, \"threshold\" : 0.001,"
      << "\"maxLevelType\" : false, \"noPoints\" : 3},\"regularizationConfig\" : {\"lambda\" : "
      << "1}}}" << std::endl;

  double mse = testDistribution(
      "datadriven/tests/datasets/densityEstimation/7D_StroSkewB2KurB4SepBiGB3F.csv", config);
  std::cout << "MSE on test " << mse << std::endl;
  BOOST_CHECK(mse <= 2e-1);
  remove(config.c_str());
}


BOOST_AUTO_TEST_CASE(Test_8D_B2KurGB3B4ClawStroSkewB1) {
  std::string samples =
      "datadriven/tests/datasets/densityEstimation/8D_B2KurGB3B4ClawStroSkewB1.csv";

  // Create config file
  std::string config = "tmpsgdeconfig.json";
  std::ofstream stream(config);
  stream << "{" << "\"dataSource\" : { \"filePath\" : \"" << samples << "\"},\"scorer\" : "
      << "{ \"testing\" : { \"metric\" : \"MSE\"}},\"fitter\" : " <<
      "{ \"type\" : \"densityEstimation\", \"gridConfig\" : { \"gridType\" : \"linear\","
      << "\"level\" : 4},\"adaptivityConfig\" : {\"numRefinements\" : 3, \"threshold\" : 0.001,"
      << "\"maxLevelType\" : false, \"noPoints\" : 3},\"regularizationConfig\" : {\"lambda\" : "
      << "1}}}" << std::endl;

  double mse = testDistribution(
      "datadriven/tests/datasets/densityEstimation/8D_B2KurGB3B4ClawStroSkewB1F.csv", config);
  std::cout << "MSE on test " << mse << std::endl;
  BOOST_CHECK(mse <= 1e-1);
  remove(config.c_str());
}


BOOST_AUTO_TEST_CASE(Test_9D_KurB1B2StroSkewSepBiClawB3B4Tri) {
  std::string samples =
      "datadriven/tests/datasets/densityEstimation/9D_KurB1B2StroSkewSepBiClawB3B4Tri.csv";

  // Create config file
  std::string config = "tmpsgdeconfig.json";
  std::ofstream stream(config);
  stream << "{" << "\"dataSource\" : { \"filePath\" : \"" << samples << "\"},\"scorer\" : "
      << "{ \"testing\" : { \"metric\" : \"MSE\"}},\"fitter\" : " <<
      "{ \"type\" : \"densityEstimation\", \"gridConfig\" : { \"gridType\" : \"linear\","
      << "\"level\" : 4},\"adaptivityConfig\" : {\"numRefinements\" : 3, \"threshold\" : 0.001,"
      << "\"maxLevelType\" : false, \"noPoints\" : 3},\"regularizationConfig\" : {\"lambda\" : "
      << "1}}}" << std::endl;

  double mse = testDistribution(
      "datadriven/tests/datasets/densityEstimation/9D_KurB1B2StroSkewSepBiClawB3B4TriF.csv",
      config);
  std::cout << "MSE on test " << mse << std::endl;
  BOOST_CHECK(mse <= 2e-1);
  remove(config.c_str());
}


BOOST_AUTO_TEST_CASE(Test_10D_B4StroSkewSepBiB2TriB1B3ClawKurG) {
  std::string samples =
      "datadriven/tests/datasets/densityEstimation/10D_B4StroSkewSepBiB2TriB1B3ClawKurG.csv";

  // Create config file
  std::string config = "tmpsgdeconfig.json";
  std::ofstream stream(config);
  stream << "{" << "\"dataSource\" : { \"filePath\" : \"" << samples << "\"},\"scorer\" : "
      << "{ \"testing\" : { \"metric\" : \"MSE\"}},\"fitter\" : " <<
      "{ \"type\" : \"densityEstimation\", \"gridConfig\" : { \"gridType\" : \"linear\","
      << "\"level\" : 3},\"adaptivityConfig\" : {\"numRefinements\" : 3, \"threshold\" : 0.001,"
      << "\"maxLevelType\" : false, \"noPoints\" : 3},\"regularizationConfig\" : {\"lambda\" : "
      << "1}}}" << std::endl;

  double mse = testDistribution(
      "datadriven/tests/datasets/densityEstimation/10D_B4StroSkewSepBiB2TriB1B3ClawKurGF.csv",
      config);
  std::cout << "MSE on test " << mse << std::endl;
  BOOST_CHECK(mse <= 3e-1);
  remove(config.c_str());
}

BOOST_AUTO_TEST_SUITE_END()


