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

double testDistributionOnOff(std::string testCSV, std::string config) {
  // Train model
  DensityEstimationMinerFactory factory;
  auto miner = std::unique_ptr<SparseGridMiner>(factory.buildMiner(config));
  miner->learn(false);
  ModelFittingBase *model = miner->getModel();
  // Test
  CSVFileSampleProvider csv;
  csv.readFile(testCSV, true);
  auto testDataset = *(csv.getAllSamples());
  DataVector predictions(testDataset.getNumberInstances());
  model->evaluate(testDataset.getData(), predictions);
  // Calculate MSE
  predictions.sub(testDataset.getTargets());
  return predictions.l2Norm() / static_cast<double>(predictions.getSize());
}

BOOST_AUTO_TEST_SUITE(testDensityEstimationOnOff)


BOOST_AUTO_TEST_CASE(Test_2D_StroSkewB2) {
  std::string samples = "datadriven/datasets/densityEstimation/2D_StroSkewB2.csv";

  // Create config file
  std::string config = "tmpsgdeconfig.json";
  std::ofstream stream(config);
  stream << "{" << "\"dataSource\" : { \"filePath\" : \"" << samples <<
      "\", \"hasTargets\" : false},\"scorer\" : "
      << "{ \"metric\" : \"NLL\"},\"fitter\" : " <<
      "{ \"type\" : \"densityEstimation\", \"gridConfig\" : { \"gridType\" : \"linear\","
      << "\"level\" : 5},\"adaptivityConfig\" : {\"numRefinements\" : 3, \"threshold\" : 0.001,"
      << "\"maxLevelType\" : false, \"noPoints\" : 3},\"regularizationConfig\" : {\"lambda\" : "
      << "1}, \"densityEstimationConfig\" : { \"densityEstimationType\" : \"decomposition\"}}}"
      << std::endl;

  double mse = testDistributionOnOff(
      "datadriven/datasets/densityEstimation/2D_StroSkewB2F.csv", config);
  std::cout << "MSE on test " << mse << std::endl;
  BOOST_CHECK(mse <= 5e-2);
  remove(config.c_str());
}


BOOST_AUTO_TEST_CASE(Test_3D_KurB4B1) {
  std::string samples = "datadriven/datasets/densityEstimation/3D_KurB4B1.csv";

  // Create config file
  std::string config = "tmpsgdeconfig.json";
  std::ofstream stream(config);
  stream << "{" << "\"dataSource\" : { \"filePath\" : \"" << samples <<
        "\", \"hasTargets\" : false},\"scorer\" : "
        << "{ \"metric\" : \"NLL\"},\"fitter\" : " <<
        "{ \"type\" : \"densityEstimation\", \"gridConfig\" : { \"gridType\" : \"linear\","
        << "\"level\" : 5},\"adaptivityConfig\" : {\"numRefinements\" : 3, \"threshold\" : 0.001,"
        << "\"maxLevelType\" : false, \"noPoints\" : 3},\"regularizationConfig\" : {\"lambda\" : "
        << "1}, \"densityEstimationConfig\" : { \"densityEstimationType\" : \"decomposition\"}}}"
        << std::endl;

  double mse = testDistributionOnOff(
      "datadrive/datasets/densityEstimation/3D_KurB4B1F.csv", config);
  std::cout << "MSE between estimation and ground truth density " << mse << std::endl;
  BOOST_CHECK(mse <= 5e-2);
  remove(config.c_str());
}



BOOST_AUTO_TEST_SUITE_END()

#endif /* USE_GSL */
