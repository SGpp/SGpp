// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifdef USE_SCALAPACK

#define BOOST_TEST_DYN_LINK
#include <boost/test/test_tools.hpp>
#include <boost/test/unit_test_suite.hpp>
#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/datadriven/datamining/base/SparseGridMiner.hpp>
#include <sgpp/datadriven/datamining/builder/DensityDerivativeEstimationMinerFactory.hpp>
#include <sgpp/datadriven/datamining/modules/dataSource/CSVFileSampleProvider.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/FitterConfigurationDensityEstimation.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingBase.hpp>

#include <cmath>
#include <fstream>
#include <iostream>
#include <random>
#include <string>

using sgpp::datadriven::CSVFileSampleProvider;
using sgpp::datadriven::DataVector;
using sgpp::datadriven::DensityDerivativeEstimationMinerFactory;
using sgpp::datadriven::ModelFittingBase;
using sgpp::datadriven::SparseGridMiner;

double testDistributionDDerivE_OnOffParallel(std::string testCSV, std::string config) {
  // Train model
  DensityDerivativeEstimationMinerFactory factory;
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

BOOST_AUTO_TEST_SUITE(testDensityDerivativeEstimationOnOffParallel)

BOOST_AUTO_TEST_CASE(Test_2D_DDeriv_Gauss) {
  std::string samples = "datadriven/datasets/densityEstimation/2D_DDeriv_Gauss_test.csv";

  // Create config file
  std::string config = "tmpsgdeconfig.json";
  std::ofstream stream(config);
  stream << "{ \"dataSource\" : { \"filePath\" : \"" << samples
         << "\", \"hasTargets\" : false}, \"scorer\" : { \"metric\" : \"NLL\" }, "
         << "\"fitter\" : { \"type\" : \"densityDerivativeEstimation\", "
         << "\"gridConfig\" : { \"gridType\" : \"bspline\", \"level\" : 3, \"maxDegree\" : 3 }, "
         << "\"adaptivityConfig\" : { \"numRefinements\" : 3, \"threshold\" : 0.001, "
         << "\"maxLevelType\" : false, \"noPoints\" : 3 }, "
         << "\"regularizationConfig\" : { \"lambda\" : 0.001 }, "
         << "\"densityEstimationConfig\" : { \"derivDim\" : 0, "
            "\"densityEstimationType\" : \"decomposition\" }, \"parallelConfig\" : {} } }"
         << std::endl;

  double mse = testDistributionDDerivE_OnOffParallel(
      "datadriven/datasets/densityEstimation/2D_DDeriv_Gauss_val.csv", config);
  if (sgpp::datadriven::BlacsProcessGrid::getCurrentProcess() == 0) {
    std::cout << "MSE on test " << mse << std::endl;
    BOOST_CHECK(mse <= 5e-2);
  }
  remove(config.c_str());
}

BOOST_AUTO_TEST_SUITE_END()

#endif /* USE_SCALAPACK */
