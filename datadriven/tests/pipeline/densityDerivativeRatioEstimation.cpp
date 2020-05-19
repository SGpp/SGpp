// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifdef USE_GSL

#include <boost/test/test_tools.hpp>
#include <boost/test/unit_test_suite.hpp>
#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/datadriven/datamining/base/SparseGridMiner.hpp>
#include <sgpp/datadriven/datamining/builder/DensityDerivativeRatioEstimationMinerFactory.hpp>
#include <sgpp/datadriven/datamining/modules/dataSource/CSVFileSampleProvider.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingBase.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/FitterConfigurationDensityEstimation.hpp>

#include <cmath>
#include <fstream>
#include <iostream>
#include <random>
#include <string>

using sgpp::datadriven::DensityDerivativeRatioEstimationMinerFactory;
using sgpp::datadriven::SparseGridMiner;
using sgpp::datadriven::ModelFittingBase;
using sgpp::datadriven::CSVFileSampleProvider;
using sgpp::datadriven::DataVector;

double testDistributionDDerivRE(std::string testCSV, std::string config) {
  // Train model
  DensityDerivativeRatioEstimationMinerFactory factory;
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

BOOST_AUTO_TEST_SUITE(testDensityDerivativeRatioEstimation)

BOOST_AUTO_TEST_CASE(Test_2D_Unif) {
  std::string samples = "datadriven/datasets/densityEstimation/2D_DDeriv_Unif_test.csv";

  // Create config file for direct density derivative ratio estimation
  std::string config = "tmpsgdeconfigDRE.json";
  std::ofstream stream(config);
  stream << "{ \"dataSource\": { \"filePath\" : \"" << samples << "\", \"hasTargets\" : false}, "
         << "\"scorer\" : { \"metric\" : \"NLL\" }, "
         << "\"fitter\" : { \"type\" : \"densityDerivativeRatioEstimation\", "
         << "\"gridConfig\" : { \"gridType\" : \"bspline\", \"level\" : 5, \"maxDegree\" : 3 }, "
         << "\"adaptivityConfig\" : { \"numRefinements\" : 3, \"threshold\" : 0.001, "
            "\"maxLevelType\" : false, \"noPoints\" : 3 }, "
         << "\"regularizationConfig\" : { \"lambda\" : 1 }, "
         << "\"densityEstimationConfig\" : { \"derivDim\" : 0 } } }" << std::endl;

  // Reuse existing dataset for computing MSE
  double mse = testDistributionDDerivRE(
      "datadriven/datasets/densityEstimation/2D_DDE_B2-64_B2-2_val.csv", config);
  std::cout << "MSE between estimation and ground truth density " << mse << std::endl;
  BOOST_CHECK(mse <= 5e-2);
  remove(config.c_str());
}

BOOST_AUTO_TEST_SUITE_END()

#endif /* USE_GSL */
