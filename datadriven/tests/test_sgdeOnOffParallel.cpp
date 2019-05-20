/*
 * Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * test_sgdeOnOffParallel.cpp
 *
 * Created on: Apr 19, 2019
 *     Author: Jan Schopohl
 */
#ifdef USE_SCALAPACK

#include <boost/test/test_tools.hpp>
#include <boost/test/unit_test_suite.hpp>
#include <sgpp/base/datatypes/DataMatrix.hpp>
#include <sgpp/base/datatypes/DataVector.hpp>
#include <sgpp/datadriven/datamining/base/SparseGridMiner.hpp>
#include <sgpp/datadriven/datamining/builder/DensityEstimationMinerFactory.hpp>
#include <sgpp/datadriven/datamining/modules/dataSource/CSVFileSampleProvider.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/FitterConfigurationDensityEstimation.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingBase.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingDensityEstimationOnOffParallel.hpp>
#include <sgpp/datadriven/scalapack/BlacsProcessGrid.hpp>

#include <mpi.h>
#include <cmath>
#include <fstream>
#include <iostream>
#include <random>
#include <string>

using sgpp::datadriven::BlacsProcessGrid;
using sgpp::datadriven::CSVFileSampleProvider;
using sgpp::datadriven::DataVector;
using sgpp::datadriven::DensityEstimationMinerFactory;
using sgpp::datadriven::ModelFittingBase;
using sgpp::datadriven::ModelFittingDensityEstimationOnOffParallel;
using sgpp::datadriven::SparseGridMiner;

// adapted sgde test for ScaLAPACK version

double testDistributionOnOffParallel(std::string testCSV, std::string config, bool parallel) {
  MPI_Barrier(MPI_COMM_WORLD);
  // Train model
  DensityEstimationMinerFactory factory;
  auto miner = std::unique_ptr<SparseGridMiner>(factory.buildMiner(config));
  miner->learn(false);
  ModelFittingBase *model = miner->getModel();

  if (parallel) {
    auto modelParallel = dynamic_cast<ModelFittingDensityEstimationOnOffParallel *>(model);
    BOOST_CHECK(modelParallel);
  }

  // Test
  CSVFileSampleProvider csv;
  csv.readFile(testCSV, true);
  auto testDataset = *(csv.getAllSamples());
  DataVector predictions(testDataset.getNumberInstances());
  model->evaluate(testDataset.getData(), predictions);

  if (BlacsProcessGrid::getCurrentProcess() == 0) {
    // Calculate MSE
    predictions.sub(testDataset.getTargets());
    double mse = predictions.l2Norm() / static_cast<double>(predictions.getSize());
    std::cout << "mse: " << mse << std::endl;
    return mse;
  }
  return 1.0;
}

BOOST_AUTO_TEST_SUITE(testDensityEstimationOnOffParallel)

BOOST_AUTO_TEST_CASE(Test_2D_StroSkewB2Parallel) {
  double tolerance = 5e-2;
  std::string samples = "datadriven/datasets/densityEstimation/2D_StroSkewB2.csv";

  // Parallel Orthoadapt
  std::string config = "tmpsgdeconfig.json";
  if (BlacsProcessGrid::getCurrentProcess() == 0) {
    std::ofstream stream0(config);
    stream0
        << "{"
        << "\"dataSource\" : { \"filePath\" : \"" << samples
        << "\", \"hasTargets\" : false},\"scorer\" : "
        << "{ \"metric\" : \"NLL\"},\"fitter\" : "
        << "{ \"type\" : \"densityEstimation\", \"gridConfig\" : { \"gridType\" : \"linear\","
        << "\"level\" : 5},\"adaptivityConfig\" : {\"numRefinements\" : 3, \"threshold\" : 0.001,"
        << "\"maxLevelType\" : false, \"noPoints\" : 3},\"regularizationConfig\" : {\"lambda\" : "
        << "1}, \"densityEstimationConfig\" : { \"densityEstimationType\" : \"decomposition\", "
        << "\"matrixDecompositionType\": \"orthoadapt\"}, \"parallelConfig\": { \"rowBlockSize\": "
        << "32, \"columnBlockSize\": 32}}}" << std::endl;
    stream0.flush();
  }

  double mseParallelOrtho = testDistributionOnOffParallel(
      "datadriven/datasets/densityEstimation/2D_StroSkewB2F.csv", config, true);

  if (BlacsProcessGrid::getCurrentProcess() == 0) {
    remove(config.c_str());
    // Serial Orthoadapt
    std::ofstream stream1(config);
    stream1
        << "{"
        << "\"dataSource\" : { \"filePath\" : \"" << samples
        << "\", \"hasTargets\" : false},\"scorer\" : "
        << "{ \"metric\" : \"NLL\"},\"fitter\" : "
        << "{ \"type\" : \"densityEstimation\", \"gridConfig\" : { \"gridType\" : \"linear\","
        << "\"level\" : 5},\"adaptivityConfig\" : {\"numRefinements\" : 3, \"threshold\" : 0.001,"
        << "\"maxLevelType\" : false, \"noPoints\" : 3},\"regularizationConfig\" : {\"lambda\" : "
        << "1}, \"densityEstimationConfig\" : { \"densityEstimationType\" : \"decomposition\", "
        << "\"matrixDecompositionType\": \"orthoadapt\"}}}" << std::endl;
    stream1.flush();
  }
  double mseOrtho = testDistributionOnOffParallel(
      "datadriven/datasets/densityEstimation/2D_StroSkewB2F.csv", config, false);

  if (BlacsProcessGrid::getCurrentProcess() == 0) {
    remove(config.c_str());
    // Parallel Cholesky
    std::ofstream stream2(config);
    stream2
        << "{"
        << "\"dataSource\" : { \"filePath\" : \"" << samples
        << "\", \"hasTargets\" : false},\"scorer\" : "
        << "{ \"metric\" : \"NLL\"},\"fitter\" : "
        << "{ \"type\" : \"densityEstimation\", \"gridConfig\" : { \"gridType\" : \"linear\","
        << "\"level\" : 5},\"adaptivityConfig\" : {\"numRefinements\" : 3, \"threshold\" : 0.001,"
        << "\"maxLevelType\" : false, \"noPoints\" : 3},\"regularizationConfig\" : {\"lambda\" : "
        << "1}, \"densityEstimationConfig\" : { \"densityEstimationType\" : \"decomposition\", "
        << "\"matrixDecompositionType\": \"chol\"}, \"parallelConfig\": { \"rowBlockSize\": "
        << "32, \"columnBlockSize\": 32}}}" << std::endl;
    stream2.flush();
  }
  double mseParallelChol = testDistributionOnOffParallel(
      "datadriven/datasets/densityEstimation/2D_StroSkewB2F.csv", config, true);

  if (BlacsProcessGrid::getCurrentProcess() == 0) {
    remove(config.c_str());
    // Serial Cholesky
    std::ofstream stream3(config);
    stream3
        << "{"
        << "\"dataSource\" : { \"filePath\" : \"" << samples
        << "\", \"hasTargets\" : false},\"scorer\" : "
        << "{ \"metric\" : \"NLL\"},\"fitter\" : "
        << "{ \"type\" : \"densityEstimation\", \"gridConfig\" : { \"gridType\" : \"linear\","
        << "\"level\" : 5},\"adaptivityConfig\" : {\"numRefinements\" : 3, \"threshold\" : 0.001,"
        << "\"maxLevelType\" : false, \"noPoints\" : 3},\"regularizationConfig\" : {\"lambda\" : "
        << "1}, \"densityEstimationConfig\" : { \"densityEstimationType\" : \"decomposition\", "
        << "\"matrixDecompositionType\": \"chol\"}}}" << std::endl;
    stream3.flush();
  }
  double mseChol = testDistributionOnOffParallel(
      "datadriven/datasets/densityEstimation/2D_StroSkewB2F.csv", config, false);

  if (BlacsProcessGrid::getCurrentProcess() == 0) {
    BOOST_CHECK(mseParallelOrtho <= tolerance);
    BOOST_CHECK(mseParallelChol <= tolerance);
    BOOST_CHECK_CLOSE(mseParallelOrtho, mseOrtho, tolerance);
    BOOST_CHECK_CLOSE(mseParallelChol, mseChol, tolerance);
    remove(config.c_str());
  }
}

BOOST_AUTO_TEST_CASE(Test_3D_KurB4B1Parallel) {
  double tolerance = 5e-2;
  std::string samples = "datadriven/datasets/densityEstimation/3D_KurB4B1.csv";

  // Parallel Orthoadapt
  std::string config = "tmpsgdeconfig.json";
  if (BlacsProcessGrid::getCurrentProcess() == 0) {
    std::ofstream stream0(config);
    stream0
        << "{"
        << "\"dataSource\" : { \"filePath\" : \"" << samples
        << "\", \"hasTargets\" : false},\"scorer\" : "
        << "{ \"metric\" : \"NLL\"},\"fitter\" : "
        << "{ \"type\" : \"densityEstimation\", \"gridConfig\" : { \"gridType\" : \"linear\","
        << "\"level\" : 5},\"adaptivityConfig\" : {\"numRefinements\" : 3, \"threshold\" : 0.001,"
        << "\"maxLevelType\" : false, \"noPoints\" : 3},\"regularizationConfig\" : {\"lambda\" : "
        << "1}, \"densityEstimationConfig\" : { \"densityEstimationType\" : \"decomposition\", "
        << "\"matrixDecompositionType\": \"orthoadapt\"}, \"parallelConfig\": { \"rowBlockSize\": "
        << "32, \"columnBlockSize\": 32}}}" << std::endl;
    stream0.flush();
  }
  double mseParallelOrtho = testDistributionOnOffParallel(
      "datadriven/datasets/densityEstimation/3D_KurB4B1F.csv", config, true);

  if (BlacsProcessGrid::getCurrentProcess() == 0) {
    remove(config.c_str());
    // Serial Orthoadapt
    std::ofstream stream1(config);
    stream1
        << "{"
        << "\"dataSource\" : { \"filePath\" : \"" << samples
        << "\", \"hasTargets\" : false},\"scorer\" : "
        << "{ \"metric\" : \"NLL\"},\"fitter\" : "
        << "{ \"type\" : \"densityEstimation\", \"gridConfig\" : { \"gridType\" : \"linear\","
        << "\"level\" : 5},\"adaptivityConfig\" : {\"numRefinements\" : 3, \"threshold\" : 0.001,"
        << "\"maxLevelType\" : false, \"noPoints\" : 3},\"regularizationConfig\" : {\"lambda\" : "
        << "1}, \"densityEstimationConfig\" : { \"densityEstimationType\" : \"decomposition\", "
        << "\"matrixDecompositionType\": \"orthoadapt\"}}}" << std::endl;
    stream1.flush();
  }
  double mseOrtho = testDistributionOnOffParallel(
      "datadriven/datasets/densityEstimation/3D_KurB4B1F.csv", config, false);

  if (BlacsProcessGrid::getCurrentProcess() == 0) {
    remove(config.c_str());
    // Parallel Cholesky
    std::ofstream stream2(config);
    stream2
        << "{"
        << "\"dataSource\" : { \"filePath\" : \"" << samples
        << "\", \"hasTargets\" : false},\"scorer\" : "
        << "{ \"metric\" : \"NLL\"},\"fitter\" : "
        << "{ \"type\" : \"densityEstimation\", \"gridConfig\" : { \"gridType\" : \"linear\","
        << "\"level\" : 5},\"adaptivityConfig\" : {\"numRefinements\" : 3, \"threshold\" : 0.001,"
        << "\"maxLevelType\" : false, \"noPoints\" : 3},\"regularizationConfig\" : {\"lambda\" : "
        << "1}, \"densityEstimationConfig\" : { \"densityEstimationType\" : \"decomposition\", "
        << "\"matrixDecompositionType\": \"chol\"}, \"parallelConfig\": { \"rowBlockSize\": "
        << "32, \"columnBlockSize\": 32}}}" << std::endl;
    stream2.flush();
  }
  double mseParallelChol = testDistributionOnOffParallel(
      "datadriven/datasets/densityEstimation/3D_KurB4B1F.csv", config, true);

  if (BlacsProcessGrid::getCurrentProcess() == 0) {
    remove(config.c_str());
    // Serial Cholesky
    std::ofstream stream3(config);
    stream3
        << "{"
        << "\"dataSource\" : { \"filePath\" : \"" << samples
        << "\", \"hasTargets\" : false},\"scorer\" : "
        << "{ \"metric\" : \"NLL\"},\"fitter\" : "
        << "{ \"type\" : \"densityEstimation\", \"gridConfig\" : { \"gridType\" : \"linear\","
        << "\"level\" : 5},\"adaptivityConfig\" : {\"numRefinements\" : 3, \"threshold\" : 0.001,"
        << "\"maxLevelType\" : false, \"noPoints\" : 3},\"regularizationConfig\" : {\"lambda\" : "
        << "1}, \"densityEstimationConfig\" : { \"densityEstimationType\" : \"decomposition\", "
        << "\"matrixDecompositionType\": \"chol\"}}}" << std::endl;
    stream3.flush();
  }
  double mseChol = testDistributionOnOffParallel(
      "datadriven/datasets/densityEstimation/3D_KurB4B1F.csv", config, false);

  if (BlacsProcessGrid::getCurrentProcess() == 0) {
    BOOST_CHECK(mseParallelOrtho <= tolerance);
    BOOST_CHECK(mseParallelChol <= tolerance);
    BOOST_CHECK_CLOSE(mseParallelOrtho, mseOrtho, tolerance);
    BOOST_CHECK_CLOSE(mseParallelChol, mseChol, tolerance);
    remove(config.c_str());
  }
}

BOOST_AUTO_TEST_SUITE_END()

#endif /* USE_SCALAPACK */
