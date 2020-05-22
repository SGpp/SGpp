// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#ifdef USE_GSL

#define BOOST_TEST_DYN_LINK
#include <boost/test/test_tools.hpp>
#include <boost/test/unit_test_suite.hpp>
#include <sgpp/datadriven/algorithm/DBMatObjectStore.hpp>
#include <sgpp/datadriven/algorithm/DBMatOfflineOrthoAdapt.hpp>
#include <sgpp/datadriven/algorithm/DBMatOnlineDEFactory.hpp>
#include <sgpp/datadriven/algorithm/DBMatOnlineDEOrthoAdapt.hpp>
#include <sgpp/datadriven/algorithm/DBMatPermutationFactory.hpp>
#include <sgpp/datadriven/algorithm/GridFactory.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/ModelFittingDensityEstimationCombi.hpp>

#include <set>
#include <vector>

BOOST_AUTO_TEST_SUITE(PermutationBlowUpTest)

BOOST_AUTO_TEST_CASE(ComponentGridOrthoTest) {
  sgpp::base::GeneralGridConfiguration baseGridConfig;
  baseGridConfig.generalType_ = sgpp::base::GeneralGridType::ComponentGrid;
  baseGridConfig.type_ = sgpp::base::GridType::Linear;
  baseGridConfig.levelVector_ = std::vector<size_t>{3, 2, 2};
  baseGridConfig.dim_ = 3;

  sgpp::base::GeneralGridConfiguration desiredGridConfig;
  desiredGridConfig.generalType_ = sgpp::base::GeneralGridType::ComponentGrid;
  desiredGridConfig.type_ = sgpp::base::GridType::Linear;
  desiredGridConfig.levelVector_ = std::vector<size_t>{2, 3, 2, 1, 1};
  desiredGridConfig.dim_ = 5;

  sgpp::base::AdaptivityConfiguration adaptivityConfig;
  adaptivityConfig.numRefinements_ = 0;

  sgpp::datadriven::RegularizationConfiguration regularizationConfig;
  regularizationConfig.lambda_ = 0;

  sgpp::datadriven::DensityEstimationConfiguration densityEstimationConfig;
  densityEstimationConfig.type_ = sgpp::datadriven::DensityEstimationType::Decomposition;
  densityEstimationConfig.decomposition_ = sgpp::datadriven::MatrixDecompositionType::OrthoAdapt;

  // build grids
  std::cout << "Build Grid" << std::endl;
  sgpp::datadriven::GridFactory gridFactory;
  std::set<std::set<size_t>> interactions = std::set<std::set<size_t>>();
  std::unique_ptr<sgpp::base::Grid> desiredGrid{
      gridFactory.createGrid(desiredGridConfig, interactions)};

  // build offline objects
  std::cout << "Instanciate Offline" << std::endl;

  std::unique_ptr<sgpp::datadriven::DBMatOfflineOrthoAdapt> desiredOff{
      dynamic_cast<sgpp::datadriven::DBMatOfflineOrthoAdapt *>(
          sgpp::datadriven::DBMatOfflineFactory::buildOfflineObject(
              desiredGridConfig, adaptivityConfig, regularizationConfig, densityEstimationConfig))};

  // object store
  std::shared_ptr<sgpp::datadriven::DBMatObjectStore> store =
      std::make_shared<sgpp::datadriven::DBMatObjectStore>();

  // build and decompose
  std::cout << "Build and decompose" << std::endl;
  // sgpp::datadriven::DBMatBaseObjectStore store2(adaptivityConfig, regularizationConfig,
  // densityEstimationConfig);
  // Add base to store
  sgpp::datadriven::GeometryConfiguration geometryConfig;
  sgpp::datadriven::DBMatPermutationFactory factory(store);
  // put base object into store
  factory.getPermutedObject(baseGridConfig, geometryConfig, adaptivityConfig, regularizationConfig,
                            densityEstimationConfig);
  // build and decompose desired offline object
  desiredOff->buildMatrix(desiredGrid.get(), regularizationConfig);
  desiredOff->decomposeMatrix(regularizationConfig, densityEstimationConfig);
  // get permuted base object from factory
  std::unique_ptr<sgpp::datadriven::DBMatOfflineOrthoAdapt> permOff(
      dynamic_cast<sgpp::datadriven::DBMatOfflineOrthoAdapt *>(
          factory.getPermutedObject(desiredGridConfig, geometryConfig, adaptivityConfig,
                                    regularizationConfig, densityEstimationConfig)));

  // build online objects
  std::unique_ptr<sgpp::datadriven::DBMatOnlineDEOrthoAdapt> online1{
      dynamic_cast<sgpp::datadriven::DBMatOnlineDEOrthoAdapt *>(
          sgpp::datadriven::DBMatOnlineDEFactory::buildDBMatOnlineDE(
              *permOff, *desiredGrid, regularizationConfig.lambda_, 0,
              sgpp::datadriven::MatrixDecompositionType::OrthoAdapt))};
  std::unique_ptr<sgpp::datadriven::DBMatOnlineDEOrthoAdapt> online2{
      dynamic_cast<sgpp::datadriven::DBMatOnlineDEOrthoAdapt *>(
          sgpp::datadriven::DBMatOnlineDEFactory::buildDBMatOnlineDE(
              *desiredOff, *desiredGrid, regularizationConfig.lambda_, 0,
              sgpp::datadriven::MatrixDecompositionType::OrthoAdapt))};

  // Generate sample dataset
  std::cout << "Generate samples" << std::endl;
  sgpp::base::DataMatrix samples(100, desiredGridConfig.dim_);
  for (size_t i = 0; i < 100; i++) {
    sgpp::base::DataVector vec(desiredGridConfig.dim_);
    for (size_t j = 0; j < desiredGridConfig.dim_; j++) {
      vec.at(j) = static_cast<double>(std::rand()) / RAND_MAX;
    }
    samples.setRow(i, vec);
  }

  // alphas
  sgpp::base::DataVector alpha1(permOff->getGridSize());
  sgpp::base::DataVector alpha2(desiredOff->getGridSize());

  // compute alphas
  online1->computeDensityFunction(alpha1, samples, *desiredGrid, densityEstimationConfig, false);
  online2->computeDensityFunction(alpha2, samples, *desiredGrid, densityEstimationConfig, false);

  for (size_t i = 0; i < alpha1.getSize(); i++) {
    BOOST_CHECK(std::abs(alpha1[i] - alpha2[i]) / std::abs(alpha2[i]) < 0.001);
  }
}

BOOST_AUTO_TEST_CASE(FullCombiSchemeOrthoTest) {
  // config
  sgpp::base::GeneralGridConfiguration gridConfig;
  gridConfig.generalType_ = sgpp::base::GeneralGridType::ComponentGrid;
  gridConfig.type_ = sgpp::base::GridType::Linear;
  gridConfig.dim_ = 10;
  gridConfig.level_ = 3;

  sgpp::base::AdaptivityConfiguration adaptivityConfig;
  adaptivityConfig.numRefinements_ = 0;

  sgpp::datadriven::RegularizationConfiguration regularizationConfig;
  regularizationConfig.lambda_ = 0;

  sgpp::datadriven::DensityEstimationConfiguration densityEstimationConfig;
  densityEstimationConfig.type_ = sgpp::datadriven::DensityEstimationType::Decomposition;
  densityEstimationConfig.decomposition_ = sgpp::datadriven::MatrixDecompositionType::OrthoAdapt;

  sgpp::datadriven::FitterConfigurationDensityEstimation config;
  config.getGridConfig() = gridConfig;
  config.getRefinementConfig() = adaptivityConfig;
  config.getDensityEstimationConfig() = densityEstimationConfig;
  config.getRegularizationConfig() = regularizationConfig;

  // generate samples
  sgpp::base::DataMatrix samples(1000, gridConfig.dim_);
  for (size_t i = 0; i < 1000; i++) {
    sgpp::base::DataVector vec(gridConfig.dim_);
    for (size_t j = 0; j < gridConfig.dim_; j++) {
      vec.at(j) = static_cast<double>(std::rand()) / RAND_MAX;
    }
    samples.setRow(i, vec);
  }

  std::shared_ptr<sgpp::datadriven::DBMatObjectStore> store(
      new sgpp::datadriven::DBMatObjectStore());

  sgpp::datadriven::ModelFittingDensityEstimationCombi modelWithPerm(config, store);
  modelWithPerm.fit(samples);

  // turn off offline permutation
  config.getDensityEstimationConfig().useOfflinePermutation_ = false;
  sgpp::datadriven::ModelFittingDensityEstimationCombi modelConventional(config);
  modelConventional.fit(samples);

  // check result
  for (size_t i = 0; i < 1000; i++) {
    sgpp::base::DataVector p(gridConfig.dim_);
    for (size_t j = 0; j < gridConfig.dim_; j++) {
      p.at(j) = static_cast<double>(std::rand()) / RAND_MAX;
    }
    double check1 = modelConventional.evaluate(p);
    double check2 = modelWithPerm.evaluate(p);

    BOOST_CHECK(std::abs(check1 - check2) / std::abs(check1) < 0.001);
  }
}

BOOST_AUTO_TEST_SUITE_END()

#endif  // USE_GSL
