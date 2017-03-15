/* Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * dataMiningConfigParserTest.cpp
 *
 *  Created on: 06.10.2016
 *      Author: Michael Lettrich
 */

#define BOOST_TEST_DYN_LINK

#include <boost/test/unit_test.hpp>
#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/datadriven/application/RegularizationConfiguration.hpp>
#include <sgpp/datadriven/datamining/configuration/DataMiningConfigParser.hpp>
#include <sgpp/datadriven/datamining/modules/dataSource/DataSourceConfig.hpp>
#include <sgpp/datadriven/datamining/modules/scoring/ScorerConfig.hpp>
#include <sgpp/solver/TypesSolver.hpp>
#include <string>

const auto datasetPath = "datadriven/tests/data/dataminingConfig.json";

BOOST_AUTO_TEST_SUITE(dataMiningConfigParserTest)

using sgpp::datadriven::DataMiningConfigParser;
using sgpp::datadriven::DataSourceConfig;
using sgpp::datadriven::DataSourceFileType;
using sgpp::datadriven::TestingConfiguration;
using sgpp::datadriven::CrossValidationConfiguration;
using sgpp::datadriven::RegularizationConfiguration;
using sgpp::datadriven::RegularizationType;
using sgpp::datadriven::ScorerShufflingType;
using sgpp::datadriven::ScorerMetricType;
using sgpp::datadriven::FitterType;
using sgpp::base::RegularGridConfiguration;
using sgpp::base::GridType;
using sgpp::solver::SLESolverConfiguration;
using sgpp::solver::SLESolverType;

BOOST_AUTO_TEST_CASE(testTopLevel) {
  DataMiningConfigParser parser{datasetPath};
  BOOST_CHECK_EQUAL(parser.hasDataSourceConfig(), true);
  BOOST_CHECK_EQUAL(parser.hasFitterConfig(), true);
  BOOST_CHECK_EQUAL(parser.hasScorerConfig(), true);
}

BOOST_AUTO_TEST_CASE(testDataSourceConfig) {
  DataMiningConfigParser parser{datasetPath};

  DataSourceConfig defaults;
  defaults.filePath = "something/false";
  defaults.fileType = DataSourceFileType::NONE;
  defaults.isCompressed = true;
  defaults.numBatches = 2;
  defaults.batchSize = 10;
  DataSourceConfig config;
  bool hasConfig;

  hasConfig = parser.getDataSourceConfig(config, defaults);

  BOOST_CHECK_EQUAL(hasConfig, true);
  BOOST_CHECK_EQUAL(std::strcmp(config.filePath.c_str(), "/path/to/some/file.arff"), 0);
  BOOST_CHECK_EQUAL(static_cast<int>(config.fileType), static_cast<int>(DataSourceFileType::ARFF));
  BOOST_CHECK_EQUAL(config.numBatches, 1);
  BOOST_CHECK_EQUAL(config.batchSize, 0);
}

BOOST_AUTO_TEST_CASE(testScorerTestingConfig) {
  DataMiningConfigParser parser{datasetPath};

  TestingConfiguration defaults;
  defaults.testingPortion = 1.0;
  defaults.shuffling = ScorerShufflingType::sequential;
  defaults.randomSeed = 40;
  defaults.metric = ScorerMetricType::mse;
  TestingConfiguration config;
  bool hasConfig;
  double tolerance = 1E-5;

  hasConfig = parser.getScorerTestingConfig(config, defaults);

  BOOST_CHECK_EQUAL(hasConfig, true);
  BOOST_CHECK_CLOSE(config.testingPortion, 0.0, tolerance);
  BOOST_CHECK_EQUAL(static_cast<int>(config.shuffling),
                    static_cast<int>(ScorerShufflingType::random));
  BOOST_CHECK_EQUAL(config.randomSeed, 42);
  BOOST_CHECK_EQUAL(static_cast<int>(config.metric), static_cast<int>(ScorerMetricType::mse));
}

BOOST_AUTO_TEST_CASE(testScorerCrossValidationConfig) {
  DataMiningConfigParser parser{datasetPath};

  CrossValidationConfiguration defaults;
  defaults.folds = 3;
  defaults.shuffling = ScorerShufflingType::sequential;
  defaults.randomSeed = 40;
  defaults.metric = ScorerMetricType::mse;
  CrossValidationConfiguration config;
  bool hasConfig;

  hasConfig = parser.getScorerCrossValidationConfig(config, defaults);

  BOOST_CHECK_EQUAL(hasConfig, true);
  BOOST_CHECK_EQUAL(config.folds, 5);
  BOOST_CHECK_EQUAL(static_cast<int>(config.shuffling),
                    static_cast<int>(ScorerShufflingType::random));
  BOOST_CHECK_EQUAL(config.randomSeed, 42);
  BOOST_CHECK_EQUAL(static_cast<int>(config.metric), static_cast<int>(ScorerMetricType::mse));
}

BOOST_AUTO_TEST_CASE(testFitterTypeConfig) {
  DataMiningConfigParser parser{datasetPath};

  FitterType defaults = FitterType::RegressionLeastSquares;
  FitterType config;
  bool hasConfig;

  hasConfig = parser.getFitterConfigType(config, defaults);

  BOOST_CHECK_EQUAL(hasConfig, true);
  BOOST_CHECK_EQUAL(static_cast<int>(config), static_cast<int>(FitterType::RegressionLeastSquares));
}

BOOST_AUTO_TEST_CASE(testFitterGridConfig) {
  DataMiningConfigParser parser{datasetPath};

  RegularGridConfiguration defaults;
  defaults.type_ = GridType::ModLinear;
  defaults.dim_ = 2;
  defaults.level_ = 4;
  defaults.maxDegree_ = 10;
  defaults.boundaryLevel_ = 5;
  defaults.filename_ = "something";
  RegularGridConfiguration config;
  bool hasConfig;

  hasConfig = parser.getFitterGridConfig(config, defaults);

  BOOST_CHECK_EQUAL(hasConfig, true);
  BOOST_CHECK_EQUAL(static_cast<int>(config.type_), static_cast<int>(GridType::Linear));
  BOOST_CHECK_EQUAL(config.dim_, 0);
  BOOST_CHECK_EQUAL(config.level_, 2);
  BOOST_CHECK_EQUAL(config.maxDegree_, 0);
  BOOST_CHECK_EQUAL(config.boundaryLevel_, 0);
  BOOST_CHECK_EQUAL(std::strcmp(config.filename_.c_str(), ""), 0);
}

BOOST_AUTO_TEST_CASE(testFitterAdaptivityConfig) {
  DataMiningConfigParser parser{datasetPath};

  AdpativityConfiguration defaults;
  defaults.numRefinements_ = 42;
  defaults.threshold_ = 42;
  defaults.maxLevelType_ = true;
  defaults.noPoints_ = 42;
  defaults.percent_ = 0.42;
  defaults.errorBasedRefinement = true;
  AdpativityConfiguration config;
  bool hasConfig;
  double tolerance = 1E-5;

  hasConfig = parser.getFitterAdaptivityConfig(config, defaults);

  BOOST_CHECK_EQUAL(hasConfig, true);
  BOOST_CHECK_EQUAL(config.numRefinements_, 0);
  BOOST_CHECK_EQUAL(config.threshold_, 0);
  BOOST_CHECK_EQUAL(config.maxLevelType_, false);
  BOOST_CHECK_EQUAL(config.noPoints_, 0);
  BOOST_CHECK_CLOSE(config.percent_, 0, tolerance);
  BOOST_CHECK_EQUAL(config.errorBasedRefinement, false);
}

BOOST_AUTO_TEST_CASE(testFitterSolverRefineConfig) {
  DataMiningConfigParser parser{datasetPath};

  SLESolverConfiguration defaults;
  defaults.type_ = SLESolverType::BiCGSTAB;
  defaults.eps_ = 10e-5;
  defaults.maxIterations_ = 42;
  defaults.threshold_ = 1;
  SLESolverConfiguration config;
  bool hasConfig;
  double tolerance = 1E-5;

  hasConfig = parser.getFitterSolverRefineConfig(config, defaults);

  BOOST_CHECK_EQUAL(hasConfig, true);
  BOOST_CHECK_EQUAL(static_cast<int>(config.type_), static_cast<int>(SLESolverType::CG));
  BOOST_CHECK_CLOSE(config.eps_, 10e-15, tolerance);
  BOOST_CHECK_EQUAL(config.maxIterations_, 100);
  BOOST_CHECK_EQUAL(config.threshold_, 1);
}

BOOST_AUTO_TEST_CASE(testFitterSolverFinalConfig) {
  DataMiningConfigParser parser{datasetPath};

  SLESolverConfiguration defaults;
  defaults.type_ = SLESolverType::BiCGSTAB;
  defaults.eps_ = 10e-5;
  defaults.maxIterations_ = 42;
  defaults.threshold_ = 1;
  SLESolverConfiguration config;
  bool hasConfig;
  double tolerance = 1E-5;

  hasConfig = parser.getFitterSolverFinalConfig(config, defaults);

  BOOST_CHECK_EQUAL(hasConfig, true);
  BOOST_CHECK_EQUAL(static_cast<int>(config.type_), static_cast<int>(SLESolverType::CG));
  BOOST_CHECK_CLOSE(config.eps_, 10e-15, tolerance);
  BOOST_CHECK_EQUAL(config.maxIterations_, 100);
  BOOST_CHECK_EQUAL(config.threshold_, 1);
}

BOOST_AUTO_TEST_CASE(testFitterRegularizationConfig) {
  DataMiningConfigParser parser{datasetPath};

  RegularizationConfiguration defaults;
  defaults.regType_ = RegularizationType::Laplace;
  defaults.lambda_ = 1;
  defaults.exponentBase_ = 2;
  defaults.l1Ratio_ = 3;
  RegularizationConfiguration config;
  bool hasConfig;
  double tolerance = 1E-5;

  hasConfig = parser.getFitterRegularizationConfig(config, defaults);

  BOOST_CHECK_EQUAL(hasConfig, true);
  BOOST_CHECK_EQUAL(static_cast<int>(config.regType_),
                    static_cast<int>(RegularizationType::Identity));
  BOOST_CHECK_CLOSE(config.lambda_, 10e-7, tolerance);
  BOOST_CHECK_CLOSE(config.exponentBase_, 3.0, tolerance);
  BOOST_CHECK_CLOSE(config.l1Ratio_, 4.0, tolerance);
}

BOOST_AUTO_TEST_SUITE_END()
