// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#define BOOST_TEST_DYN_LINK

#include <boost/test/unit_test.hpp>

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/datadriven/configuration/GeometryConfiguration.hpp>
#include <sgpp/datadriven/configuration/ParallelConfiguration.hpp>
#include <sgpp/datadriven/configuration/RegularizationConfiguration.hpp>
#include <sgpp/datadriven/datamining/configuration/DataMiningConfigParser.hpp>
#include <sgpp/datadriven/datamining/modules/dataSource/DataSourceConfig.hpp>
#include <sgpp/datadriven/datamining/modules/dataSource/DataTransformationConfig.hpp>
#include <sgpp/datadriven/datamining/modules/scoring/ScorerConfig.hpp>
#include <sgpp/datadriven/datamining/modules/visualization/VisualizationGeneralConfig.hpp>
#include <sgpp/datadriven/datamining/modules/visualization/VisualizationParameters.hpp>
#include <sgpp/solver/TypesSolver.hpp>

#include <string>
#include <vector>

const auto datasetPath = "datadriven/tests/dataminingConfig.json";

BOOST_AUTO_TEST_SUITE(dataMiningConfigParserTest)

using sgpp::base::GridType;
using sgpp::base::RegularGridConfiguration;
using sgpp::base::AdaptivityConfiguration;
using sgpp::datadriven::DataMiningConfigParser;
using sgpp::datadriven::DataSourceConfig;
using sgpp::datadriven::DataSourceFileType;
using sgpp::datadriven::DataSourceShufflingType;
using sgpp::datadriven::DataTransformationType;
using sgpp::datadriven::FitterType;
using sgpp::datadriven::ParallelConfiguration;
using sgpp::datadriven::RegularizationConfiguration;
using sgpp::datadriven::RegularizationType;
using sgpp::datadriven::ScorerConfiguration;
using sgpp::datadriven::ScorerMetricType;
using sgpp::solver::SLESolverConfiguration;
using sgpp::solver::SLESolverType;
using sgpp::datadriven::VisualizationGeneralConfig;
using sgpp::datadriven::VisualizationParameters;

using sgpp::datadriven::VisualizationFileType;

BOOST_AUTO_TEST_CASE(testTopLevel) {
  DataMiningConfigParser parser{datasetPath};
  BOOST_CHECK_EQUAL(parser.hasDataSourceConfig(), true);
  BOOST_CHECK_EQUAL(parser.hasFitterConfig(), true);
  BOOST_CHECK_EQUAL(parser.hasScorerConfig(), true);
  BOOST_CHECK_EQUAL(parser.hasVisualizationConfig(), true);
}

BOOST_AUTO_TEST_CASE(testDataSourceConfig) {
  DataMiningConfigParser parser{datasetPath};

  DataSourceConfig defaults;
  defaults.filePath = "something/false";
  defaults.fileType = DataSourceFileType::NONE;
  defaults.isCompressed = true;
  defaults.numBatches = 2;
  defaults.batchSize = 10;
  defaults.epochs = 1;
  defaults.shuffling = DataSourceShufflingType::sequential;
  defaults.validationPortion = 0.1;
  defaults.randomSeed = 1337;

  defaults.testFilePath = "something/testFalse";
  defaults.testFileType = DataSourceFileType::NONE;
  defaults.testNumBatches = 1;
  defaults.testBatchSize = 4;
  defaults.testIsCompressed = true;

  DataSourceConfig config;
  bool hasConfig;
  bool hasDataTransformationConfig;

  hasConfig = parser.getDataSourceConfig(config, defaults);
  hasDataTransformationConfig = parser.hasDataTransformationConfig();

  BOOST_CHECK_EQUAL(hasConfig, true);
  BOOST_CHECK_EQUAL(hasDataTransformationConfig, true);
  BOOST_CHECK_EQUAL(std::strcmp(config.filePath.c_str(), "/path/to/some/file.arff"), 0);
  BOOST_CHECK_EQUAL(static_cast<int>(config.fileType), static_cast<int>(DataSourceFileType::ARFF));
  BOOST_CHECK_EQUAL(config.numBatches, 1);
  BOOST_CHECK_EQUAL(config.batchSize, 0);
  BOOST_CHECK_EQUAL(static_cast<int>(config.dataTransformationConfig.type),
                    static_cast<int>(DataTransformationType::ROSENBLATT));
  BOOST_CHECK_EQUAL(config.dataTransformationConfig.rosenblattConfig.solverMaxIterations, 1000);
  BOOST_CHECK_EQUAL(config.validationPortion, 0.634);
  BOOST_CHECK_EQUAL(config.epochs, 12);
  BOOST_CHECK_EQUAL(static_cast<int>(config.shuffling),
                    static_cast<int>(DataSourceShufflingType::random));
  BOOST_CHECK_EQUAL(std::strcmp(config.testFilePath.c_str(), "/path/to/some/testFile.arff"), 0);
  BOOST_CHECK_EQUAL(static_cast<int>(config.testFileType),
                    static_cast<int>(DataSourceFileType::ARFF));
  BOOST_CHECK_EQUAL(config.testNumBatches, 2);
  BOOST_CHECK_EQUAL(config.testBatchSize, 16);
  BOOST_CHECK_EQUAL(config.testIsCompressed, false);
}

BOOST_AUTO_TEST_CASE(testScorerConfig) {
  DataMiningConfigParser parser{datasetPath};

  ScorerConfiguration defaults;
  defaults.metric = ScorerMetricType::nll;
  ScorerConfiguration config;
  bool hasConfig;

  hasConfig = parser.getScorerConfig(config, defaults);

  BOOST_CHECK_EQUAL(hasConfig, true);
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

  AdaptivityConfiguration defaults;
  defaults.numRefinements_ = 42;
  defaults.threshold_ = 42;
  defaults.maxLevelType_ = true;
  defaults.noPoints_ = 42;
  defaults.percent_ = 0.42;
  defaults.errorBasedRefinement = true;
  AdaptivityConfiguration config;
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
  defaults.eps_ = 1e-4;
  defaults.maxIterations_ = 42;
  defaults.threshold_ = 1;
  SLESolverConfiguration config;
  bool hasConfig;
  double tolerance = 1E-5;

  hasConfig = parser.getFitterSolverRefineConfig(config, defaults);

  BOOST_CHECK_EQUAL(hasConfig, true);
  BOOST_CHECK_EQUAL(static_cast<int>(config.type_), static_cast<int>(SLESolverType::CG));
  BOOST_CHECK_CLOSE(config.eps_, 1e-14, tolerance);
  BOOST_CHECK_EQUAL(config.maxIterations_, 100);
  BOOST_CHECK_EQUAL(config.threshold_, 1);
}

BOOST_AUTO_TEST_CASE(testFitterSolverFinalConfig) {
  DataMiningConfigParser parser{datasetPath};

  SLESolverConfiguration defaults;
  defaults.type_ = SLESolverType::BiCGSTAB;
  defaults.eps_ = 1e-4;
  defaults.maxIterations_ = 42;
  defaults.threshold_ = 1;
  SLESolverConfiguration config;
  bool hasConfig;
  double tolerance = 1E-5;

  hasConfig = parser.getFitterSolverFinalConfig(config, defaults);

  BOOST_CHECK_EQUAL(hasConfig, true);
  BOOST_CHECK_EQUAL(static_cast<int>(config.type_), static_cast<int>(SLESolverType::CG));
  BOOST_CHECK_CLOSE(config.eps_, 1e-14, tolerance);
  BOOST_CHECK_EQUAL(config.maxIterations_, 100);
  BOOST_CHECK_EQUAL(config.threshold_, 1);
}

BOOST_AUTO_TEST_CASE(testFitterRegularizationConfig) {
  DataMiningConfigParser parser{datasetPath};

  RegularizationConfiguration defaults;
  defaults.type_ = RegularizationType::Laplace;
  defaults.lambda_ = 1;
  defaults.exponentBase_ = 2;
  defaults.l1Ratio_ = 3;
  RegularizationConfiguration config;
  bool hasConfig;
  double tolerance = 1E-5;

  hasConfig = parser.getFitterRegularizationConfig(config, defaults);

  BOOST_CHECK_EQUAL(hasConfig, true);
  BOOST_CHECK_EQUAL(static_cast<int>(config.type_), static_cast<int>(RegularizationType::Identity));
  BOOST_CHECK_CLOSE(config.lambda_, 1e-6, tolerance);
  BOOST_CHECK_CLOSE(config.exponentBase_, 3.0, tolerance);
  BOOST_CHECK_CLOSE(config.l1Ratio_, 4.0, tolerance);
}

BOOST_AUTO_TEST_CASE(testFitterGeometryConfig) {
  DataMiningConfigParser parser{datasetPath};

  sgpp::datadriven::GeometryConfiguration defaults;
  defaults.dim = std::vector<std::vector<int64_t>>();
  sgpp::datadriven::GeometryConfiguration config;
  bool hasConfig;
  std::vector<std::vector<int64_t>> dim = {{1, 2}, {3, 4}};

  hasConfig = parser.getGeometryConfig(config, defaults);

  BOOST_CHECK_EQUAL(hasConfig, true);

  BOOST_CHECK_EQUAL(static_cast<int>(config.stencils.at(0).stencilType),
                    static_cast<int>(sgpp::datadriven::StencilType::DirectNeighbour));
  BOOST_CHECK_EQUAL(config.stencils[0].applyOnLayers.size(), 2);
  BOOST_CHECK_EQUAL(config.stencils[0].applyOnLayers[0], 0);
  BOOST_CHECK_EQUAL(config.stencils[0].applyOnLayers[1], 1);
  BOOST_CHECK_EQUAL(config.stencils[0].colorIndex, 0);

  BOOST_CHECK_EQUAL(static_cast<int>(config.stencils.at(1).stencilType),
                    static_cast<int>(sgpp::datadriven::StencilType::Block));
  BOOST_CHECK_EQUAL(config.stencils[1].applyOnLayers.size(), 1);
  BOOST_CHECK_EQUAL(config.stencils[1].applyOnLayers[0], 0);
  BOOST_CHECK_EQUAL(config.stencils[1].colorIndex, 1);
  BOOST_CHECK_EQUAL(config.stencils[1].blockLenght, 2);

  BOOST_CHECK_EQUAL(config.dim.size(), dim.size());
  for (size_t i = 0; i < config.dim.size(); i++) {
    BOOST_CHECK_EQUAL_COLLECTIONS(config.dim[i].begin(), config.dim[i].end(),
                                  dim[i].begin(), dim[i].end());
  }
}

BOOST_AUTO_TEST_CASE(testParallelConfig) {
  DataMiningConfigParser parser{datasetPath};

  ParallelConfiguration defaults;

  defaults.scalapackEnabled_ = false;
  defaults.processRows_ = 0;
  defaults.processCols_ = 0;
  defaults.rowBlockSize_ = 0;
  defaults.columnBlockSize_ = 0;
  ParallelConfiguration config;
  bool hasConfig;

  hasConfig = parser.getFitterParallelConfig(config, defaults);

  BOOST_CHECK_EQUAL(hasConfig, true);
  BOOST_ASSERT(config.scalapackEnabled_);
  BOOST_CHECK_EQUAL(config.processRows_, 4);
  BOOST_CHECK_EQUAL(config.processCols_, 1);
  BOOST_CHECK_EQUAL(config.rowBlockSize_, 64);
  BOOST_CHECK_EQUAL(config.columnBlockSize_, 128);
}

BOOST_AUTO_TEST_CASE(testVisualizationGeneralConfig) {
  DataMiningConfigParser parser{datasetPath};
  VisualizationGeneralConfig defaults;

  defaults.algorithm = "otherAlgorithm";
  defaults.targetDirectory = "./targetDirectory";
  defaults.targetFileType = VisualizationFileType::json;
  defaults.numBatches = 5;
  defaults.crossValidation = false;

  VisualizationGeneralConfig config;
  bool hasConfig;
  bool hasGeneralVisualizationConfig;

  hasGeneralVisualizationConfig = parser.hasVisualizationGeneralConfig();
  hasConfig = parser.getVisualizationGeneralConfig(config, defaults);

  BOOST_CHECK_EQUAL(hasConfig, true);
  BOOST_CHECK_EQUAL(hasGeneralVisualizationConfig, true);

  BOOST_CHECK_EQUAL(std::strcmp(config.algorithm.c_str(), "tsne"), 0);
  BOOST_CHECK_EQUAL(std::strcmp(config.targetDirectory.c_str(), "./output"), 0);
  BOOST_CHECK_EQUAL(static_cast<int>(config.targetFileType),
    static_cast<int>(VisualizationFileType::json));
  BOOST_CHECK_EQUAL(defaults.numBatches, 5);
  BOOST_CHECK_EQUAL(defaults.crossValidation, false);
}

BOOST_AUTO_TEST_CASE(testVisualizationParameters) {
  DataMiningConfigParser parser{datasetPath};
  VisualizationParameters defaults;

  defaults.perplexity = 22;
  defaults.theta = 0.1;
  defaults.targetDimension = 2;
  defaults.seed = 50;
  defaults.numberCores = 3;
  defaults.maxNumberIterations = 200;

  VisualizationParameters config;
  bool hasVisualizationParameters;
  bool hasConfig;

  hasVisualizationParameters = parser.hasVisualizationParametersConfig();
  hasConfig = parser.getVisualizationParameters(config, defaults);

  BOOST_CHECK_EQUAL(hasConfig, true);
  BOOST_CHECK_EQUAL(hasVisualizationParameters, true);
  BOOST_CHECK_EQUAL(config.perplexity, 30);
  BOOST_CHECK_EQUAL(config.theta, 0.5);
  BOOST_CHECK_EQUAL(config.targetDimension, 2);
  BOOST_CHECK_EQUAL(config.seed, 150);
  BOOST_CHECK_EQUAL(config.numberCores, 3);
  BOOST_CHECK_EQUAL(config.maxNumberIterations, 500);
}

BOOST_AUTO_TEST_SUITE_END()
