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

const auto datasetPath = "datadriven/tests/pipeline/config_configParser.json";
const auto multiDatasetPath = "datadriven/tests/pipeline/config_multiDatasetConfigParser.json";

BOOST_AUTO_TEST_SUITE(dataMiningConfigParserTest)

using sgpp::base::AdaptivityConfiguration;
using sgpp::base::AdaptivityThresholdType;
using sgpp::base::GridType;
using sgpp::base::RegularGridConfiguration;
using sgpp::datadriven::DataMiningConfigParser;
using sgpp::datadriven::DataSourceConfig;
using sgpp::datadriven::DataSourceFileType;
using sgpp::datadriven::DataSourceShufflingType;
using sgpp::datadriven::DataTransformationType;
using sgpp::datadriven::FitterType;
using sgpp::datadriven::ParallelConfiguration;
using sgpp::datadriven::RegularizationConfiguration;
using sgpp::datadriven::RegularizationMetricType;
using sgpp::datadriven::RegularizationType;
using sgpp::datadriven::ScorerConfiguration;
using sgpp::datadriven::ScorerMetricType;
using sgpp::datadriven::VisualizationGeneralConfig;
using sgpp::datadriven::VisualizationParameters;
using sgpp::solver::SLESolverConfiguration;
using sgpp::solver::SLESolverType;

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
  defaults.filePath_ = "something/false";
  defaults.fileType_ = DataSourceFileType::NONE;
  defaults.isCompressed_ = true;
  defaults.numBatches_ = 2;
  defaults.batchSize_ = 10;
  defaults.epochs_ = 1;
  defaults.shuffling_ = DataSourceShufflingType::sequential;
  defaults.validationPortion_ = 0.1;
  defaults.randomSeed_ = 1337;

  defaults.testFilePath_ = "something/testFalse";
  defaults.testFileType_ = DataSourceFileType::NONE;
  defaults.testNumBatches_ = 1;
  defaults.testBatchSize_ = 4;
  defaults.testIsCompressed_ = true;

  DataSourceConfig config;
  bool hasConfig;
  bool hasDataTransformationConfig;

  hasConfig = parser.getDataSourceConfig(config, defaults);
  hasDataTransformationConfig = parser.hasDataTransformationConfig();

  BOOST_CHECK_EQUAL(hasConfig, true);
  BOOST_CHECK_EQUAL(hasDataTransformationConfig, true);
  BOOST_CHECK_EQUAL(std::strcmp(config.filePath_.c_str(), "/path/to/some/file.arff"), 0);
  BOOST_CHECK_EQUAL(static_cast<int>(config.fileType_), static_cast<int>(DataSourceFileType::ARFF));
  BOOST_CHECK_EQUAL(config.numBatches_, 1);
  BOOST_CHECK_EQUAL(config.batchSize_, 0);
  BOOST_CHECK_EQUAL(static_cast<int>(config.dataTransformationConfig_.type_),
                    static_cast<int>(DataTransformationType::ROSENBLATT));
  BOOST_CHECK_EQUAL(config.dataTransformationConfig_.rosenblattConfig_.solverMaxIterations_, 1000);
  BOOST_CHECK_EQUAL(config.validationPortion_, 0.634);
  BOOST_CHECK_EQUAL(config.epochs_, 12);
  BOOST_CHECK_EQUAL(static_cast<int>(config.shuffling_),
                    static_cast<int>(DataSourceShufflingType::random));
  BOOST_CHECK_EQUAL(std::strcmp(config.testFilePath_.c_str(), "/path/to/some/testFile.arff"), 0);
  BOOST_CHECK_EQUAL(static_cast<int>(config.testFileType_),
                    static_cast<int>(DataSourceFileType::ARFF));
  BOOST_CHECK_EQUAL(config.testNumBatches_, 2);
  BOOST_CHECK_EQUAL(config.testBatchSize_, 16);
  BOOST_CHECK_EQUAL(config.testIsCompressed_, false);
}

BOOST_AUTO_TEST_CASE(testMultiDataSourceConfig) {
  DataMiningConfigParser parser{multiDatasetPath};

  std::vector<DataSourceConfig> defaults(2);
  defaults[0].filePath_ = "something/false";
  defaults[0].fileType_ = DataSourceFileType::NONE;
  defaults[0].isCompressed_ = true;
  defaults[0].numBatches_ = 2;
  defaults[0].batchSize_ = 10;
  defaults[0].epochs_ = 1;
  defaults[0].shuffling_ = DataSourceShufflingType::sequential;
  defaults[0].validationPortion_ = 0.1;
  defaults[0].randomSeed_ = 1337;

  std::vector<DataSourceConfig> config(2);
  bool hasConfig;
  bool hasDataTransformationConfig;

  hasConfig = parser.getMultiDataSourceConfig(config, defaults);
  hasDataTransformationConfig = parser.hasDataTransformationConfig();

  BOOST_CHECK_EQUAL(hasConfig, true);
  BOOST_CHECK_EQUAL(hasDataTransformationConfig, true);
  BOOST_CHECK_EQUAL(std::strcmp(config[0].filePath_.c_str(), "/path/to/some/file1.arff"), 0);
  BOOST_CHECK_EQUAL(static_cast<int>(config[0].fileType_),
                    static_cast<int>(DataSourceFileType::ARFF));
  BOOST_CHECK_EQUAL(config[0].numBatches_, 1);
  BOOST_CHECK_EQUAL(config[0].batchSize_, 0);
  BOOST_CHECK_EQUAL(static_cast<int>(config[0].dataTransformationConfig_.type_),
                    static_cast<int>(DataTransformationType::ROSENBLATT));
  BOOST_CHECK_EQUAL(config[0].dataTransformationConfig_.rosenblattConfig_.solverMaxIterations_,
                    1000);
  BOOST_CHECK_EQUAL(config[0].validationPortion_, 0.634);
  BOOST_CHECK_EQUAL(config[0].epochs_, 12);
  BOOST_CHECK_EQUAL(static_cast<int>(config[0].shuffling_),
                    static_cast<int>(DataSourceShufflingType::random));

  BOOST_CHECK_EQUAL(std::strcmp(config[1].filePath_.c_str(), "/path/to/some/file2.arff"), 0);
  BOOST_CHECK_EQUAL(static_cast<int>(config[1].fileType_),
                    static_cast<int>(DataSourceFileType::ARFF));
  BOOST_CHECK_EQUAL(config[1].numBatches_, 1);
  BOOST_CHECK_EQUAL(config[1].batchSize_, 0);
  BOOST_CHECK_EQUAL(static_cast<int>(config[1].dataTransformationConfig_.type_),
                    static_cast<int>(DataTransformationType::ROSENBLATT));
  BOOST_CHECK_EQUAL(config[1].dataTransformationConfig_.rosenblattConfig_.solverMaxIterations_,
                    1000);
  BOOST_CHECK_EQUAL(config[1].validationPortion_, 0.634);
  BOOST_CHECK_EQUAL(config[1].epochs_, 12);
  BOOST_CHECK_EQUAL(static_cast<int>(config[1].shuffling_),
                    static_cast<int>(DataSourceShufflingType::random));
}

BOOST_AUTO_TEST_CASE(testScorerConfig) {
  DataMiningConfigParser parser{datasetPath};

  ScorerConfiguration defaults;
  defaults.metric_ = ScorerMetricType::nll;
  ScorerConfiguration config;
  bool hasConfig;

  hasConfig = parser.getScorerConfig(config, defaults);

  BOOST_CHECK_EQUAL(hasConfig, true);
  BOOST_CHECK_EQUAL(static_cast<int>(config.metric_), static_cast<int>(ScorerMetricType::mse));
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
  defaults.thresholdType_ = AdaptivityThresholdType::Absolute;
  defaults.refinementThreshold_ = 42;
  defaults.coarseningThreshold_ = 42;
  defaults.maxLevelType_ = true;
  defaults.numRefinementPoints_ = 42;
  defaults.numCoarseningPoints_ = 42;
  defaults.coarsenInitialPoints_ = true;
  defaults.percent_ = 0.42;
  defaults.errorBasedRefinement_ = true;
  AdaptivityConfiguration config;
  bool hasConfig;
  double tolerance = 1E-5;

  hasConfig = parser.getFitterAdaptivityConfig(config, defaults);

  BOOST_CHECK_EQUAL(hasConfig, true);
  BOOST_CHECK_EQUAL(config.numRefinements_, 0);
  BOOST_CHECK_EQUAL(static_cast<int>(config.thresholdType_),
                    static_cast<int>(AdaptivityThresholdType::Relative));
  BOOST_CHECK_EQUAL(config.refinementThreshold_, 0);
  BOOST_CHECK_EQUAL(config.coarseningThreshold_, 1.0);
  BOOST_CHECK_EQUAL(config.maxLevelType_, false);
  BOOST_CHECK_EQUAL(config.numRefinementPoints_, 0);
  BOOST_CHECK_EQUAL(config.numCoarseningPoints_, 0);
  BOOST_CHECK_EQUAL(config.coarsenInitialPoints_, false);
  BOOST_CHECK_CLOSE(config.percent_, 0, tolerance);
  BOOST_CHECK_EQUAL(config.errorBasedRefinement_, false);
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
  defaults.optimizeLambda_ = false;
  defaults.optimizerTolerance_ = 1.0;
  defaults.convergenceThreshold_ = 1.0;
  defaults.intervalA_ = 1.0;
  defaults.intervalB_ = 0.0;
  RegularizationConfiguration config;
  bool hasConfig;
  double tolerance = 1E-5;

  hasConfig = parser.getFitterRegularizationConfig(config, defaults);

  BOOST_CHECK_EQUAL(hasConfig, true);
  BOOST_CHECK_EQUAL(static_cast<int>(config.type_), static_cast<int>(RegularizationType::Identity));
  BOOST_CHECK_CLOSE(config.lambda_, 1e-6, tolerance);
  BOOST_CHECK_CLOSE(config.exponentBase_, 3.0, tolerance);
  BOOST_CHECK_CLOSE(config.l1Ratio_, 4.0, tolerance);
  BOOST_CHECK_EQUAL(config.optimizeLambda_, true);
  BOOST_CHECK_CLOSE(config.optimizerTolerance_, 1e-10, tolerance);
  BOOST_CHECK_CLOSE(config.convergenceThreshold_, 1e-11, tolerance);
  BOOST_CHECK_CLOSE(config.intervalA_, 1e-3, tolerance);
  BOOST_CHECK_CLOSE(config.intervalB_, 0.5, tolerance);
  BOOST_CHECK_EQUAL(static_cast<int>(config.regularizationMetric_),
                    static_cast<int>(RegularizationMetricType::accuracy));
}

BOOST_AUTO_TEST_CASE(testFitterGeometryConfig) {
  DataMiningConfigParser parser{datasetPath};

  sgpp::datadriven::GeometryConfiguration defaults;
  defaults.dim_ = std::vector<std::vector<int64_t>>();
  sgpp::datadriven::GeometryConfiguration config;
  bool hasConfig;
  std::vector<std::vector<int64_t>> dim = {{1, 2}, {3, 4}};

  hasConfig = parser.getGeometryConfig(config, defaults);

  BOOST_CHECK_EQUAL(hasConfig, true);

  BOOST_CHECK_EQUAL(static_cast<int>(config.stencils_.at(0).stencilType_),
                    static_cast<int>(sgpp::datadriven::StencilType::DirectNeighbour));
  BOOST_CHECK_EQUAL(config.stencils_[0].applyOnLayers_.size(), 2);
  BOOST_CHECK_EQUAL(config.stencils_[0].applyOnLayers_[0], 0);
  BOOST_CHECK_EQUAL(config.stencils_[0].applyOnLayers_[1], 1);
  BOOST_CHECK_EQUAL(config.stencils_[0].colorIndex_, 0);

  BOOST_CHECK_EQUAL(static_cast<int>(config.stencils_.at(1).stencilType_),
                    static_cast<int>(sgpp::datadriven::StencilType::Block));
  BOOST_CHECK_EQUAL(config.stencils_[1].applyOnLayers_.size(), 1);
  BOOST_CHECK_EQUAL(config.stencils_[1].applyOnLayers_[0], 0);
  BOOST_CHECK_EQUAL(config.stencils_[1].colorIndex_, 1);
  BOOST_CHECK_EQUAL(config.stencils_[1].blockLenght_, 2);

  BOOST_CHECK_EQUAL(config.dim_.size(), dim.size());
  for (size_t i = 0; i < config.dim_.size(); i++) {
    BOOST_CHECK_EQUAL_COLLECTIONS(config.dim_[i].begin(), config.dim_[i].end(), dim[i].begin(),
                                  dim[i].end());
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

  defaults.algorithm_ = std::vector<std::string>({"otherAlgorithm"});
  defaults.targetDirectory_ = "./targetDirectory";
  defaults.targetFileType_ = VisualizationFileType::json;
  defaults.numBatches_ = 5;
  defaults.crossValidation_ = false;

  VisualizationGeneralConfig config;
  bool hasConfig;
  bool hasGeneralVisualizationConfig;
  std::vector<std::string> expectedAlgorithm = std::vector<std::string>({"tsne", "heatmaps"});

  hasGeneralVisualizationConfig = parser.hasVisualizationGeneralConfig();

  hasConfig = parser.getVisualizationGeneralConfig(config, defaults);

  BOOST_CHECK_EQUAL(hasConfig, true);
  BOOST_CHECK_EQUAL(hasGeneralVisualizationConfig, true);

  BOOST_CHECK_EQUAL_COLLECTIONS(config.algorithm_.begin(), config.algorithm_.end(),
                                expectedAlgorithm.begin(), expectedAlgorithm.end());
  BOOST_CHECK_EQUAL(std::strcmp(config.targetDirectory_.c_str(), "./output"), 0);
  BOOST_CHECK_EQUAL(static_cast<int>(config.targetFileType_),
                    static_cast<int>(VisualizationFileType::json));
  BOOST_CHECK_EQUAL(defaults.numBatches_, 5);
  BOOST_CHECK_EQUAL(defaults.crossValidation_, false);
}

BOOST_AUTO_TEST_CASE(testVisualizationParameters) {
  DataMiningConfigParser parser{datasetPath};
  VisualizationParameters defaults;

  defaults.perplexity_ = 22;
  defaults.theta_ = 0.1;
  defaults.targetDimension_ = 2;
  defaults.seed_ = 50;
  defaults.numberCores_ = 3;
  defaults.maxNumberIterations_ = 200;

  VisualizationParameters config;
  bool hasVisualizationParameters;
  bool hasConfig;

  hasVisualizationParameters = parser.hasVisualizationParametersConfig();
  hasConfig = parser.getVisualizationParameters(config, defaults);

  BOOST_CHECK_EQUAL(hasConfig, true);
  BOOST_CHECK_EQUAL(hasVisualizationParameters, true);
  BOOST_CHECK_EQUAL(config.perplexity_, 30);
  BOOST_CHECK_EQUAL(config.theta_, 0.5);
  BOOST_CHECK_EQUAL(config.targetDimension_, 2);
  BOOST_CHECK_EQUAL(config.seed_, 150);
  BOOST_CHECK_EQUAL(config.numberCores_, 3);
  BOOST_CHECK_EQUAL(config.maxNumberIterations_, 500);
}

BOOST_AUTO_TEST_SUITE_END()
