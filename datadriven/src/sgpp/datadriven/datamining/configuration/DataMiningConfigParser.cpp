// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/datadriven/datamining/configuration/DataMiningConfigParser.hpp>

#include <sgpp/base/exception/data_exception.hpp>
#include <sgpp/base/exception/file_exception.hpp>

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/grid/LevelIndexTypes.hpp>
#include <sgpp/base/grid/GeneralGridTypeParser.hpp>
#include <sgpp/base/grid/GridTypeParser.hpp>
#include <sgpp/base/grid/AdaptivityThresholdTypeParser.hpp>
#include <sgpp/base/grid/CoarseningFunctorTypeParser.hpp>
#include <sgpp/base/grid/RefinementFunctorTypeParser.hpp>

#include <sgpp/base/tools/json/JSON.hpp>
#include <sgpp/base/tools/json/ListNode.hpp>
#include <sgpp/base/tools/json/json_exception.hpp>

#include <sgpp/datadriven/configuration/GeometryConfiguration.hpp>
#include <sgpp/datadriven/configuration/RegularizationConfiguration.hpp>
#include <sgpp/datadriven/configuration/DensityEstimationTypeParser.hpp>
#include <sgpp/datadriven/configuration/GeometryConfigurationParser.hpp>
#include <sgpp/datadriven/configuration/MatrixDecompositionTypeParser.hpp>
#include <sgpp/datadriven/configuration/RegularizationTypeParser.hpp>

#include <sgpp/datadriven/datamining/modules/dataSource/DataSourceFileTypeParser.hpp>
#include <sgpp/datadriven/datamining/modules/dataSource/DataTransformationTypeParser.hpp>
#include <sgpp/datadriven/datamining/modules/dataSource/shuffling/DataSourceShufflingTypeParser.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/FitterTypeParser.hpp>
#include <sgpp/datadriven/datamining/modules/scoring/ScorerMetricTypeParser.hpp>
#include <sgpp/datadriven/datamining/modules/visualization/VisualizationTypesParser.hpp>

#include <sgpp/solver/TypesSolver.hpp>
#include <sgpp/solver/SLESolverTypeParser.hpp>

#include <string>
#include <vector>
#include <iostream>
#include <map>

using json::DictNode;
using json::JSON;
using json::json_exception;
using sgpp::base::AdaptivityConfiguration;
using sgpp::base::data_exception;
using sgpp::base::file_exception;
using sgpp::base::GeneralGridConfiguration;
using sgpp::solver::SLESolverConfiguration;

namespace sgpp {
namespace datadriven {

const std::string DataMiningConfigParser::dataSource = "dataSource";
const std::string DataMiningConfigParser::scorer = "scorer";
const std::string DataMiningConfigParser::fitter = "fitter";
const std::string DataMiningConfigParser::visualization = "visualization";

DataMiningConfigParser::DataMiningConfigParser(const std::string &filepath) : configFile(nullptr) {
  try {
    configFile = std::make_unique<JSON>(filepath);
  } catch (json_exception &exception) {
    std::cout << exception.what() << std::endl;
    throw file_exception("Cannot open JSON file.");
  }
}

DataMiningConfigParser::~DataMiningConfigParser() {}

bool DataMiningConfigParser::hasDataSourceConfig() const {
  return configFile->contains(dataSource);
}

bool DataMiningConfigParser::hasDataTransformationConfig() const {
  bool hasDataTransformationConfig =
      hasDataSourceConfig() ? (*configFile)[dataSource].contains("dataTransformation") : false;
  return hasDataTransformationConfig;
}

bool DataMiningConfigParser::hasScorerConfig() const { return configFile->contains(scorer); }

bool DataMiningConfigParser::hasFitterConfig() const { return configFile->contains(fitter); }

bool DataMiningConfigParser::hasParallelConfig() const {
  bool hasParallelConfig =
      hasFitterConfig() ? (*configFile)[fitter].contains("parallelConfig") : false;
  return hasParallelConfig;
}

bool DataMiningConfigParser::hasFitterConfigCrossValidation() const {
  bool hasFitterCrossValidationConfig =
      hasFitterConfig() ? (*configFile)[fitter].contains("crossValidation") : false;
  return hasFitterCrossValidationConfig;
}

bool DataMiningConfigParser::hasGeometryConfig() const {
  return configFile->contains("geometryConfig");
}

bool DataMiningConfigParser::hasVisualizationConfig() const {
  return configFile->contains(visualization);
}

bool DataMiningConfigParser::hasVisualizationGeneralConfig() const {
  bool hasVisualizationParameters =
      hasVisualizationConfig() ? (*configFile)[visualization].contains("generalConfig") : false;
  return hasVisualizationParameters;
}

bool DataMiningConfigParser::hasVisualizationParametersConfig() const {
  bool hasVisualizationParameters =
      hasVisualizationConfig() ? (*configFile)[visualization].contains("parameters") : false;
  return hasVisualizationParameters;
}

bool DataMiningConfigParser::getDataSourceConfig(DataSourceConfig &config,
                                                 const DataSourceConfig &defaults) const {
  bool hasDataSource = hasDataSourceConfig();

  if (hasDataSource) {
    auto dataSourceConfig = static_cast<DictNode *>(&(*configFile)[dataSource]);

    config.filePath_ = parseString(*dataSourceConfig, "filePath", defaults.filePath_, "dataSource");
    config.isCompressed_ =
        parseBool(*dataSourceConfig, "compression", defaults.isCompressed_, "dataSource");
    config.numBatches_ =
        parseUInt(*dataSourceConfig, "numBatches", defaults.numBatches_, "dataSource");
    config.batchSize_ =
        parseUInt(*dataSourceConfig, "batchSize", defaults.batchSize_, "dataSource");
    config.hasTargets_ =
        parseBool(*dataSourceConfig, "hasTargets", defaults.hasTargets_, "dataSource");
    config.validationPortion_ = parseDouble(*dataSourceConfig, "validationPortion",
                                            defaults.validationPortion_, "dataSource");
    // if negative we want UINT_MAX here, so all should be fine
    config.readinCutoff_ = static_cast<size_t>(
        parseInt(*dataSourceConfig, "readinCutoff", defaults.readinCutoff_, "dataSource"));
    config.readinClasses_ =
        parseDoubleArray(*dataSourceConfig, "readinClasses", defaults.readinClasses_, "dataSource");
    config.readinColumns_ =
        parseUIntArray(*dataSourceConfig, "readinColumns", defaults.readinColumns_, "dataSource");

    // parse file type
    if (dataSourceConfig->contains("fileType")) {
      config.fileType_ = DataSourceFileTypeParser::parse((*dataSourceConfig)["fileType"].get());
    } else {
      std::cout << "# Did not find " << dataSource << "[fileType]. Setting default value "
                << DataSourceFileTypeParser::toString(defaults.fileType_) << "." << std::endl;
      config.fileType_ = defaults.fileType_;
    }

    // parse dataTransformationConfig
    bool hasDataTransformation = hasDataTransformationConfig();

    if (hasDataTransformation) {
      auto dataTransformationConfig =
          static_cast<DictNode *>(&(*configFile)[dataSource]["dataTransformation"]);
      parseDataTransformationConfig(*dataTransformationConfig, config.dataTransformationConfig_,
                                    defaults.dataTransformationConfig_, "dataTransformation");
    } else {
      std::cout << "# Could not find specification of dataSource[dataTransformationConfig]. "
                   "Falling back to default values."
                << std::endl;
      config.dataTransformationConfig_ = defaults.dataTransformationConfig_;
    }

    // parse the shuffling
    if (dataSourceConfig->contains("shuffling")) {
      config.shuffling_ =
          DataSourceShufflingTypeParser::parse((*dataSourceConfig)["shuffling"].get());
    } else {
      std::cout << "# Did not find dataSource[shuffling]. Setting default value "
                << DataSourceShufflingTypeParser::toString(defaults.shuffling_) << "." << std::endl;
      config.shuffling_ = defaults.shuffling_;
    }

    config.randomSeed_ =
        parseUInt(*dataSourceConfig, "randomSeed", defaults.randomSeed_, "dataSource");
    config.epochs_ = parseUInt(*dataSourceConfig, "epochs", defaults.epochs_, "dataSource");

    // Parse info for test data
    config.testFilePath_ =
        parseString(*dataSourceConfig, "testFilePath", defaults.filePath_, "dataSource");

    // parse file type of test data
    if (dataSourceConfig->contains("testFileType")) {
      config.testFileType_ =
          DataSourceFileTypeParser::parse((*dataSourceConfig)["testFileType"].get());
    } else {
      std::cout << "# Did not find " << dataSource << "[testFileType]. Setting default value "
                << DataSourceFileTypeParser::toString(defaults.testFileType_) << "." << std::endl;

      config.testFileType_ = defaults.testFileType_;
    }

    config.testIsCompressed_ =
        parseBool(*dataSourceConfig, "testCompression", defaults.testIsCompressed_, "dataSource");
    config.testNumBatches_ =
        parseUInt(*dataSourceConfig, "testNumBatches", defaults.testNumBatches_, "dataSource");
    config.testBatchSize_ =
        parseUInt(*dataSourceConfig, "testBatchSize", defaults.testBatchSize_, "dataSource");
    config.testHasTargets_ =
        parseBool(*dataSourceConfig, "testHasTargets", defaults.testHasTargets_, "dataSource");

    config.testReadinCutoff_ = static_cast<size_t>(
        parseInt(*dataSourceConfig, "testReadinCutoff", defaults.testReadinCutoff_, "dataSource"));
    config.testReadinClasses_ = parseDoubleArray(*dataSourceConfig, "testReadinClasses",
                                                 defaults.testReadinClasses_, "dataSource");
    config.testReadinColumns_ = parseUIntArray(*dataSourceConfig, "testReadinColumns",
                                               defaults.testReadinColumns_, "dataSource");

  } else {
    std::cout << "# Could not find specification of dataSource. Falling Back to default values."
              << std::endl;
    config = defaults;
  }
  return hasDataSource;
}

bool DataMiningConfigParser::getMultiDataSourceConfig(
    std::vector<DataSourceConfig> &config, const std::vector<DataSourceConfig> &defaults) const {
  bool hasDataSource = hasDataSourceConfig();

  if (hasDataSource) {
    auto dataSourceConfig = dynamic_cast<json::DictNode *>(&(*configFile)[dataSource]);

    // Fill in all parameters for first dataset (except the filePath)
    config[0].isCompressed_ =
        parseBool(*dataSourceConfig, "compression", defaults[0].isCompressed_, "dataSource");
    config[0].numBatches_ =
        parseUInt(*dataSourceConfig, "numBatches", defaults[0].numBatches_, "dataSource");
    config[0].batchSize_ =
        parseUInt(*dataSourceConfig, "batchSize", defaults[0].batchSize_, "dataSource");
    config[0].hasTargets_ =
        parseBool(*dataSourceConfig, "hasTargets", defaults[0].hasTargets_, "dataSource");
    config[0].validationPortion_ = parseDouble(*dataSourceConfig, "validationPortion",
                                               defaults[0].validationPortion_, "dataSource");
    // if negative we want UINT_MAX here, so all should be fine
    config[0].readinCutoff_ = static_cast<size_t>(
        parseInt(*dataSourceConfig, "readinCutoff", defaults[0].readinCutoff_, "dataSource"));
    config[0].readinClasses_ = parseDoubleArray(*dataSourceConfig, "readinClasses",
                                                defaults[0].readinClasses_, "dataSource");
    config[0].readinColumns_ = parseUIntArray(*dataSourceConfig, "readinColumns",
                                              defaults[0].readinColumns_, "dataSource");

    // parse file type
    if (dataSourceConfig->contains("fileType")) {
      config[0].fileType_ = DataSourceFileTypeParser::parse((*dataSourceConfig)["fileType"].get());
    } else {
      std::cout << "# Did not find " << dataSource << "[fileType]. Setting default value "
                << DataSourceFileTypeParser::toString(defaults[0].fileType_) << "." << std::endl;
      config[0].fileType_ = defaults[0].fileType_;
    }

    // parse dataTransformationConfig
    bool hasDataTransformation = hasDataTransformationConfig();

    if (hasDataTransformation) {
      auto dataTransformationConfig =
          static_cast<DictNode *>(&(*configFile)[dataSource]["dataTransformation"]);
      parseDataTransformationConfig(*dataTransformationConfig, config[0].dataTransformationConfig_,
                                    defaults[0].dataTransformationConfig_, "dataTransformation");
    } else {
      std::cout << "# Could not find specification of dataSource[dataTransformationConfig]. "
                   "Falling back to default values."
                << std::endl;
      config[0].dataTransformationConfig_ = defaults[0].dataTransformationConfig_;
    }

    // parse the shuffling
    if (dataSourceConfig->contains("shuffling")) {
      config[0].shuffling_ =
          DataSourceShufflingTypeParser::parse((*dataSourceConfig)["shuffling"].get());
    } else {
      std::cout << "# Did not find dataSource[shuffling]. Setting default value "
                << DataSourceShufflingTypeParser::toString(defaults[0].shuffling_) << "."
                << std::endl;
      config[0].shuffling_ = defaults[0].shuffling_;
    }

    config[0].randomSeed_ =
        parseUInt(*dataSourceConfig, "randomSeed", defaults[0].randomSeed_, "dataSource");
    config[0].epochs_ = parseUInt(*dataSourceConfig, "epochs", defaults[0].epochs_, "dataSource");

    // Parse the filePath field for all config instances
    json::ListNode &paths = dynamic_cast<json::ListNode &>((*dataSourceConfig)["filePath"]);

    // First copy all parameters to all other config instances...
    for (size_t i = 1; i < paths.size(); i++) {
      config[i] = config[0];
    }
    // ... then fill in the correct paths
    for (size_t i = 0; i < paths.size(); i++) {
      config[i].filePath_ = paths[i].get();
    }
  } else {
    std::cout << "# Could not find specification of dataSource. Falling Back to default values."
              << std::endl;
    config = defaults;
  }
  return hasDataSource;
}

bool DataMiningConfigParser::getScorerConfig(ScorerConfiguration &config,
                                             const ScorerConfiguration &defaults) const {
  bool hasScorer = hasScorerConfig();

  if (hasScorer) {
    auto scorerConfig = static_cast<DictNode *>(&(*configFile)[scorer]);

    // parse metric type
    if (scorerConfig->contains("metric")) {
      config.metric_ = ScorerMetricTypeParser::parse((*scorerConfig)["metric"].get());
    } else {
      std::cout << "# Did not find scorer[metric]. Setting default value "
                << ScorerMetricTypeParser::toString(defaults.metric_) << "." << std::endl;
    }
  } else {
    std::cout << "# Could not find specification of scorer. Falling Back to default values."
              << std::endl;
    config = defaults;
  }
  return hasScorer;
}

// TODO(lettrich): is this consistent with the rest of the parsing?
bool DataMiningConfigParser::getFitterConfigType(FitterType &config,
                                                 const FitterType &defaults) const {
  bool hasFitterConfig = this->hasFitterConfig();

  if (hasFitterConfig) {
    auto fitterConfig = static_cast<DictNode *>(&(*configFile)[fitter]);
    if (fitterConfig->contains("type")) {
      config = FitterTypeParser::parse((*fitterConfig)["type"].get());
    } else {
      std::cout << "# Could not find specification of fitter[type]. Falling Back to default values."
                << std::endl;
      config = defaults;
    }
  }

  return hasFitterConfig;
}

bool DataMiningConfigParser::getFitterGridConfig(GeneralGridConfiguration &config,
                                                 const GeneralGridConfiguration &defaults) const {
  bool hasFitterGridConfig =
      hasFitterConfig() ? (*configFile)[fitter].contains("gridConfig") : false;

  if (hasFitterGridConfig) {
    auto fitterConfig = static_cast<DictNode *>(&(*configFile)[fitter]["gridConfig"]);
    config.dim_ = parseUInt(*fitterConfig, "dim", defaults.dim_, "gridConfig");
    config.level_ =
        static_cast<int>(parseInt(*fitterConfig, "level", defaults.level_, "gridConfig"));
    config.levelVector_ = static_cast<std::vector<size_t>>(
        parseUIntArray(*fitterConfig, "levelVector", defaults.levelVector_, "gridConfig"));
    config.maxDegree_ = parseUInt(*fitterConfig, "maxDegree", defaults.maxDegree_, "gridConfig");
    config.boundaryLevel_ = static_cast<unsigned int>(
        parseUInt(*fitterConfig, "boundaryLevel", defaults.boundaryLevel_, "gridConfig"));
    config.filename_ = parseString(*fitterConfig, "fileName", defaults.filename_, "gridConfig");

    // parse general grid type
    if (fitterConfig->contains("generalGridType")) {
      if ((*fitterConfig)["generalGridType"].size() == 1) {
        config.generalType_ =
            base::GeneralGridTypeParser::parse((*fitterConfig)["generalGridType"].get());
      } else {
        config.generalType_ =
            base::GeneralGridTypeParser::parse((*fitterConfig)["generalGridType"]["value"].get());
      }
    } else {
      std::cout << "# Did not find gridConfig[generalGridType]. Setting default value."
                << std::endl;
    }

    // parse  grid type
    if (fitterConfig->contains("gridType")) {
      if ((*fitterConfig)["gridType"].size() == 1) {
        config.type_ = base::GridTypeParser::parse((*fitterConfig)["gridType"].get());
      } else {
        config.type_ = base::GridTypeParser::parse((*fitterConfig)["gridType"]["value"].get());
      }
    } else {
      std::cout << "# Did not find gridConfig[gridType]. Setting default value "
                << base::GridTypeParser::toString(defaults.type_) << "." << std::endl;
      config.type_ = defaults.type_;
    }
  } else {
    std::cout
        << "# Could not find specification of fitter[gridConfig]. Falling Back to default values."
        << std::endl;
    config = defaults;
  }
  return hasFitterGridConfig;
}

bool DataMiningConfigParser::getFitterAdaptivityConfig(
    sgpp::base::AdaptivityConfiguration &config,
    const sgpp::base::AdaptivityConfiguration &defaults) const {
  bool hasFitterAdaptivityConfig =
      hasFitterConfig() ? (*configFile)[fitter].contains("adaptivityConfig") : false;

  if (hasFitterAdaptivityConfig) {
    auto adaptivityConfig = static_cast<DictNode *>(&(*configFile)[fitter]["adaptivityConfig"]);
    config.numRefinements_ = parseUInt(*adaptivityConfig, "numRefinements",
                                       defaults.numRefinements_, "adaptivityConfig");
    config.refinementThreshold_ = parseDouble(*adaptivityConfig, "refinementThreshold",
                                              defaults.refinementThreshold_, "adaptivityConfig");
    config.coarseningThreshold_ = parseDouble(*adaptivityConfig, "coarseningThreshold",
                                              defaults.coarseningThreshold_, "adaptivityConfig");
    config.maxLevelType_ =
        parseBool(*adaptivityConfig, "maxLevelType", defaults.maxLevelType_, "adaptivityConfig");
    config.numRefinementPoints_ = parseUInt(*adaptivityConfig, "numRefinementPoints",
                                            defaults.numRefinementPoints_, "adaptivityConfig");
    config.numCoarseningPoints_ = parseUInt(*adaptivityConfig, "numCoarseningPoints",
                                            defaults.numRefinementPoints_, "adaptivityConfig");
    config.coarsenInitialPoints_ = parseBool(*adaptivityConfig, "coarsenInitialPoints",
                                             defaults.coarsenInitialPoints_, "adaptivityConfig");
    config.percent_ =
        parseDouble(*adaptivityConfig, "percent", defaults.percent_, "adaptivityConfig");
    config.errorBasedRefinement_ = parseBool(*adaptivityConfig, "errorBasedRefinement",
                                             defaults.errorBasedRefinement_, "adaptivityConfig");
    config.errorConvergenceThreshold_ =
        parseDouble(*adaptivityConfig, "errorConvergenceThreshold",
                    defaults.errorConvergenceThreshold_, "adaptivityConfig");
    config.errorBufferSize_ = parseUInt(*adaptivityConfig, "errorBufferSize",
                                        defaults.errorBufferSize_, "adaptivityConfig");
    config.errorMinInterval_ = parseUInt(*adaptivityConfig, "errorMinInterval",
                                         defaults.errorMinInterval_, "adaptivityConfig");
    config.refinementPeriod_ = parseUInt(*adaptivityConfig, "refinementPeriod",
                                         defaults.refinementPeriod_, "adaptivityConfig");
    config.precomputeEvaluations_ = parseBool(*adaptivityConfig, "precomputeEvaluations",
                                              defaults.precomputeEvaluations_, "adaptivityConfig");
    config.levelPenalize_ =
        parseBool(*adaptivityConfig, "penalizeLevels", defaults.levelPenalize_, "adaptivityConfig");

    // Parse scaling coefficients if present
    if (adaptivityConfig->contains("scalingCoefficients")) {
      json::ListNode &coefs =
          dynamic_cast<json::ListNode &>((*adaptivityConfig)["scalingCoefficients"]);
      for (size_t i = 0; i < coefs.size(); i++) {
        config.scalingCoefficients_.push_back(coefs[i].getDouble());
      }
    }

    // Parse refinement indicator
    if (adaptivityConfig->contains("refinementIndicator")) {
      config.refinementFunctorType_ = base::RefinementFunctorTypeParser::parse(
          (*adaptivityConfig)["refinementIndicator"].get());
    } else {
      std::cout << "# Did not find adaptivityConfig[refinementIndicator]. Setting default value "
                << base::RefinementFunctorTypeParser::toString(defaults.refinementFunctorType_)
                << "." << std::endl;
      config.refinementFunctorType_ = defaults.refinementFunctorType_;
    }

    // Parse coarsening indicator
    if (adaptivityConfig->contains("coarseningIndicator")) {
      config.coarseningFunctorType_ = base::CoarseningFunctorTypeParser::parse(
          (*adaptivityConfig)["coarseningIndicator"].get());
    } else {
      std::cout << "# Did not find adaptivityConfig[coarseningIndicator]. Setting default value "
                << base::CoarseningFunctorTypeParser::toString(defaults.coarseningFunctorType_)
                << "." << std::endl;
      config.coarseningFunctorType_ = defaults.coarseningFunctorType_;
    }

    // Parse threshold type
    if (adaptivityConfig->contains("thresholdType")) {
      config.thresholdType_ =
          base::AdaptivityThresholdTypeParser::parse((*adaptivityConfig)["thresholdType"].get());
    } else {
      std::cout << "# Did not find adaptivityConfig[thresholdType]. Setting default value "
                << base::AdaptivityThresholdTypeParser::toString(defaults.thresholdType_) << "."
                << std::endl;
      config.thresholdType_ = defaults.thresholdType_;
    }

  } else {
    std::cout << "# Could not find specification of fitter[adaptivityConfig]. Falling Back to "
                 "default values."
              << std::endl;
    config = defaults;
  }
  return hasFitterAdaptivityConfig;
}

bool DataMiningConfigParser::getFitterCrossvalidationConfig(
    CrossvalidationConfiguration &config, const CrossvalidationConfiguration &defaults) const {
  bool hasFitterCrossvalidationConfig =
      hasFitterConfig() ? (*configFile)[fitter].contains("crossValidation") : false;

  if (hasFitterCrossvalidationConfig) {
    auto crossvalidationConfig = static_cast<DictNode *>(&(*configFile)[fitter]["crossValidation"]);
    config.enable_ =
        parseBool(*crossvalidationConfig, "enable", defaults.enable_, "crossValidation");
    config.kfold_ = parseUInt(*crossvalidationConfig, "kFold", defaults.kfold_, "crossValidation");
    config.seed_ = static_cast<int>(
        parseInt(*crossvalidationConfig, "randomSeed", defaults.seed_, "crossvalidationConfig"));
    config.shuffle_ =
        parseBool(*crossvalidationConfig, "shuffle", defaults.shuffle_, "crossValidation");
    config.silent_ =
        parseBool(*crossvalidationConfig, "silent", defaults.silent_, "crossValidation");
    config.lambda_ =
        parseDouble(*crossvalidationConfig, "lambda", defaults.lambda_, "crossValidation");
    config.lambdaStart_ = parseDouble(*crossvalidationConfig, "lambdaStart", defaults.lambdaStart_,
                                      "crossValidation");
    config.lambdaEnd_ =
        parseDouble(*crossvalidationConfig, "lambdaEnd", defaults.lambdaEnd_, "crossValidation");
    config.lambdaSteps_ =
        parseUInt(*crossvalidationConfig, "lambdaSteps", defaults.lambdaSteps_, "crossValidation");
    config.logScale_ =
        parseBool(*crossvalidationConfig, "logScale", defaults.logScale_, "crossValidation");
  } else {
    std::cout << "# Could not find specification of fitter[crossvalidationConfig]. Falling Back to "
                 "default values."
              << std::endl;
    config = defaults;
  }
  return hasFitterCrossvalidationConfig;
}

bool DataMiningConfigParser::getFitterDensityEstimationConfig(
    DensityEstimationConfiguration &config, const DensityEstimationConfiguration &defaults) const {
  bool hasFitterDensityEstimationConfig =
      hasFitterConfig() ? (*configFile)[fitter].contains("densityEstimationConfig") : false;

  if (hasFitterDensityEstimationConfig) {
    auto densityEstimationConfig =
        static_cast<DictNode *>(&(*configFile)[fitter]["densityEstimationConfig"]);
    config.iCholSweepsDecompose_ =
        parseUInt(*densityEstimationConfig, "iCholSweepsDecompose", defaults.iCholSweepsDecompose_,
                  "densityEstimationConfig");
    config.iCholSweepsRefine_ = parseUInt(*densityEstimationConfig, "iCholSweepsRefine",
                                          defaults.iCholSweepsRefine_, "densityEstimationConfig");
    config.iCholSweepsUpdateLambda_ =
        parseUInt(*densityEstimationConfig, "iCholSweepsUpdateLambda",
                  defaults.iCholSweepsUpdateLambda_, "densityEstimationConfig");
    config.iCholSweepsSolver_ = parseUInt(*densityEstimationConfig, "iCholSweepsSolver",
                                          defaults.iCholSweepsSolver_, "densityEstimationConfig");

    config.normalize_ = parseBool(*densityEstimationConfig, "normalize", defaults.normalize_,
                                  "densityEstimationConfig");

    config.useOfflinePermutation_ =
        parseBool(*densityEstimationConfig, "useOfflinePermutation",
                  defaults.useOfflinePermutation_, "densityEstimationConfig");

    // parse  density estimation type
    if (densityEstimationConfig->contains("densityEstimationType")) {
      config.type_ = DensityEstimationTypeParser::parse(
          (*densityEstimationConfig)["densityEstimationType"].get());
    } else {
      std::cout
          << "# Did not find densityEstimationConfig[densityEstimationType]. Setting default value "
          << DensityEstimationTypeParser::toString(defaults.type_) << "." << std::endl;
      config.type_ = defaults.type_;
    }

    // parse matrix decomposition type
    if (densityEstimationConfig->contains("matrixDecompositionType")) {
      config.decomposition_ = MatrixDecompositionTypeParser::parse(
          (*densityEstimationConfig)["matrixDecompositionType"].get());
    } else {
      std::cout << "# Did not find densityEstimationConfig[matrixDecompositionType]. Setting "
                   "default value "
                << MatrixDecompositionTypeParser::toString(defaults.decomposition_) << "."
                << std::endl;
      config.decomposition_ = defaults.decomposition_;
    }
  } else {
    std::cout << "# Could not find specification of fitter[densityEstimationConfig]. Falling Back "
                 "to default values."
              << std::endl;
    config = defaults;
  }
  return hasFitterDensityEstimationConfig;
}

bool DataMiningConfigParser::getFitterSolverRefineConfig(
    SLESolverConfiguration &config, const SLESolverConfiguration &defaults) const {
  bool hasFitterSolverRefineConfig =
      hasFitterConfig() ? (*configFile)[fitter].contains("solverRefineConfig") : false;

  if (hasFitterSolverRefineConfig) {
    auto solverConfig = static_cast<DictNode *>(&(*configFile)[fitter]["solverRefineConfig"]);

    parseSLESolverConfig(*solverConfig, config, defaults, "solverRefineConfig");
  } else {
    std::cout << "# Could not find specification of fitter[solverRefineConfig]. Falling Back to "
                 "default values."
              << std::endl;
    config = defaults;
  }
  return hasFitterSolverRefineConfig;
}

bool DataMiningConfigParser::getFitterSolverFinalConfig(
    SLESolverConfiguration &config, const SLESolverConfiguration &defaults) const {
  bool hasFitterSolverFinalConfig =
      hasFitterConfig() ? (*configFile)[fitter].contains("solverFinalConfig") : false;

  if (hasFitterSolverFinalConfig) {
    auto solverConfig = static_cast<DictNode *>(&(*configFile)[fitter]["solverFinalConfig"]);

    parseSLESolverConfig(*solverConfig, config, defaults, "solverFinalConfig");
  } else {
    std::cout << "# Could not find specification of fitter[solverFinalConfig]. Falling Back to "
                 "default values."
              << std::endl;
    config = defaults;
  }
  return hasFitterSolverFinalConfig;
}

bool DataMiningConfigParser::getFitterRegularizationConfig(
    RegularizationConfiguration &config, const RegularizationConfiguration &defaults) const {
  bool hasRegularizationConfig =
      hasFitterConfig() ? (*configFile)[fitter].contains("regularizationConfig") : false;

  if (hasRegularizationConfig) {
    auto regularizationConfig =
        static_cast<DictNode *>(&(*configFile)[fitter]["regularizationConfig"]);

    // parse  regularization type
    if (regularizationConfig->contains("regularizationType")) {
      config.type_ =
          RegularizationTypeParser::parse((*regularizationConfig)["regularizationType"].get());
    } else {
      std::cout << "# Did not find regularizationConfig[regularizationType]. Setting default value "
                << RegularizationTypeParser::toString(defaults.type_) << "." << std::endl;
      config.type_ = defaults.type_;
    }

    config.lambda_ =
        parseDouble(*regularizationConfig, "lambda", defaults.lambda_, "regularizationConfig");

    config.exponentBase_ = parseDouble(*regularizationConfig, "exponentBase",
                                       defaults.exponentBase_, "regularizationConfig");

    config.l1Ratio_ =
        parseDouble(*regularizationConfig, "l1Ratio", defaults.l1Ratio_, "regularizationConfig");

    config.optimizeLambda_ = parseBool(*regularizationConfig, "optimizeLambda",
                                       defaults.optimizeLambda_, "regularizationConfig");

    config.optimizerTolerance_ = parseDouble(*regularizationConfig, "optimizerTolerance",
                                             defaults.optimizerTolerance_, "regularizationConfig");

    config.convergenceThreshold_ =
        parseDouble(*regularizationConfig, "convergenceThreshold", defaults.convergenceThreshold_,
                    "regularizationConfig");

    config.intervalA_ = parseDouble(*regularizationConfig, "intervalA", defaults.intervalA_,
                                    "regularizationConfig");
    config.intervalB_ = parseDouble(*regularizationConfig, "intervalB", defaults.intervalB_,
                                    "regularizationConfig");

    if (regularizationConfig->contains("regularizationMetric")) {
      config.regularizationMetric_ = ScorerMetricTypeParser::parseRegularizationMetric(
          (*regularizationConfig)["regularizationMetric"].get());
    } else {
      std::cout << "# Did not find scorer[metric]. Setting default value "
                << ScorerMetricTypeParser::regularizationMetricToString(
                       defaults.regularizationMetric_)
                << "." << std::endl;
    }
  }

  return hasRegularizationConfig;
}

bool DataMiningConfigParser::getVisualizationGeneralConfig(
    VisualizationGeneralConfig &config, const VisualizationGeneralConfig &defaults) const {
  bool hasVisualization = hasVisualizationConfig();

  if (!hasVisualization) {
    config.execute_ = false;
  } else {
    config.execute_ = true;
  }
  bool hasGeneralConfig = hasVisualizationGeneralConfig();

  if (hasGeneralConfig) {
    auto visualizationGeneralConfig =
        static_cast<DictNode *>(&(*configFile)[visualization]["generalConfig"]);

    std::cout << "Starting reading visualization " << std::endl;
    config.algorithm_ = parseStringArray(*visualizationGeneralConfig, "algorithm",
                                         defaults.algorithm_, "visualization");

    config.targetDirectory_ = parseString(*visualizationGeneralConfig, "targetDirectory",
                                          defaults.targetDirectory_, "visualization");

    // parse file type
    if (visualizationGeneralConfig->contains("targetFileType")) {
      config.targetFileType_ = VisualizationTypesParser::parseFileType(
          (*visualizationGeneralConfig)["targetFileType"].get());
    } else {
      std::cout << "# Did not find " << dataSource << "[fileType]. Setting default value "
                << VisualizationTypesParser::toString(defaults.targetFileType_) << "." << std::endl;
      config.targetFileType_ = defaults.targetFileType_;
    }

    config.numBatches_ =
        parseUInt(*visualizationGeneralConfig, "numBatches", defaults.numBatches_, "visualization");
  } else {
    std::cout << "# Could not find specification of visualization[generalConfig]. Falling Back to "
                 "default values."
              << std::endl;
    config = defaults;
  }

  // This is out of the if, since the default value is not directly given by the
  // user
  bool hasFitterCrossvalidationConfig =
      hasFitterConfig() ? (*configFile)[fitter].contains("crossValidation") : false;
  if (hasFitterCrossvalidationConfig) {
    auto crossvalidationConfig = static_cast<DictNode *>(&(*configFile)[fitter]["crossValidation"]);
    config.crossValidation_ =
        parseBool(*crossvalidationConfig, "enable", defaults.crossValidation_, "crossValidation");
  } else {
    config.crossValidation_ = false;
  }
  return hasGeneralConfig;
}

bool DataMiningConfigParser::getVisualizationParameters(
    VisualizationParameters &config, const VisualizationParameters &defaults) const {
  bool hasVisualizationParameters = hasVisualizationParametersConfig();

  if (hasVisualizationParameters) {
    auto visualizationParameters =
        static_cast<DictNode *>(&(*configFile)[visualization]["parameters"]);

    config.perplexity_ =
        parseDouble(*visualizationParameters, "perplexity", defaults.perplexity_, "visualization");

    config.theta_ =
        parseDouble(*visualizationParameters, "theta", defaults.theta_, "visualization");

    config.seed_ = parseUInt(*visualizationParameters, "seed", defaults.seed_, "visualization");

    config.maxNumberIterations_ = parseUInt(*visualizationParameters, "maxNumberIterations",
                                            defaults.maxNumberIterations_, "visualization");

    config.targetDimension_ = parseUInt(*visualizationParameters, "targetDimension",
                                        defaults.targetDimension_, "visualization");

    config.numberCores_ =
        parseUInt(*visualizationParameters, "numberCores", defaults.numberCores_, "visualization");
  } else {
    std::cout << "# Could not find specification of visualization parameters. Falling Back to "
                 "default values."
              << std::endl;
    config = defaults;
  }
  return hasVisualizationParameters;
}

std::string DataMiningConfigParser::parseString(DictNode &dict, const std::string &key,
                                                const std::string &defaultValue,
                                                const std::string &parentDict) const {
  if (dict.contains(key)) {
    try {
      return dict[key].get();
    } catch (json_exception &) {
      try {
        return dict[key]["value"].get();
      } catch (json_exception &) {
        std::string errorMsg = "# Failed to parse string " + parentDict + "[" + key +
                               "] from string";  // + dict[key].get() + ".";
        throw data_exception(errorMsg.c_str());
      }
    }
  } else {
    std::cout << "# Did not find " << parentDict << "[" << key << "]. Setting default value "
              << defaultValue << "." << std::endl;
    return defaultValue;
  }
}

double DataMiningConfigParser::parseDouble(DictNode &dict, const std::string &key,
                                           double defaultValue,
                                           const std::string &parentDict) const {
  if (dict.contains(key)) {
    try {
      return dict[key].getDouble();
    } catch (json_exception &) {
      try {
        if (dict[key].contains("logscale")) {
          if (dict[key]["logscale"].getBool()) {
            return pow(10, dict[key]["value"].getDouble());
          }
        }
        return dict[key]["value"].getDouble();
      } catch (json_exception &) {
        std::string errorMsg = "# Failed to parse double " + parentDict + "[" + key +
                               "] from string";  // + dict[key].get() + ".";
        throw data_exception(errorMsg.c_str());
      }
    }
  } else {
    std::cout << "# Did not find " << parentDict << "[" << key << "]. Setting default value "
              << defaultValue << "." << std::endl;
    return defaultValue;
  }
}

size_t DataMiningConfigParser::parseUInt(DictNode &dict, const std::string &key,
                                         size_t defaultValue, const std::string &parentDict) const {
  if (dict.contains(key)) {
    try {
      return dict[key].getUInt();
    } catch (json_exception &) {
      try {
        return dict[key]["value"].getUInt();
      } catch (json_exception &) {
        std::string errorMsg = "# Failed to parse unsigned integer " + parentDict + "[" + key +
                               "] from string";  // + dict[key].get() + ".";
        throw data_exception(errorMsg.c_str());
      }
    }
  } else {
    std::cout << "# Did not find " << parentDict << "[" << key << "]. Setting default value "
              << defaultValue << "." << std::endl;
    return defaultValue;
  }
}

bool DataMiningConfigParser::parseBool(DictNode &dict, const std::string &key, bool defaultValue,
                                       const std::string &parentNode) const {
  if (dict.contains(key)) {
    try {
      return dict[key].getBool();
    } catch (json_exception &) {
      std::string errorMsg = "# Failed to parse bool " + parentNode + "[" + key + "] from string" +
                             dict[key].get() + ".";
      throw data_exception(errorMsg.c_str());
    }
  } else {
    std::cout << "# Did not find " << parentNode << "[" << key << "]. Setting default value "
              << defaultValue << "." << std::endl;
    return defaultValue;
  }
}

int64_t DataMiningConfigParser::parseInt(DictNode &dict, const std::string &key,
                                         int64_t defaultValue,
                                         const std::string &parentNode) const {
  if (dict.contains(key)) {
    try {
      return dict[key].getInt();
    } catch (json_exception &) {
      try {
        return dict[key]["value"].getInt();
      } catch (json_exception &) {
        std::string errorMsg = "# Failed to parse integer " + parentNode + "[" + key +
                               "] from string";  // + dict[key].get() + ".";
        throw data_exception(errorMsg.c_str());
      }
    }
  } else {
    std::cout << "# Did not find " << parentNode << "[" << key << "]. Setting default value "
              << defaultValue << "." << std::endl;
    return defaultValue;
  }
}

std::vector<int64_t> DataMiningConfigParser::parseIntArray(DictNode &dict, const std::string &key,
                                                           std::vector<int64_t> defaultValue,
                                                           const std::string &parentNode) const {
  if (dict.contains(key)) {
    try {
      std::vector<int64_t> array;
      for (size_t i = 0; i < dict[key].size(); ++i) {
        array.push_back(dict[key][i].getInt());
      }
      return array;
    } catch (json_exception &) {
      std::string errorMsg = "# Failed to parse integer array" + parentNode + "[" + key +
                             "] from string" + dict[key].get() + ".";
      throw data_exception(errorMsg.c_str());
    }
  } else {
    std::cout << "# Did not find " << parentNode << "[" << key << "]. Setting to default value."
              << std::endl;
    return defaultValue;
  }
}

std::vector<std::vector<int64_t>> DataMiningConfigParser::parseArrayOfIntArrays(
    DictNode &dict, const std::string &key, std::vector<std::vector<int64_t>> defaultValue,
    const std::string &parentNode) const {
  if (dict.contains(key)) {
    try {
      std::vector<std::vector<int64_t>> array;
      for (size_t i = 0; i < dict[key].size(); ++i) {
        std::vector<int64_t> entry;
        for (size_t j = 0; j < dict[key][i].size(); j++) {
          entry.push_back(dict[key][i][j].getInt());
        }

        array.push_back(entry);
      }
      return array;
    } catch (json_exception &e) {
      std::string errorMsg = "# Failed to parse array of integer arrays" + parentNode + "[" + key +
                             "] from string" + dict[key].get() + ": " + e.what();
      throw data_exception(errorMsg.c_str());
    }
  } else {
    std::cout << "# Did not find " << parentNode << "[" << key << "]. Setting to default value."
              << std::endl;
    return defaultValue;
  }
}

// (Sebastian) Adapted from parseIntArray without much thought
std::vector<double> DataMiningConfigParser::parseDoubleArray(DictNode &dict, const std::string &key,
                                                             std::vector<double> defaultValue,
                                                             const std::string &parentNode) const {
  if (dict.contains(key)) {
    try {
      std::vector<double> array;
      for (size_t i = 0; i < dict[key].size(); ++i) {
        array.push_back(dict[key][i].getDouble());
      }
      return array;
    } catch (json_exception &) {
      std::string errorMsg = "# Failed to parse double array" + parentNode + "[" + key +
                             "] from string" + dict[key].get() + ".";
      throw data_exception(errorMsg.c_str());
    }
  } else {
    std::cout << "# Did not find " << parentNode << "[" << key << "]. Setting to default value."
              << std::endl;
    return defaultValue;
  }
}

std::vector<std::string> DataMiningConfigParser::parseStringArray(
    json::JSON::DictNode &dict, const std::string &key, std::vector<std::string> defaultValue,
    const std::string &parentNode) const {
  if (dict.contains(key)) {
    try {
      std::vector<std::string> array;
      for (size_t i = 0; i < dict[key].size(); ++i) {
        array.push_back(dict[key][i].get());
      }
      return array;
    } catch (json_exception &) {
      std::string errorMsg = "# Failed to parse double array" + parentNode + "[" + key +
                             "] from string" + dict[key].get() + ".";
      throw data_exception(errorMsg.c_str());
    }
  } else {
    std::cout << "# Did not find " << parentNode << "[" << key << "]. Setting to default value."
              << std::endl;
    return defaultValue;
  }
}

// (Sebastian) Adapted from parseIntArray without much thought
std::vector<size_t> DataMiningConfigParser::parseUIntArray(DictNode &dict, const std::string &key,
                                                           std::vector<size_t> defaultValue,
                                                           const std::string &parentNode) const {
  if (dict.contains(key)) {
    try {
      std::vector<size_t> array;
      for (size_t i = 0; i < dict[key].size(); ++i) {
        array.push_back(dict[key][i].getUInt());
      }
      return array;
    } catch (json_exception &) {
      std::string errorMsg = "# Failed to parse uint array" + parentNode + "[" + key +
                             "] from string" + dict[key].get() + ".";
      throw data_exception(errorMsg.c_str());
    }
  } else {
    std::cout << "# Did not find " << parentNode << "[" << key << "]. Setting to default value."
              << std::endl;
    return defaultValue;
  }
}

void DataMiningConfigParser::parseSLESolverConfig(DictNode &dict, SLESolverConfiguration &config,
                                                  const SLESolverConfiguration &defaults,
                                                  const std::string &parentNode) const {
  config.eps_ = parseDouble(dict, "eps", defaults.eps_, parentNode);
  config.maxIterations_ = parseUInt(dict, "maxIterations", defaults.maxIterations_, parentNode);
  config.threshold_ = parseDouble(dict, "threshold", defaults.threshold_, parentNode);

  // parse  CG type
  if (dict.contains("solverType")) {
    config.type_ = solver::SLESolverTypeParser::parse(dict["solverType"].get());
  } else {
    std::cout << "# Did not find " << parentNode << "[solverType]. Setting default value "
              << solver::SLESolverTypeParser::toString(defaults.type_) << "." << std::endl;
    config.type_ = defaults.type_;
  }
}

void DataMiningConfigParser::getHyperparameters(std::map<std::string, ContinuousParameter> &conpar,
                                                std::map<std::string, DiscreteParameter> &dispar,
                                                std::map<std::string, DiscreteParameter> &catpar,
                                                std::vector<base::GridType> &basisFunctions) const {
  try {
    if ((*configFile)[fitter]["gridConfig"]["level"]["optimize"].getBool()) {
      int64_t min = (*configFile)[fitter]["gridConfig"]["level"]["min"].getInt();
      int64_t max = (*configFile)[fitter]["gridConfig"]["level"]["max"].getInt();
      if (max > min) {
        dispar["level"] = DiscreteParameter("level", static_cast<int>(min), static_cast<int>(max));
      }
    }
  } catch (json_exception &) {
  }
  try {
    if ((*configFile)[fitter]["gridConfig"]["gridType"]["optimize"].getBool()) {
      size_t nOptions = (*configFile)[fitter]["gridConfig"]["gridType"]["options"].size();
      if (nOptions > 1) {
        for (size_t i = 0; i < nOptions; ++i) {
          basisFunctions.push_back(base::GridTypeParser::parse(
              (*configFile)[fitter]["gridConfig"]["gridType"]["options"][i].get()));
        }
        catpar["basisFunction"] =
            DiscreteParameter("basisFunction", 0, static_cast<int>(nOptions - 1));
      }
    }
  } catch (json_exception &) {
  }
  try {
    if ((*configFile)[fitter]["adaptivityConfig"]["noPoints"]["optimize"].getBool()) {
      int64_t min = (*configFile)[fitter]["adaptivityConfig"]["noPoints"]["min"].getInt();
      int64_t max = (*configFile)[fitter]["adaptivityConfig"]["noPoints"]["max"].getInt();
      if (max > min) {
        dispar["noPoints"] =
            DiscreteParameter("noPoints", static_cast<int>(min), static_cast<int>(max));
      }
    }
  } catch (json_exception &) {
  }
  try {
    if ((*configFile)[fitter]["adaptivityConfig"]["threshold"]["optimize"].getBool()) {
      double min = (*configFile)[fitter]["adaptivityConfig"]["threshold"]["min"].getDouble();
      double max = (*configFile)[fitter]["adaptivityConfig"]["threshold"]["max"].getDouble();
      int64_t bits = (*configFile)[fitter]["adaptivityConfig"]["threshold"]["bits"].getInt();
      bool logscale = (*configFile)[fitter]["adaptivityConfig"]["threshold"]["logscale"].getBool();
      if (max > min) {
        conpar["threshold"] =
            ContinuousParameter(static_cast<size_t>(bits), "threshold", min, max, logscale);
      }
    }
  } catch (json_exception &) {
  }
  try {
    if ((*configFile)[fitter]["regularizationConfig"]["lambda"]["optimize"].getBool()) {
      double min = (*configFile)[fitter]["regularizationConfig"]["lambda"]["min"].getDouble();
      double max = (*configFile)[fitter]["regularizationConfig"]["lambda"]["max"].getDouble();
      int64_t bits = (*configFile)[fitter]["regularizationConfig"]["lambda"]["bits"].getInt();
      bool logscale = (*configFile)[fitter]["adaptivityConfig"]["threshold"]["logscale"].getBool();
      if (max > min) {
        conpar["lambda"] =
            ContinuousParameter(static_cast<size_t>(bits), "lambda", min, max, logscale);
      }
    }
  } catch (json_exception &) {
  }
}

void DataMiningConfigParser::getHPOConfig(HPOConfig &config) {
  if (configFile->contains("hpo")) {
    auto node = static_cast<DictNode *>(&(*configFile)["hpo"]);
    config.setSeed(parseInt(*node, "randomSeed", config.getSeed(), "hpo"));
    config.setNTrainSamples(parseInt(*node, "trainSize", config.getNTrainSamples(), "hpo"));
    if (node->contains("harmonica")) {
      auto harmonica = static_cast<DictNode *>(&(*node)["harmonica"]);
      config.setLambda(parseDouble(*harmonica, "lambda", config.getLambda(), "hpo"));
      config.setStages(parseIntArray(*harmonica, "stages", config.getStages(), "hpo"));
      config.setConstraints(
          parseIntArray(*harmonica, "constraints", config.getConstraints(), "hpo"));
    } else {
      std::cout
          << "# Could not find specification of hpo[harmonica]. Falling Back to default values."
          << std::endl;
    }
    if (node->contains("bayesianOptimization")) {
      auto bo = static_cast<DictNode *>(&(*node)["bayesianOptimization"]);
      config.setNRandom(parseInt(*bo, "nRandom", config.getNRandom(), "hpo"));
      config.setNRuns(parseInt(*bo, "nRuns", config.getNRuns(), "hpo"));
    } else {
      std::cout << "# Could not find specification of hpo[bayesianOptimization]. Falling Back to "
                   "default values."
                << std::endl;
    }
  } else {
    std::cout << "# Could not find specification of hpo. Falling Back to default values."
              << std::endl;
  }
}

std::string DataMiningConfigParser::getHPOMethod(std::string defaultValue) const {
  if (configFile->contains("hpo")) {
    auto node = static_cast<DictNode *>(&(*configFile)["hpo"]);
    return parseString(*node, "method", defaultValue, "hpo");
  } else {
    std::cout << "# Did not find hpo[method]. Setting default value " << defaultValue << "."
              << std::endl;
  }
  return defaultValue;
}

void DataMiningConfigParser::parseDataTransformationConfig(DictNode &dict,
                                                           DataTransformationConfig &config,
                                                           const DataTransformationConfig &defaults,
                                                           const std::string &parentNode) const {
  // Parse transformation type
  if (dict.contains("type")) {
    config.type_ = DataTransformationTypeParser::parse(dict["type"].get());
  } else {
    std::cout << "# Did not find [dataTransformationType]. Setting default value "
              << DataTransformationTypeParser::toString(defaults.type_) << "." << std::endl;
    config.type_ = defaults.type_;
  }

  // If type Rosenblatt parse RosenblattTransformationConfig
  if (config.type_ == DataTransformationType::ROSENBLATT) {
    auto rosenblattTransformationConfig = static_cast<DictNode *>(
        &(*configFile)[dataSource]["dataTransformation"]["rosenblattConfig"]);
    parseRosenblattTransformationConfig(*rosenblattTransformationConfig, config.rosenblattConfig_,
                                        defaults.rosenblattConfig_, "rosenblattConfig");
  } else {
    std::cout << "# Could not find specification of "
                 "dataSource[dataTransformationConfig][rosenblattConfig]. Falling back to default "
                 "values."
              << std::endl;
    config.rosenblattConfig_ = defaults.rosenblattConfig_;
  }
}

void DataMiningConfigParser::parseRosenblattTransformationConfig(
    DictNode &dict, RosenblattTransformationConfig &config,
    const RosenblattTransformationConfig &defaults, const std::string &parentNode) const {
  config.numSamples_ = parseUInt(dict, "numSamples", defaults.numSamples_, parentNode);
  config.gridLevel_ = parseUInt(dict, "gridLevel", defaults.gridLevel_, parentNode);

  config.solverMaxIterations_ =
      parseUInt(dict, "solverMaxIterations", defaults.solverMaxIterations_, parentNode);
  config.solverEps_ = parseDouble(dict, "solverEps", defaults.solverEps_, parentNode);
  config.solverThreshold_ =
      parseDouble(dict, "solverThreshold", defaults.solverThreshold_, parentNode);
}

bool DataMiningConfigParser::getFitterDatabaseConfig(
    datadriven::DatabaseConfiguration &config,
    const datadriven::DatabaseConfiguration &defaults) const {
  bool hasDatabaseConfig =
      hasFitterConfig() ? (*configFile)[fitter].contains("databaseConfig") : false;

  if (hasDatabaseConfig) {
    auto databaseConfig = static_cast<DictNode *>(&(*configFile)[fitter]["databaseConfig"]);

    // Parse filepath
    if (databaseConfig->contains("filePath")) {
      config.filePath_ = (*databaseConfig)["filePath"].get();
    } else {
      std::cout << "# Did not find databaseConfig[filepath]. No database loaded" << std::endl;
      config.filePath_ = defaults.filePath_;
    }
  }

  return hasDatabaseConfig;
}

bool DataMiningConfigParser::getFitterLearnerConfig(
    datadriven::LearnerConfiguration &config,
    const datadriven::LearnerConfiguration &defaults) const {
  bool hasLearnerConfig =
      hasFitterConfig() ? (*configFile)[fitter].contains("learnerConfig") : false;

  if (hasLearnerConfig) {
    auto learnerConfig = static_cast<DictNode *>(&(*configFile)[fitter]["learnerConfig"]);

    config.learningRate_ =
        parseDouble(*learnerConfig, "learningRate", defaults.learningRate_, "learnerConfig");
    config.usePrior_ = parseBool(*learnerConfig, "usePrior", defaults.usePrior_, "learnerConfig");
  }

  return hasLearnerConfig;
}

bool DataMiningConfigParser::getGeometryConfig(
    datadriven::GeometryConfiguration &config,
    const datadriven::GeometryConfiguration &defaults) const {
  bool hasGeometryConfig =
      hasFitterConfig() ? (*configFile)[fitter].contains("geometryConfig") : false;

  if (hasGeometryConfig) {
    std::cout << "Has geometry config" << std::endl;
    auto geometryConfig = static_cast<DictNode *>(&(*configFile)[fitter]["geometryConfig"]);

    config.dim_ = parseArrayOfIntArrays(*geometryConfig, "dim", defaults.dim_, "geometryConfig");

    // check if global color available
    int64_t colorIndexDefault = parseInt(*geometryConfig, "colorIndex", -1, "geometryConfig");
    std::vector<size_t> layerDefault;
    for (size_t i = 0; i < config.dim_.size(); i++) {
      layerDefault.push_back(i);
    }

    config.stencils_ = std::vector<sgpp::datadriven::StencilConfiguration>();

    if ((*geometryConfig).contains("stencils")) {
      size_t nStencils = (*geometryConfig)["stencils"].size();
      for (size_t i = 0; i < nStencils; ++i) {
        StencilConfiguration stencil;
        auto stencilConfig = static_cast<DictNode *>(&(*geometryConfig)["stencils"][i]);
        stencil.applyOnLayers_ =
            parseUIntArray(*stencilConfig, "applyOnLayers", layerDefault, "stencils");
        stencil.colorIndex_ = parseInt(*stencilConfig, "colorIndex", colorIndexDefault, "stencils");
        stencil.stencilType_ = GeometryConfigurationParser::parseStencil(
            (*geometryConfig)["stencils"][i]["stencil"].get());
        if (stencil.stencilType_ == sgpp::datadriven::StencilType::Block) {
          stencil.blockLenght_ = parseInt(*stencilConfig, "blockLenght", 2, "stencils");
        } else {
          stencil.blockLenght_ = 0;
        }
        config.stencils_.push_back(stencil);
      }
    } else {
      config.stencils_ = defaults.stencils_;
    }

    // Validate configuration
    size_t numberOfLayers = config.dim_.size();
    if (numberOfLayers > 0) {
      size_t numberOfAxes = config.dim_[0].size();
      for (const std::vector<int64_t> &res : config.dim_) {
        if (numberOfAxes != res.size())
          throw data_exception("Each layer has to have identical number of axes");
      }
      for (const sgpp::datadriven::StencilConfiguration &stencilConf : config.stencils_) {
        if (stencilConf.colorIndex_ != -1 &&
            static_cast<size_t>(stencilConf.colorIndex_) >= numberOfAxes) {
          throw data_exception("ColorIndex is not a valid index for an axis:");
        }
        for (size_t layerIndex : stencilConf.applyOnLayers_) {
          if (layerIndex >= numberOfAxes) {
            throw data_exception("There is an invalid index contained in ApplyOnLayers");
          }
        }
      }
    }
  }
  return hasGeometryConfig;
}

bool DataMiningConfigParser::getFitterParallelConfig(
    datadriven::ParallelConfiguration &config,
    const datadriven::ParallelConfiguration &defaults) const {
  bool hasParallelConfig =
      hasFitterConfig() ? (*configFile)[fitter].contains("parallelConfig") : false;

  if (hasParallelConfig) {
    auto parallelConfig = static_cast<DictNode *>(&(*configFile)[fitter]["parallelConfig"]);

    config.scalapackEnabled_ = true;

    config.processRows_ = static_cast<int>(
        parseInt(*parallelConfig, "processRows", defaults.processRows_, "parallelConfig"));
    config.processCols_ = static_cast<int>(
        parseInt(*parallelConfig, "processColumns", defaults.processCols_, "parallelConfig"));

    config.rowBlockSize_ =
        parseUInt(*parallelConfig, "rowBlockSize", defaults.rowBlockSize_, "parallelConfig");
    config.columnBlockSize_ =
        parseUInt(*parallelConfig, "columnBlockSize", defaults.columnBlockSize_, "parallelConfig");
  }

  return hasParallelConfig;
}

} /* namespace datadriven */
} /* namespace sgpp */
