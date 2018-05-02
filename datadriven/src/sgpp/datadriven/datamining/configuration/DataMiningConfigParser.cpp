/* Copyright (C) 2008-today The SG++ project
 * This file is part of the SG++ project. For conditions of distribution and
 * use, please see the copyright notice provided with SG++ or at
 * sgpp.sparsegrids.org
 *
 * DataMiningConfigParser.cpp
 *
 *  Created on: Aug 14, 2016
 *  	Author: Michael Lettrich
 */

#include <sgpp/base/exception/data_exception.hpp>
#include <sgpp/base/exception/file_exception.hpp>
#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/grid/LevelIndexTypes.hpp>
#include <sgpp/base/tools/json/JSON.hpp>
#include <sgpp/base/tools/json/json_exception.hpp>
#include <sgpp/datadriven/configuration/RegularizationConfiguration.hpp>
#include <sgpp/datadriven/datamining/configuration/DataMiningConfigParser.hpp>
#include <sgpp/datadriven/datamining/configuration/DensityEstimationTypeParser.hpp>
#include <sgpp/datadriven/datamining/configuration/GridTypeParser.hpp>
#include <sgpp/datadriven/datamining/configuration/MatrixDecompositionTypeParser.hpp>
#include <sgpp/datadriven/datamining/configuration/RegularizationTypeParser.hpp>
#include <sgpp/datadriven/datamining/configuration/SLESolverTypeParser.hpp>
#include <sgpp/datadriven/datamining/modules/dataSource/DataSourceFileTypeParser.hpp>
#include <sgpp/datadriven/datamining/modules/dataSource/DataTransformationTypeParser.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/FitterTypeParser.hpp>
#include <sgpp/datadriven/datamining/modules/scoring/ScorerMetricTypeParser.hpp>
#include <sgpp/datadriven/datamining/modules/scoring/ScorerShufflingTypeParser.hpp>
#include <sgpp/solver/TypesSolver.hpp>

#include <string>

using json::JSON;
using json::json_exception;
using sgpp::base::file_exception;
using sgpp::base::data_exception;
using json::DictNode;
using sgpp::solver::SLESolverConfiguration;

namespace sgpp {
namespace datadriven {

const std::string DataMiningConfigParser::dataSource = "dataSource";
const std::string DataMiningConfigParser::scorer = "scorer";
const std::string DataMiningConfigParser::fitter = "fitter";

DataMiningConfigParser::DataMiningConfigParser(const std::string& filepath) : configFile(nullptr) {
  try {
    configFile = std::make_unique<JSON>(filepath);
  } catch (json_exception& exception) {
    std::cout << exception.what() << std::endl;
    std::string errorMsg = "can not open file: \"" + filepath + "\"";
    throw file_exception(errorMsg.c_str());
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

bool DataMiningConfigParser::hasScorerConfigCrossValidation() const {
  bool hasScorerCrossValidationConfig =
      hasScorerConfig() ? (*configFile)[scorer].contains("crossValidation") : false;
  return hasScorerCrossValidationConfig;
}

bool DataMiningConfigParser::hasScorerConfigTesting() const {
  bool hasScorerTestingConfig =
      hasScorerConfig() ? (*configFile)[scorer].contains("testing") : false;
  return hasScorerTestingConfig;
}

bool DataMiningConfigParser::hasFitterConfig() const { return configFile->contains(fitter); }

bool DataMiningConfigParser::getDataSourceConfig(DataSourceConfig& config,
                                                 const DataSourceConfig& defaults) const {
  bool hasDataSource = hasDataSourceConfig();

  if (hasDataSource) {
    auto dataSourceConfig = static_cast<DictNode*>(&(*configFile)[dataSource]);

    config.filePath = parseString(*dataSourceConfig, "filePath", defaults.filePath, "dataSource");
    config.isCompressed =
        parseBool(*dataSourceConfig, "compression", defaults.isCompressed, "dataSource");
    config.numBatches =
        parseUInt(*dataSourceConfig, "numBatches", defaults.numBatches, "dataSource");
    config.batchSize = parseUInt(*dataSourceConfig, "batchSize", defaults.batchSize, "dataSource");

    // parse file type
    if (dataSourceConfig->contains("fileType")) {
      config.fileType = DataSourceFileTypeParser::parse((*dataSourceConfig)["fileType"].get());
    } else {
      std::cout << "# Did not find " << dataSource << "[fileType]. Setting default value "
                << DataSourceFileTypeParser::toString(defaults.fileType) << "." << std::endl;
      config.fileType = defaults.fileType;
    }

    // parse dataTransformationConfig
    bool hasDataTransformation = hasDataTransformationConfig();

    if (hasDataTransformation) {
      auto dataTransformationConfig =
          static_cast<DictNode*>(&(*configFile)[dataSource]["dataTransformation"]);
      parseDataTransformationConfig(*dataTransformationConfig, config.dataTransformationConfig,
          defaults.dataTransformationConfig, "dataTransformation");
    } else {
      std::cout << "# Could not find specification of dataSource[dataTransformationConfig]. "
          "Falling back to default values." << std::endl;
      config.dataTransformationConfig = defaults.dataTransformationConfig;
    }

  } else {
    std::cout << "# Could not find specification of dataSource. Falling Back to default values."
              << std::endl;
    config = defaults;
  }
  return hasDataSource;
}

bool DataMiningConfigParser::getScorerTestingConfig(TestingConfiguration& config,
                                                    const TestingConfiguration& defaults) const {
  bool hasScorerTestingConfig = hasScorerConfigTesting();

  if (hasScorerTestingConfig) {
    auto scorerTestingConfig = static_cast<DictNode*>(&(*configFile)[scorer]["testing"]);

    config.testingPortion =
        parseDouble(*scorerTestingConfig, "testingPortion", defaults.testingPortion, "testing");
    // parse shuffling type
    if (scorerTestingConfig->contains("shuffling")) {
      config.shuffling =
          ScorerShufflingTypeParser::parse((*scorerTestingConfig)["shuffling"].get());
    } else {
      std::cout << "# Did not find testing[shuffling]. Setting default value "
                << ScorerShufflingTypeParser::toString(defaults.shuffling) << "." << std::endl;
      config.shuffling = defaults.shuffling;
    }

    config.randomSeed =
        parseInt(*scorerTestingConfig, "randomSeed", defaults.randomSeed, "testing");
    // parse metric type
    if (scorerTestingConfig->contains("metric")) {
      config.metric = ScorerMetricTypeParser::parse((*scorerTestingConfig)["metric"].get());
    } else {
      std::cout << "# Did not find testing[metric]. Setting default value "
                << ScorerMetricTypeParser::toString(defaults.metric) << "." << std::endl;
    }

  } else {
    std::cout
        << "# Could not find specification  of scorer[testing]. Falling Back to default values."
        << std::endl;
    config = defaults;
  }
  return hasScorerTestingConfig;
}

bool DataMiningConfigParser::getScorerCrossValidationConfig(
    CrossValidationConfiguration& config, const CrossValidationConfiguration& defaults) const {
  bool hasScorerCrossValidationConfig = hasScorerConfigCrossValidation();

  if (hasScorerCrossValidationConfig) {
    auto scorerTestingConfig = static_cast<DictNode*>(&(*configFile)[scorer]["crossValidation"]);

    config.folds = parseUInt(*scorerTestingConfig, "folds", defaults.folds, "crossValidation");
    // parse shuffling type
    if (scorerTestingConfig->contains("shuffling")) {
      config.shuffling =
          ScorerShufflingTypeParser::parse((*scorerTestingConfig)["shuffling"].get());
    } else {
      std::cout << "# Did not find crossValidation[shuffling]. Setting default value "
                << ScorerShufflingTypeParser::toString(defaults.shuffling) << "." << std::endl;
      config.shuffling = defaults.shuffling;
    }
    config.randomSeed =
        parseInt(*scorerTestingConfig, "randomSeed", defaults.randomSeed, "crossValidation");
    // parse metric type
    if (scorerTestingConfig->contains("metric")) {
      config.metric = ScorerMetricTypeParser::parse((*scorerTestingConfig)["metric"].get());
    } else {
      std::cout << "# Did not find crossValidation[metric]. Setting default value "
                << ScorerMetricTypeParser::toString(defaults.metric) << "." << std::endl;
    }

  } else {
    std::cout << "# Could not find specification  of scorer[crossValidation]. Falling Back to "
                 "default values."
              << std::endl;
    config = defaults;
  }
  return hasScorerCrossValidationConfig;
}

// TODO(lettrich): is this consistent with the rest of the parsing?
bool DataMiningConfigParser::getFitterConfigType(FitterType& config,
                                                 const FitterType& defaults) const {
  bool hasFitterConfig = this->hasFitterConfig();

  if (hasFitterConfig) {
    auto fitterConfig = static_cast<DictNode*>(&(*configFile)[fitter]);
    if (fitterConfig->contains("type")) {
      config = FitterTypeParser::parse((*fitterConfig)["type"].get());
    } else {
      std::cout << "# Could not find specification  of fitter[type]. Falling Back to default "
                   "values."
                << std::endl;
      config = defaults;
    }
  }

  return hasFitterConfig;
}

bool DataMiningConfigParser::getFitterGridConfig(RegularGridConfiguration& config,
                                                 const RegularGridConfiguration& defaults) const {
  bool hasFitterGridConfig =
      hasFitterConfig() ? (*configFile)[fitter].contains("gridConfig") : false;

  if (hasFitterGridConfig) {
    auto fitterConfig = static_cast<DictNode*>(&(*configFile)[fitter]["gridConfig"]);
    config.dim_ = parseUInt(*fitterConfig, "dim", defaults.dim_, "gridConfig");
    config.level_ =
        static_cast<int>(parseInt(*fitterConfig, "level", defaults.level_, "gridConfig"));
    config.maxDegree_ = parseUInt(*fitterConfig, "maxDegree", defaults.maxDegree_, "gridConfig");
    config.boundaryLevel_ = static_cast<unsigned int>(
        parseUInt(*fitterConfig, "boundaryLevel", defaults.boundaryLevel_, "gridConfig"));
    config.filename_ = parseString(*fitterConfig, "fileName", defaults.filename_, "gridConfig");

    // parse  grid type
    if (fitterConfig->contains("gridType")) {
      config.type_ = GridTypeParser::parse((*fitterConfig)["gridType"].get());
    } else {
      std::cout << "# Did not find gridConfig[gridType]. Setting default value "
                << GridTypeParser::toString(defaults.type_) << "." << std::endl;
      config.type_ = defaults.type_;
    }

  } else {
    std::cout << "# Could not find specification  of fitter[gridConfig]. Falling Back to default "
                 "values."
              << std::endl;
    config = defaults;
  }
  return hasFitterGridConfig;
}

bool DataMiningConfigParser::getFitterAdaptivityConfig(
    AdpativityConfiguration& config, const AdpativityConfiguration& defaults) const {
  bool hasFitterAdaptivityConfig =
      hasFitterConfig() ? (*configFile)[fitter].contains("adaptivityConfig") : false;

  if (hasFitterAdaptivityConfig) {
    auto adaptivityConfig = static_cast<DictNode*>(&(*configFile)[fitter]["adaptivityConfig"]);
    config.numRefinements_ = parseUInt(*adaptivityConfig, "numRefinements",
                                       defaults.numRefinements_, "adaptivityConfig");
    config.threshold_ =
        parseDouble(*adaptivityConfig, "threshold", defaults.threshold_, "adaptivityConfig");
    config.maxLevelType_ =
        parseBool(*adaptivityConfig, "maxLevelType", defaults.maxLevelType_, "adaptivityConfig");
    config.noPoints_ =
        parseUInt(*adaptivityConfig, "noPoints", defaults.noPoints_, "adaptivityConfig");
    config.percent_ =
        parseDouble(*adaptivityConfig, "percent", defaults.percent_, "adaptivityConfig");
    config.errorBasedRefinement = parseBool(*adaptivityConfig, "errorBasedRefinement",
                                            defaults.errorBasedRefinement, "adaptivityConfig");
  } else {
    std::cout << "# Could not find specification  of fitter[adaptivityConfig]. Falling Back to "
                 "default values."
              << std::endl;
    config = defaults;
  }
  return hasFitterAdaptivityConfig;
}

bool DataMiningConfigParser::getFitterCrossvalidationConfig(
    CrossvalidationConfiguration& config, const CrossvalidationConfiguration& defaults) const {
  bool hasFitterCrossvalidationConfig =
      hasFitterConfig() ? (*configFile)[fitter].contains("crossvalidationConfig") : false;

  if (hasFitterCrossvalidationConfig) {
    auto crossvalidationConfig =
    static_cast<DictNode*>(&(*configFile)[fitter]["crossvalidationConfig"]);
    config.enable_ =
        parseBool(*crossvalidationConfig, "enable", defaults.enable_, "crossvalidationConfig");
    config.kfold_ = parseUInt(*crossvalidationConfig, "kFold",
                                       defaults.kfold_, "crossvalidationConfig");
    config.seed_ =
        static_cast<int>(parseInt(*crossvalidationConfig, "randomSeed",
                                  defaults.seed_, "crossvalidationConfig"));
    config.shuffle_ =
        parseBool(*crossvalidationConfig, "shuffle", defaults.shuffle_, "crossvalidationConfig");
    config.silent_ =
        parseBool(*crossvalidationConfig, "silent", defaults.silent_, "crossvalidationConfig");
    config.lambda_ =
        parseDouble(*crossvalidationConfig, "lambda", defaults.lambda_, "crossvalidationConfig");
    config.lambdaStart_ =
        parseDouble(*crossvalidationConfig, "lambdaStart",
                    defaults.lambdaStart_, "crossvalidationConfig");
    config.lambdaEnd_ =
        parseDouble(*crossvalidationConfig, "lambdaEnd",
                    defaults.lambdaEnd_, "crossvalidationConfig");
    config.lambdaSteps_ =
        parseUInt(*crossvalidationConfig, "lambdaSteps",
                  defaults.lambdaSteps_, "crossvalidationConfig");
    config.logScale_ =
        parseBool(*crossvalidationConfig, "logScale", defaults.logScale_, "crossvalidationConfig");
  } else {
    std::cout << "# Could not find specification  of fitter[crossvalidationConfig]. Falling "
                 "Back to default values."
              << std::endl;
    config = defaults;
  }
  return hasFitterCrossvalidationConfig;
}

bool DataMiningConfigParser::getFitterDensityEstimationConfig(
    DensityEstimationConfiguration& config, const DensityEstimationConfiguration& defaults) const {
  bool hasFitterDensityEstimationConfig =
      hasFitterConfig() ? (*configFile)[fitter].contains("densityEstimationConfig") : false;

  if (hasFitterDensityEstimationConfig) {
    auto densityEstimationConfig =
    static_cast<DictNode*>(&(*configFile)[fitter]["densityEstimationConfig"]);
    config.iCholSweepsDecompose_ = parseUInt(*densityEstimationConfig,
                                             "iCholSweepsDecompose",
                                             defaults.iCholSweepsDecompose_,
                                             "densityEstimationConfig");
    config.iCholSweepsRefine_ = parseUInt(*densityEstimationConfig, "iCholSweepsRefine",
                                          defaults.iCholSweepsRefine_, "densityEstimationConfig");
    config.iCholSweepsUpdateLambda_ =
    parseUInt(*densityEstimationConfig, "iCholSweepsUpdateLambda",
              defaults.iCholSweepsUpdateLambda_, "densityEstimationConfig");
    config.iCholSweepsSolver_ = parseUInt(*densityEstimationConfig, "iCholSweepsSolver",
                                          defaults.iCholSweepsSolver_, "densityEstimationConfig");

    // parse  density estimation type
    if (densityEstimationConfig->contains("densityEstimationType")) {
      config.type_ =
      DensityEstimationTypeParser::parse((*densityEstimationConfig)["densityEstimationType"].get());
    } else {
      std::cout <<
      "# Did not find densityEstimationConfig[densityEstimationType]. Setting default value "
                << DensityEstimationTypeParser::toString(defaults.type_) << "." << std::endl;
      config.type_ = defaults.type_;
    }

    // parse matrix decomposition type
    if (densityEstimationConfig->contains("matrixDecompositionType")) {
      config.decomposition_ =
      MatrixDecompositionTypeParser::
      parse((*densityEstimationConfig)["matrixDecompositionEstimationType"].get());
    } else {
      std::cout <<
      "# Did not find densityEstimationConfig[matrixDecompositionType]. Setting default value "
                << MatrixDecompositionTypeParser::toString(defaults.decomposition_)
                << "." << std::endl;
      config.decomposition_ = defaults.decomposition_;
    }

  } else {
    std::cout <<
    "# Could not find specification  of fitter[densityEstimationConfig]. Falling Back to "
                 "default values."
              << std::endl;
    config = defaults;
  }
  return hasFitterDensityEstimationConfig;
}

bool DataMiningConfigParser::getFitterSolverRefineConfig(
    SLESolverConfiguration& config, const SLESolverConfiguration& defaults) const {
  bool hasFitterSolverRefineConfig =
      hasFitterConfig() ? (*configFile)[fitter].contains("solverRefineConfig") : false;

  if (hasFitterSolverRefineConfig) {
    auto solverConfig = static_cast<DictNode*>(&(*configFile)[fitter]["solverRefineConfig"]);

    parseSLESolverConfig(*solverConfig, config, defaults, "solverRefineConfig");
  } else {
    std::cout << "# Could not find specification  of fitter[solverRefineConfig]. Falling Back to "
                 "default values."
              << std::endl;
    config = defaults;
  }
  return hasFitterSolverRefineConfig;
}

bool DataMiningConfigParser::getFitterSolverFinalConfig(
    SLESolverConfiguration& config, const SLESolverConfiguration& defaults) const {
  bool hasFitterSolverFinalConfig =
      hasFitterConfig() ? (*configFile)[fitter].contains("solverFinalConfig") : false;

  if (hasFitterSolverFinalConfig) {
    auto solverConfig = static_cast<DictNode*>(&(*configFile)[fitter]["solverFinalConfig"]);

    parseSLESolverConfig(*solverConfig, config, defaults, "solverFinalConfig");
  } else {
    std::cout << "# Could not find specification  of fitter[solverFinalConfig]. Falling Back to "
                 "default values."
              << std::endl;
    config = defaults;
  }
  return hasFitterSolverFinalConfig;
}

bool DataMiningConfigParser::getFitterRegularizationConfig(
    RegularizationConfiguration& config, const RegularizationConfiguration& defaults) const {
  bool hasRegularizationConfig =
      hasFitterConfig() ? (*configFile)[fitter].contains("regularizationConfig") : false;

  if (hasRegularizationConfig) {
    auto regularizationConfig =
        static_cast<DictNode*>(&(*configFile)[fitter]["regularizationConfig"]);

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
  }

  return hasRegularizationConfig;
}

std::string DataMiningConfigParser::parseString(DictNode& dict, const std::string& key,
                                                const std::string& defaultValue,
                                                const std::string& parentDict) const {
  if (dict.contains(key)) {
    try {
      return dict[key].get();
    } catch (json_exception& e) {
      std::string errorMsg = "# Failed to parse string " + parentDict + "[" + key +
                             "] from string" + dict[key].get() + ".";
      throw data_exception(errorMsg.c_str());
    }
  } else {
    std::cout << "# Did not find " << parentDict << "[" << key << "]. Setting default value "
              << defaultValue << "." << std::endl;
    return defaultValue;
  }
}

double DataMiningConfigParser::parseDouble(DictNode& dict, const std::string& key,
                                           double defaultValue,
                                           const std::string& parentDict) const {
  if (dict.contains(key)) {
    try {
      return dict[key].getDouble();
    } catch (json_exception& e) {
      std::string errorMsg = "# Failed to parse double " + parentDict + "[" + key +
                             "] from string" + dict[key].get() + ".";
      throw data_exception(errorMsg.c_str());
    }
  } else {
    std::cout << "# Did not find " << parentDict << "[" << key << "]. Setting default value "
              << defaultValue << "." << std::endl;
    return defaultValue;
  }
}

size_t DataMiningConfigParser::parseUInt(DictNode& dict, const std::string& key,
                                         size_t defaultValue, const std::string& parentDict) const {
  if (dict.contains(key)) {
    try {
      return dict[key].getUInt();
    } catch (json_exception& e) {
      std::string errorMsg = "# Failed to parse unsigned integer " + parentDict + "[" + key +
                             "] from string" + dict[key].get() + ".";
      throw data_exception(errorMsg.c_str());
    }
  } else {
    std::cout << "# Did not find " << parentDict << "[" << key << "]. Setting default value "
              << defaultValue << "." << std::endl;
    return defaultValue;
  }
}

bool DataMiningConfigParser::parseBool(DictNode& dict, const std::string& key, bool defaultValue,
                                       const std::string& parentNode) const {
  if (dict.contains(key)) {
    try {
      return dict[key].getBool();
    } catch (json_exception& e) {
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

int64_t DataMiningConfigParser::parseInt(DictNode& dict, const std::string& key,
                                         int64_t defaultValue,
                                         const std::string& parentNode) const {
  if (dict.contains(key)) {
    try {
      return dict[key].getInt();
    } catch (json_exception& e) {
      std::string errorMsg = "# Failed to parse integer " + parentNode + "[" + key +
                             "] from string" + dict[key].get() + ".";
      throw data_exception(errorMsg.c_str());
    }
  } else {
    std::cout << "# Did not find " << parentNode << "[" << key << "]. Setting default value "
              << defaultValue << "." << std::endl;
    return defaultValue;
  }
}

void DataMiningConfigParser::parseSLESolverConfig(DictNode& dict, SLESolverConfiguration& config,
                                                  const SLESolverConfiguration& defaults,
                                                  const std::string& parentNode) const {
  config.eps_ = parseDouble(dict, "eps", defaults.eps_, parentNode);
  config.maxIterations_ = parseUInt(dict, "maxIterations", defaults.maxIterations_, parentNode);
  config.threshold_ = parseDouble(dict, "threshold", defaults.threshold_, parentNode);

  // parse  CG type
  if (dict.contains("solverType")) {
    config.type_ = SLESolverTypeParser::parse(dict["solverType"].get());
  } else {
    std::cout << "# Did not find " << parentNode << "[solverType]. Setting default value "
              << SLESolverTypeParser::toString(defaults.type_) << "." << std::endl;
    config.type_ = defaults.type_;
  }
}

void DataMiningConfigParser::parseDataTransformationConfig(
    DictNode& dict, DataTransformationConfig& config,
    const DataTransformationConfig& defaults,
    const std::string& parentNode) const {

  // Parse transformation type
  if (dict.contains("type")) {
    config.type = DataTransformationTypeParser::parse(dict["type"].get());
  } else {
    std::cout << "# Did not find [dataTransformationType]. Setting default value "
        << DataTransformationTypeParser::toString(defaults.type) << "." << std::endl;
    config.type = defaults.type;
  }

  // If type Rosenblatt parse RosenblattTransformationConfig
  if (config.type == DataTransformationType::ROSENBLATT) {
    auto rosenblattTransformationConfig = static_cast<DictNode*>
       (&(*configFile)[dataSource]["dataTransformation"]["rosenblattConfig"]);
    parseRosenblattTransformationConfig(*rosenblattTransformationConfig,
        config.rosenblattConfig, defaults.rosenblattConfig, "rosenblattConfig");
  } else {
    std::cout << "# Could not find specification of dataSource[dataTransformationConfig]"
        "[rosenblattConfig]. Falling back to default values." << std::endl;
    config.rosenblattConfig = defaults.rosenblattConfig;
  }
}

void DataMiningConfigParser::parseRosenblattTransformationConfig(
    DictNode& dict, RosenblattTransformationConfig& config,
    const RosenblattTransformationConfig& defaults,
    const std::string& parentNode) const {

  config.numSamples = parseUInt(dict, "numSamples", defaults.numSamples, parentNode);
  config.gridLevel = parseUInt(dict, "gridLevel", defaults.gridLevel, parentNode);

  config.solverMaxIterations = parseUInt(dict, "solverMaxIterations",
      defaults.solverMaxIterations, parentNode);
  config.solverEps = parseDouble(dict, "solverEps",
      defaults.solverEps, parentNode);
  config.solverThreshold = parseDouble(dict, "solverThreshold",
      defaults.solverThreshold, parentNode);
}


bool DataMiningConfigParser::getFitterDatabaseConfig(
    datadriven::DatabaseConfiguration config, const datadriven::DatabaseConfiguration& defaults)
const {
  bool hasDatabaseConfig =
      hasFitterConfig() ? (*configFile)[fitter].contains("database") : false;
  auto databaseConfig = static_cast<DictNode*>(&(*configFile)[fitter]["database"]);

  // Parse filepath
  if (databaseConfig->contains("filepath")) {
    config.filepath = (*databaseConfig)["filepath"].get();
  } else {
    std::cout << "# Did not find databaseConfig[filepath]. No database loaded" << std::endl;
    config.filepath = defaults.filepath;
  }

  return hasDatabaseConfig;
}
} /* namespace datadriven */
} /* namespace sgpp */
