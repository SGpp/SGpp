/*
 * DataMiningParser.cpp
 *
 *  Created on: Aug 14, 2016
 *      Author: michael
 */

#include "DataMiningConfigParser.hpp"
#include <sgpp/base/exception/file_exception.hpp>
#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/tools/GridTypeParser.hpp>
#include <sgpp/base/tools/json/JSON.hpp>
#include <sgpp/base/tools/json/json_exception.hpp>
#include <sgpp/datadriven/application/RegularizationConfiguration.hpp>
#include <sgpp/solver/TypesSolver.hpp>

#include <string>

using json::JSON;
using json::json_exception;
using sgpp::base::file_exception;
using json::DictNode;
using sgpp::base::GridTypeParser;
using sgpp::solver::SLESolverConfiguration;
using sgpp::solver::SolverTypeParser;

namespace sgpp {
namespace datadriven {

const std::string DataMiningConfigParser::dataSource = "dataSource";
const std::string DataMiningConfigParser::scorer = "scorer";
const std::string DataMiningConfigParser::fitter = "fitter";

DataMiningConfigParser::DataMiningConfigParser(const std::string& filepath) : configFile(nullptr) {
  try {
    configFile = std::make_shared<JSON>(filepath);
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

bool DataMiningConfigParser::hasScorerConfig() const { return configFile->contains(scorer); }

bool DataMiningConfigParser::hasFitterConfig() const { return configFile->contains(fitter); }

bool DataMiningConfigParser::getDataSourceConfig(DataSourceConfig& config,
                                                 const DataSourceConfig& defaults) const {
  bool hasDataSource = hasDataSourceConfig();

  if (hasDataSource) {
    auto parser = DataSourceFileTypeParser();
    auto dataSourceConfig = static_cast<DictNode*>(&(*configFile)[dataSource]);

    config.filePath = parseString(*dataSourceConfig, "filePath", defaults.filePath, "dataSource");
    config.isCompressed =
        parseBool(*dataSourceConfig, "compression", defaults.isCompressed, "dataSource");
    config.numBatches =
        parseUint(*dataSourceConfig, "numBatches", defaults.numBatches, "dataSource");
    config.batchSize = parseUint(*dataSourceConfig, "batchSize", defaults.batchSize, "dataSource");

    // parse file type
    if (dataSourceConfig->contains("fileType")) {
      auto parser = DataSourceFileTypeParser();
      config.fileType = parser((*dataSourceConfig)["fileType"].get());
    }

  } else {
    std::cout << "# Could not find specification  of dataSource. Falling Back to default values."
              << std::endl;
    config = defaults;
  }
  return hasDataSource;
}

bool DataMiningConfigParser::getScorerTestingConfig(TestingConfiguration& config,
                                                    const TestingConfiguration& defaults) const {
  bool hasScorerTestingConfig = hasScorerConfig() ? configFile->contains("testing") : false;

  if (hasScorerTestingConfig) {
    auto scorerTestingConfig = static_cast<DictNode*>(&(*configFile)[scorer]);

    config.testingPortion =
        parseDouble(*scorerTestingConfig, "testingPortion", defaults.testingPortion, "testing");
    // parse shuffling type
    if (scorerTestingConfig->contains("shuffling")) {
      config.shuffling =
          ScorerShufflingTypeParser::parse((*scorerTestingConfig)["shuffling"].get());
    }
    config.randomSeed = static_cast<int64_t>(
        parseInt(*scorerTestingConfig, "randomSeed", defaults.randomSeed, "testing"));
    // parse metric type
    if (scorerTestingConfig->contains("shuffling")) {
      config.metric = ScorerMetricParser::parse((*scorerTestingConfig)["metric"].get());
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
  bool hasScorerCrossValidationConfig =
      hasScorerConfig() ? configFile->contains("crossValidation") : false;

  if (hasScorerCrossValidationConfig) {
    auto scorerTestingConfig = static_cast<DictNode*>(&(*configFile)[scorer]);

    config.folds =
        parseUint(*scorerTestingConfig, "testingPortion", defaults.folds, "crossValidation");
    // parse shuffling type
    if (scorerTestingConfig->contains("shuffling")) {
      config.shuffling =
          ScorerShufflingTypeParser::parse((*scorerTestingConfig)["shuffling"].get());
    }
    config.randomSeed =
        parseInt(*scorerTestingConfig, "randomSeed", defaults.randomSeed, "crossValidation");
    // parse metric type
    if (scorerTestingConfig->contains("shuffling")) {
      config.metric = ScorerMetricParser::parse((*scorerTestingConfig)["metric"].get());
    }

  } else {
    std::cout << "# Could not find specification  of scorer[crossValidation]. Falling Back to "
                 "default values."
              << std::endl;
    config = defaults;
  }
  return hasScorerCrossValidationConfig;
}

bool DataMiningConfigParser::getFitterConfigType(FitterType& config,
                                                 const FitterType& defaults) const {
  bool hasFitterConfig = this->hasFitterConfig();

  if (hasFitterConfig) {
    auto fitterConfig = static_cast<DictNode*>(&(*configFile)[fitter]);
    if (fitterConfig->contains("type")) {
      config = FitterTypeParser::parse((*fitterConfig)["type"].get());
    } else {
      config = defaults;
    }
  } else {
    std::cout << "# Could not find specification  of fitter[type]. Falling Back to default "
                 "values."
              << std::endl;
    config = defaults;
  }

  return hasFitterConfig;
}

bool DataMiningConfigParser::getFitterGridConfig(RegularGridConfiguration& config,
                                                 const RegularGridConfiguration& defaults) const {
  bool hasFitterGridConfig = hasFitterConfig() ? configFile->contains("gridConfig") : false;

  if (hasFitterGridConfig) {
    auto fitterConfig = static_cast<DictNode*>(&(*configFile)[fitter]["gridConfig"]);
    config.dim_ = parseUint(*fitterConfig, "dim", defaults.dim_, "gridConfig");
    config.level_ = parseInt(*fitterConfig, "level", defaults.level_, "gridConfig");
    config.maxDegree_ = parseUint(*fitterConfig, "maxDegree", defaults.maxDegree_, "gridConfig");
    config.boundaryLevel_ =
        parseUint(*fitterConfig, "boundaryLevel", defaults.boundaryLevel_, "gridConfig");
    config.filename_ = parseString(*fitterConfig, "fileName", defaults.filename_, "gridConfig");

    // parse  grid type
    if (fitterConfig->contains("gridType")) {
      auto parser = GridTypeParser();
      config.type_ = parser((*fitterConfig)["gridType"].get());
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
      hasFitterConfig() ? configFile->contains("adaptivityConfig") : false;

  if (hasFitterAdaptivityConfig) {
    auto adaptivityConfig = static_cast<DictNode*>(&(*configFile)[fitter]["adaptivityConfig"]);
    config.numRefinements_ = parseUint(*adaptivityConfig, "numRefinements",
                                       defaults.numRefinements_, "adaptivityConfig");
    config.threshold_ =
        parseDouble(*adaptivityConfig, "threshold", defaults.threshold_, "adaptivityConfig");
    config.maxLevelType_ =
        parseBool(*adaptivityConfig, "maxLevelType", defaults.maxLevelType_, "adaptivityConfig");
    config.noPoints_ =
        parseUint(*adaptivityConfig, "noPoints", defaults.noPoints_, "adaptivityConfig");
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

bool DataMiningConfigParser::getFitterSolverRefineConfig(
    SLESolverConfiguration& config, const SLESolverConfiguration& defaults) const {
  bool hasFitterSolverRefineConfig =
      hasFitterConfig() ? configFile->contains("solverRefineConfig") : false;

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
      hasFitterConfig() ? configFile->contains("solverFinalConfig") : false;

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
      hasFitterConfig() ? configFile->contains("regularizationConfig") : false;

  if (hasRegularizationConfig) {
    auto regularizationConfig =
        static_cast<DictNode*>(&(*configFile)[fitter]["regularizationConfig"]);

    // parse  grid type
    if (regularizationConfig->contains("regularizationType")) {
      auto parser = RegularizationTypeParser();
      config.regType_ = parser((*regularizationConfig)["gridType"].get());
    }
  }

  return hasRegularizationConfig;
}

bool DataMiningConfigParser::getFitterLambda(double& lambda, double defaultValue) const {
  bool hasRegularizationConfig =
      hasFitterConfig() ? configFile->contains("regularizationConfig") : false;

  if (hasRegularizationConfig) {
    auto RegularizationConfig =
        static_cast<DictNode*>(&(*configFile)[fitter]["regularizationConfig"]);
    lambda = parseDouble(*RegularizationConfig, "lambda", defaultValue, "regularizationConfig");
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
      throw sgpp::base::data_exception(errorMsg.c_str());
    }
  } else {
    std::cout << "# Did not find  " << parentDict << "[" << key << "]. Setting default value "
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
      throw sgpp::base::data_exception(errorMsg.c_str());
    }
  } else {
    std::cout << "# Did not find  " << parentDict << "[" << key << "]. Setting default value "
              << defaultValue << "." << std::endl;
    return defaultValue;
  }
}

size_t DataMiningConfigParser::parseUint(DictNode& dict, const std::string& key,
                                         size_t defaultValue, const std::string& parentDict) const {
  if (dict.contains(key)) {
    try {
      return dict[key].getUInt();
    } catch (json_exception& e) {
      std::string errorMsg = "# Failed to parse unsigned integer " + parentDict + "[" + key +
                             "] from string" + dict[key].get() + ".";
      throw sgpp::base::data_exception(errorMsg.c_str());
    }
  } else {
    std::cout << "# Did not find  " << parentDict << "[" << key << "]. Setting default value "
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
      throw sgpp::base::data_exception(errorMsg.c_str());
    }
  } else {
    std::cout << "# Did not find " << parentNode << "[" << key << "]. Setting default value "
              << defaultValue << "." << std::endl;
    return defaultValue;
  }
}

int DataMiningConfigParser::parseInt(DictNode& dict, const std::string& key, int defaultValue,
                                     const std::string& parentNode) const {
  if (dict.contains(key)) {
    try {
      return static_cast<int>(dict[key].getInt());
    } catch (json_exception& e) {
      std::string errorMsg = "# Failed to parse integer integer " + parentNode + "[" + key +
                             "] from string" + dict[key].get() + ".";
      throw sgpp::base::data_exception(errorMsg.c_str());
    }
  } else {
    std::cout << "# Did not find  " << parentNode << "[" << key << "]. Setting default value "
              << defaultValue << "." << std::endl;
    return defaultValue;
  }
}

void DataMiningConfigParser::parseSLESolverConfig(DictNode& dict, SLESolverConfiguration& config,
                                                  const SLESolverConfiguration& defaults,
                                                  const std::string& parentNode) const {
  config.eps_ = parseDouble(dict, "eps", defaults.eps_, parentNode);
  config.maxIterations_ = parseUint(dict, "maxIterations", defaults.maxIterations_, parentNode);
  config.threshold_ = parseDouble(dict, "threshold", defaults.threshold_, parentNode);

  // parse  CG type
  if (dict.contains("solverType")) {
    auto parser = SolverTypeParser();
    config.type_ = parser(dict["solverType"].get());
  }
}

} /* namespace datadriven */
} /* namespace sgpp */
