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

#include <sgpp/datadriven/datamining/configuration/DataMiningConfigParser.hpp>

#include <sgpp/base/exception/data_exception.hpp>
#include <sgpp/base/exception/file_exception.hpp>
#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/grid/LevelIndexTypes.hpp>
#include <sgpp/base/tools/json/JSON.hpp>
#include <sgpp/base/tools/json/ListNode.hpp>
#include <sgpp/base/tools/json/json_exception.hpp>
#include <sgpp/datadriven/configuration/RegularizationConfiguration.hpp>
#include <sgpp/datadriven/datamining/configuration/DensityEstimationTypeParser.hpp>
#include <sgpp/datadriven/datamining/configuration/GeneralGridTypeParser.hpp>
#include <sgpp/datadriven/datamining/configuration/GridTypeParser.hpp>
#include <sgpp/datadriven/configuration/GeometryConfiguration.hpp>
#include <sgpp/datadriven/datamining/configuration/GeometryConfigurationParser.hpp>
#include <sgpp/datadriven/datamining/configuration/MatrixDecompositionTypeParser.hpp>
#include <sgpp/datadriven/datamining/configuration/RefinementFunctorTypeParser.hpp>
#include <sgpp/datadriven/datamining/configuration/RegularizationTypeParser.hpp>
#include <sgpp/datadriven/datamining/configuration/SLESolverTypeParser.hpp>
#include <sgpp/datadriven/datamining/modules/dataSource/DataSourceFileTypeParser.hpp>
#include <sgpp/datadriven/datamining/modules/dataSource/DataTransformationTypeParser.hpp>
#include <sgpp/datadriven/datamining/modules/dataSource/shuffling/DataSourceShufflingTypeParser.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/FitterTypeParser.hpp>
#include <sgpp/datadriven/datamining/modules/scoring/ScorerMetricTypeParser.hpp>
#include <sgpp/solver/TypesSolver.hpp>
#include <sgpp/datadriven/datamining/modules/visualization/VisualizationTypesParser.hpp>

#include <map>
#include <string>
#include <vector>

using json::DictNode;
using json::JSON;
using json::json_exception;
using sgpp::base::data_exception;
using sgpp::base::file_exception;
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
       hasVisualizationConfig() ? (*configFile)[visualization]
           .contains("generalConfig") : false;
  return hasVisualizationParameters;
}

bool DataMiningConfigParser::hasVisualizationParametersConfig() const {
  bool hasVisualizationParameters =
       hasVisualizationConfig() ? (*configFile)[visualization]
           .contains("parameters") : false;
  return hasVisualizationParameters;
}

bool DataMiningConfigParser::getDataSourceConfig(DataSourceConfig &config,
                                                 const DataSourceConfig &defaults) const {
  bool hasDataSource = hasDataSourceConfig();

  if (hasDataSource) {
    auto dataSourceConfig = static_cast<DictNode *>(&(*configFile)[dataSource]);

    config.filePath = parseString(*dataSourceConfig, "filePath", defaults.filePath, "dataSource");
    config.isCompressed =
        parseBool(*dataSourceConfig, "compression", defaults.isCompressed, "dataSource");
    config.numBatches =
        parseUInt(*dataSourceConfig, "numBatches", defaults.numBatches, "dataSource");
    config.batchSize = parseUInt(*dataSourceConfig, "batchSize", defaults.batchSize, "dataSource");
    config.hasTargets =
        parseBool(*dataSourceConfig, "hasTargets", defaults.hasTargets, "dataSource");
    config.validationPortion = parseDouble(*dataSourceConfig, "validationPortion",
                                           defaults.validationPortion, "dataSource");
    // if negative we want UINT_MAX here, so all should be fine
    config.readinCutoff = static_cast<size_t>(
        parseInt(*dataSourceConfig, "readinCutoff", defaults.readinCutoff, "dataSource"));
    config.readinClasses =
        parseDoubleArray(*dataSourceConfig, "readinClasses", defaults.readinClasses, "dataSource");
    config.readinColumns =
        parseUIntArray(*dataSourceConfig, "readinColumns", defaults.readinColumns, "dataSource");

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
          static_cast<DictNode *>(&(*configFile)[dataSource]["dataTransformation"]);
      parseDataTransformationConfig(*dataTransformationConfig, config.dataTransformationConfig,
                                    defaults.dataTransformationConfig, "dataTransformation");
    } else {
      std::cout << "# Could not find specification of dataSource[dataTransformationConfig]. "
                   "Falling back to default values."
                << std::endl;
      config.dataTransformationConfig = defaults.dataTransformationConfig;
    }

    // parse theshuffling
    if (dataSourceConfig->contains("shuffling")) {
      config.shuffling =
          DataSourceShufflingTypeParser::parse((*dataSourceConfig)["shuffling"].get());
    } else {
      std::cout << "# Did not find dataSource[shuffling]. Setting default value "
                << DataSourceShufflingTypeParser::toString(defaults.shuffling) << "." << std::endl;
      config.shuffling = defaults.shuffling;
    }

    config.randomSeed =
        parseUInt(*dataSourceConfig, "randomSeed", defaults.randomSeed, "dataSource");
    config.epochs = parseUInt(*dataSourceConfig, "epochs", defaults.epochs, "dataSource");

    // Parse info for test data
    config.testFilePath = parseString(*dataSourceConfig, "testFilePath",
           defaults.filePath, "dataSource");

    // parse file type of test data
    if (dataSourceConfig->contains("testFileType")) {
      config.testFileType = DataSourceFileTypeParser::parse(
      (*dataSourceConfig)["testFileType"].get());
    } else {
      std::cout << "# Did not find " << dataSource <<
        "[testFileType]. Setting default value "
        << DataSourceFileTypeParser::toString(defaults.testFileType)
        << "." << std::endl;

        config.testFileType = defaults.testFileType;
    }

    config.testIsCompressed = parseBool(*dataSourceConfig, "testCompression",
            defaults.testIsCompressed, "dataSource");
    config.testNumBatches = parseUInt(*dataSourceConfig, "testNumBatches",
            defaults.testNumBatches, "dataSource");
    config.testBatchSize = parseUInt(*dataSourceConfig, "testBatchSize",
            defaults.testBatchSize, "dataSource");
    config.testHasTargets = parseBool(*dataSourceConfig, "testHasTargets",
            defaults.testHasTargets, "dataSource");

    config.testReadinCutoff = static_cast<size_t>(parseInt(*dataSourceConfig,
            "testReadinCutoff", defaults.testReadinCutoff, "dataSource"));
    config.testReadinClasses = parseDoubleArray(*dataSourceConfig,
            "testReadinClasses", defaults.testReadinClasses, "dataSource");
    config.testReadinColumns = parseUIntArray(*dataSourceConfig,
            "testReadinColumns", defaults.testReadinColumns, "dataSource");

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
      config.metric = ScorerMetricTypeParser::parse((*scorerConfig)["metric"].get());
    } else {
      std::cout << "# Did not find scorer[metric]. Setting default value "
                << ScorerMetricTypeParser::toString(defaults.metric) << "." << std::endl;
    }
  } else {
    std::cout << "# Could not find specification  of scorer. Falling Back to default values."
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
      std::cout << "# Could not find specification  of fitter[type]. Falling Back to default "
                   "values."
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
            GeneralGridTypeParser::parse((*fitterConfig)["generalGridType"].get());
      } else {
        config.generalType_ =
            GeneralGridTypeParser::parse((*fitterConfig)["generalGridType"]["value"].get());
      }
    } else {
      std::cout << "# Did not find gridConfig[generalGridType]. Setting default value."
                << std::endl;
    }

    // parse  grid type
    if (fitterConfig->contains("gridType")) {
      if ((*fitterConfig)["gridType"].size() == 1) {
        config.type_ = GridTypeParser::parse((*fitterConfig)["gridType"].get());
      } else {
        config.type_ = GridTypeParser::parse((*fitterConfig)["gridType"]["value"].get());
      }
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
    AdaptivityConfiguration &config, const AdaptivityConfiguration &defaults) const {
  bool hasFitterAdaptivityConfig =
      hasFitterConfig() ? (*configFile)[fitter].contains("adaptivityConfig") : false;

  if (hasFitterAdaptivityConfig) {
    auto adaptivityConfig = static_cast<DictNode *>(&(*configFile)[fitter]["adaptivityConfig"]);
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
    config.errorConvergenceThreshold =
        parseDouble(*adaptivityConfig, "errorConvergenceThreshold",
                    defaults.errorConvergenceThreshold, "adaptivityConfig");
    config.errorBufferSize = parseUInt(*adaptivityConfig, "errorBufferSize",
                                       defaults.errorBufferSize, "adaptivityConfig");
    config.errorMinInterval = parseUInt(*adaptivityConfig, "errorMinInterval",
                                        defaults.errorMinInterval, "adaptivityConfig");
    config.refinementPeriod = parseUInt(*adaptivityConfig, "refinementPeriod",
                                        defaults.refinementPeriod, "adaptivityConfig");
    config.precomputeEvaluations = parseBool(*adaptivityConfig, "precomputeEvaluations",
                                             defaults.precomputeEvaluations, "adaptivityConfig");
    config.levelPenalize =
        parseBool(*adaptivityConfig, "penalizeLevels", defaults.levelPenalize, "adaptivityConfig");

    // Parse scaling coefficients if present
    if (adaptivityConfig->contains("scalingCoefficients")) {
      json::ListNode &coefs =
          dynamic_cast<json::ListNode &>((*adaptivityConfig)["scalingCoefficients"]);
      for (size_t i = 0; i < coefs.size(); i++) {
        config.scalingCoefficients.push_back(coefs[i].getDouble());
      }
    }

    // Parse refinement indicator
    if (adaptivityConfig->contains("refinementIndicator")) {
      config.refinementFunctorType =
          RefinementFunctorTypeParser::parse((*adaptivityConfig)["refinementIndicator"].get());
    } else {
      std::cout << "# Did not find adaptivityConfig[refinementIndicator]. Setting default "
                << "value " << RefinementFunctorTypeParser::toString(defaults.refinementFunctorType)
                << "." << std::endl;
      config.refinementFunctorType = defaults.refinementFunctorType;
    }
  } else {
    std::cout << "# Could not find specification  of fitter[adaptivityConfig]. Falling Back to "
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
    std::cout << "# Could not find specification  of fitter[crossvalidationConfig]. Falling "
                 "Back to default values."
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

    // parse  density estimation type
    if (densityEstimationConfig->contains("densityEstimationType")) {
      config.type_ = DensityEstimationTypeParser::parse(
          (*densityEstimationConfig)["densityEstimationType"].get());
    } else {
      std::cout << "# Did not find densityEstimationConfig[densityEstimationType]."
                   " Setting default value "
                << DensityEstimationTypeParser::toString(defaults.type_) << "." << std::endl;
      config.type_ = defaults.type_;
    }

    // parse matrix decomposition type
    if (densityEstimationConfig->contains("matrixDecompositionType")) {
      config.decomposition_ = MatrixDecompositionTypeParser::parse(
          (*densityEstimationConfig)["matrixDecompositionType"].get());
    } else {
      std::cout << "# Did not find densityEstimationConfig[matrixDecompositionType]."
                   " Setting default value "
                << MatrixDecompositionTypeParser::toString(defaults.decomposition_) << "."
                << std::endl;
      config.decomposition_ = defaults.decomposition_;
    }
  } else {
    std::cout
        << "# Could not find specification  of fitter[densityEstimationConfig]. Falling Back to "
           "default values."
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
    std::cout << "# Could not find specification  of fitter[solverRefineConfig]. Falling Back to "
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
    std::cout << "# Could not find specification  of fitter[solverFinalConfig]. Falling Back to "
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
  }

  return hasRegularizationConfig;
}

bool DataMiningConfigParser::getVisualizationGeneralConfig(
  VisualizationGeneralConfig &config, const VisualizationGeneralConfig &defaults) const {
  bool hasGeneralConfig = hasVisualizationGeneralConfig();

  if (hasGeneralConfig) {
    auto visualizationGeneralConfig =
     static_cast<DictNode *>(&(*configFile)[visualization]["generalConfig"]);

    std::cout << "Starting reading visualization " << std::endl;
    config.algorithm = parseString(*visualizationGeneralConfig, "algorithm",
     defaults.algorithm, "visualization");

    config.targetDirectory = parseString(*visualizationGeneralConfig, "targetDirectory",
     defaults.targetDirectory, "visualization");

    // parse file type
    if (visualizationGeneralConfig->contains("targetFileType")) {
      config.targetFileType = VisualizationTypesParser::
        parseFileType((*visualizationGeneralConfig)["targetFileType"].get());
    } else {
      std::cout << "# Did not find " << dataSource << "[fileType]. Setting default value "
                << VisualizationTypesParser::toString(defaults.targetFileType) << "." << std::endl;
      config.targetFileType = defaults.targetFileType;
    }

    // Checks if the specified number of batches for visualization exceeds the one
    // of the data source. If that's the case set the visualization number of batches
    // to the one of the data source.
    auto dataSourceConfig = static_cast<DictNode *>(&(*configFile)[dataSource]);

    config.numBatches = parseUInt(*visualizationGeneralConfig, "numBatches",
      defaults.numBatches, "visualization");


  } else {
    std::cout << "# Could not find specification of "
     "visualization general config. Falling Back to default values." << std::endl;
    config = defaults;
  }

  // This is out of the if, since the default value is not directly given by the user
  bool hasFitterCrossvalidationConfig =
     hasFitterConfig() ? (*configFile)[fitter].contains("crossValidation") : false;
  if (hasFitterCrossvalidationConfig) {
     auto crossvalidationConfig = static_cast<DictNode *>(&(*configFile)[fitter]
                                                                         ["crossValidation"]);
     config.crossValidation = parseBool(*crossvalidationConfig, "enable",
       defaults.crossValidation, "crossValidation");
   }
   else {
     config.crossValidation = false;
  }
  return hasGeneralConfig;
}

bool DataMiningConfigParser::getVisualizationParameters(
  VisualizationParameters &config, const VisualizationParameters &defaults) const {
  bool hasVisualizationParameters = hasVisualizationParametersConfig();

  if (hasVisualizationParameters) {
    auto visualizationParameters =
      static_cast<DictNode *>(&(*configFile)[visualization]["parameters"]);

    config.perplexity = parseDouble(*visualizationParameters, "perplexity",
      defaults.perplexity, "visualization");

    config.theta = parseDouble(*visualizationParameters, "theta",
      defaults.theta, "visualization");

    config.seed = parseUInt(*visualizationParameters, "seed",
      defaults.seed, "visualization");

    config.maxNumberIterations = parseUInt(*visualizationParameters,
      "maxNumberIterations", defaults.maxNumberIterations, "visualization");

    config.targetDimension = parseUInt(*visualizationParameters, "targetDimension",
      defaults.targetDimension, "visualization");

    config.numberCores = parseUInt(*visualizationParameters,
      "numberCores", defaults.numberCores, "visualization");
  } else {
    std::cout << "# Could not find specification of visualization parameters. "
      "Falling Back to default values." << std::endl;
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
    } catch (json_exception &e) {
      try {
        return dict[key]["value"].get();
      } catch (json_exception &e) {
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
    } catch (json_exception &e) {
      try {
        if (dict[key].contains("logscale")) {
          if (dict[key]["logscale"].getBool()) {
            return pow(10, dict[key]["value"].getDouble());
          }
        }
        return dict[key]["value"].getDouble();
      } catch (json_exception &e) {
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
    } catch (json_exception &e) {
      try {
        return dict[key]["value"].getUInt();
      } catch (json_exception &e) {
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
    } catch (json_exception &e) {
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
    } catch (json_exception &e) {
      try {
        return dict[key]["value"].getInt();
      } catch (json_exception &e) {
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
    } catch (json_exception &e) {
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
    } catch (json_exception &e) {
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
    } catch (json_exception &e) {
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
    config.type_ = SLESolverTypeParser::parse(dict["solverType"].get());
  } else {
    std::cout << "# Did not find " << parentNode << "[solverType]. Setting default value "
              << SLESolverTypeParser::toString(defaults.type_) << "." << std::endl;
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
  } catch (json_exception &e) {
  }
  try {
    if ((*configFile)[fitter]["gridConfig"]["gridType"]["optimize"].getBool()) {
      size_t nOptions = (*configFile)[fitter]["gridConfig"]["gridType"]["options"].size();
      if (nOptions > 1) {
        for (size_t i = 0; i < nOptions; ++i) {
          basisFunctions.push_back(GridTypeParser::parse(
              (*configFile)[fitter]["gridConfig"]["gridType"]["options"][i].get()));
        }
        catpar["basisFunction"] =
            DiscreteParameter("basisFunction", 0, static_cast<int>(nOptions - 1));
      }
    }
  } catch (json_exception &e) {
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
  } catch (json_exception &e) {
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
  } catch (json_exception &e) {
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
  } catch (json_exception &e) {
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
      std::cout << "# Could not find specification  of hpo[harmonica]. Falling Back to "
                   "default values."
                << std::endl;
    }
    if (node->contains("bayesianOptimization")) {
      auto bo = static_cast<DictNode *>(&(*node)["bayesianOptimization"]);
      config.setNRandom(parseInt(*bo, "nRandom", config.getNRandom(), "hpo"));
      config.setNRuns(parseInt(*bo, "nRuns", config.getNRuns(), "hpo"));
    } else {
      std::cout << "# Could not find specification  of hpo[bayesianOptimization]. Falling Back to "
                   "default values."
                << std::endl;
    }
  } else {
    std::cout << "# Could not find specification  of hpo. Falling Back to "
                 "default values."
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
    config.type = DataTransformationTypeParser::parse(dict["type"].get());
  } else {
    std::cout << "# Did not find [dataTransformationType]. Setting default value "
              << DataTransformationTypeParser::toString(defaults.type) << "." << std::endl;
    config.type = defaults.type;
  }

  // If type Rosenblatt parse RosenblattTransformationConfig
  if (config.type == DataTransformationType::ROSENBLATT) {
    auto rosenblattTransformationConfig = static_cast<DictNode *>(
        &(*configFile)[dataSource]["dataTransformation"]["rosenblattConfig"]);
    parseRosenblattTransformationConfig(*rosenblattTransformationConfig, config.rosenblattConfig,
                                        defaults.rosenblattConfig, "rosenblattConfig");
  } else {
    std::cout << "# Could not find specification of dataSource[dataTransformationConfig]"
                 "[rosenblattConfig]. Falling back to default values."
              << std::endl;
    config.rosenblattConfig = defaults.rosenblattConfig;
  }
}

void DataMiningConfigParser::parseRosenblattTransformationConfig(
    DictNode &dict, RosenblattTransformationConfig &config,
    const RosenblattTransformationConfig &defaults, const std::string &parentNode) const {
  config.numSamples = parseUInt(dict, "numSamples", defaults.numSamples, parentNode);
  config.gridLevel = parseUInt(dict, "gridLevel", defaults.gridLevel, parentNode);

  config.solverMaxIterations =
      parseUInt(dict, "solverMaxIterations", defaults.solverMaxIterations, parentNode);
  config.solverEps = parseDouble(dict, "solverEps", defaults.solverEps, parentNode);
  config.solverThreshold =
      parseDouble(dict, "solverThreshold", defaults.solverThreshold, parentNode);
}


bool DataMiningConfigParser::getFitterDatabaseConfig(
    datadriven::DatabaseConfiguration &config,
    const datadriven::DatabaseConfiguration &defaults) const {
  bool hasDatabaseConfig = hasFitterConfig() ? (*configFile)[fitter].contains("database") : false;

  if (hasDatabaseConfig) {
    auto databaseConfig = static_cast<DictNode *>(&(*configFile)[fitter]["database"]);

    // Parse filepath
    if (databaseConfig->contains("filepath")) {
      config.filepath = (*databaseConfig)["filepath"].get();
    } else {
      std::cout << "# Did not find databaseConfig[filepath]. No database loaded" << std::endl;
      config.filepath = defaults.filepath;
    }
  }

  return hasDatabaseConfig;
}

bool DataMiningConfigParser::getFitterLearnerConfig(
    datadriven::LearnerConfiguration &config,
    const datadriven::LearnerConfiguration &defaults) const {
  bool hasLearnerConfig = hasFitterConfig() ? (*configFile)[fitter].contains("learner") : false;

  if (hasLearnerConfig) {
    auto learnerConfig = static_cast<DictNode *>(&(*configFile)[fitter]["learner"]);

    config.beta = parseDouble(*learnerConfig, "beta", defaults.beta, "learnerConfig");
    config.usePrior = parseBool(*learnerConfig, "usePrior", defaults.usePrior, "learnerConfig");
  }

  return hasLearnerConfig;
}

bool DataMiningConfigParser::getGeometryConfig(
    datadriven::GeometryConfiguration &config,
    const datadriven::GeometryConfiguration &defaults) const {
  bool hasGeometryConfig = hasFitterConfig() ?
    (*configFile)[fitter].contains("geometryConfig") : false;

  if (hasGeometryConfig) {
    std::cout << "Has geometry config" << std::endl;
    auto geometryConfig = static_cast<DictNode *>(&(*configFile)[fitter]["geometryConfig"]);

    config.dim = parseIntArray(*geometryConfig, "dim", defaults.dim, "geometryConfig");

    // parse  density estimation type
    if (geometryConfig->contains("stencil")) {
      config.stencilType = GeometryConfigurationParser::parse(
          (*geometryConfig)["stencil"].get());
    } else {
      std::cout << "# Did not find geometryConfig[stencil]."
                   " Setting default value ";
      config.stencilType = defaults.stencilType;
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
