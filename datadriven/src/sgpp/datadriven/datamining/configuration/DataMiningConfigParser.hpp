// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/tools/json/JSON.hpp>
#include <sgpp/datadriven/configuration/GeometryConfiguration.hpp>
#include <sgpp/datadriven/datamining/modules/dataSource/DataTransformationConfig.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/FitterConfiguration.hpp>
#include <sgpp/datadriven/datamining/modules/hpo/HPOConfig.hpp>
#include <sgpp/datadriven/datamining/modules/hpo/parameters/ContinuousParameter.hpp>
#include <sgpp/datadriven/datamining/modules/hpo/parameters/DiscreteParameter.hpp>
#include <sgpp/datadriven/datamining/modules/scoring/ScorerConfig.hpp>
#include <sgpp/datadriven/datamining/modules/visualization/VisualizationGeneralConfig.hpp>
#include <sgpp/datadriven/datamining/modules/visualization/VisualizationParameters.hpp>
#include <sgpp/solver/TypesSolver.hpp>

#include <map>
#include <string>
#include <vector>

namespace sgpp {
namespace datadriven {

// forward declarations to break dependency cycles
struct DataSourceConfig;
enum class FitterType;

class DataMiningConfigParser {
 public:
  explicit DataMiningConfigParser(const std::string &filepath);
  virtual ~DataMiningConfigParser();

  bool hasDataSourceConfig() const;
  bool hasDataTransformationConfig() const;
  bool hasScorerConfig() const;
  bool hasFitterConfig() const;
  bool hasGeometryConfig() const;
  bool hasParallelConfig() const;
  bool hasVisualizationConfig() const;
  bool hasVisualizationGeneralConfig() const;
  bool hasVisualizationParametersConfig() const;

  void getHyperparameters(std::map<std::string, ContinuousParameter> &conpar,
                          std::map<std::string, DiscreteParameter> &dispar,
                          std::map<std::string, DiscreteParameter> &catpar,
                          std::vector<base::GridType> &basisFunctions) const;
  void getHPOConfig(HPOConfig &config);
  std::vector<int64_t> parseIntArray(json::JSON::DictNode &dict, const std::string &key,
                                     std::vector<int64_t> defaultValue,
                                     const std::string &parentNode) const;
  std::string getHPOMethod(std::string defaultValue) const;

  /**
   * Checks whether the fitter configuration contains a cross validation configuration
   * @return if the fitter configuration contains a cross validation configuration
   */
  bool hasFitterConfigCrossValidation() const;

  bool getDataSourceConfig(DataSourceConfig &config, const DataSourceConfig &defaults) const;

  /**
   * Expands on getDataSourceConfig by allowing multiple dataSource instances. The first instance is
   * used as a template for filling in the parameters for all dataSource instances.
   */
  bool getMultiDataSourceConfig(std::vector<DataSourceConfig> &config,
                                const std::vector<DataSourceConfig> &defaults) const;

  /**
   * Reads the configuration for the scorer
   * @param config the configuration instance to initialize
   * @param defaults a set of configurations initialized with default values
   * @return if the configuration file contained a scorer configuration
   */
  bool getScorerConfig(ScorerConfiguration &config, const ScorerConfiguration &defaults) const;

  bool getFitterConfigType(FitterType &fitter, const FitterType &defaults) const;

  bool getFitterGridConfig(base::GeneralGridConfiguration &config,
                           const base::GeneralGridConfiguration &defaults) const;

  bool getFitterAdaptivityConfig(base::AdaptivityConfiguration &config,
                                 const base::AdaptivityConfiguration &defaults) const;

  bool getFitterCrossvalidationConfig(CrossvalidationConfiguration &config,
                                      const CrossvalidationConfiguration &defaults) const;
  bool getFitterDensityEstimationConfig(DensityEstimationConfiguration &config,
                                        const DensityEstimationConfiguration &defaults) const;
  bool getFitterSolverRefineConfig(solver::SLESolverConfiguration &config,
                                   const solver::SLESolverConfiguration &defaults) const;
  bool getFitterSolverFinalConfig(solver::SLESolverConfiguration &config,
                                  const solver::SLESolverConfiguration &defaults) const;
  bool getFitterRegularizationConfig(RegularizationConfiguration &config,
                                     const RegularizationConfiguration &defaults) const;
  bool getFitterLambda(double &lambda, double defaultValue) const;

  /**
   * Returns the database configuration of the fitter if it exists
   * @param config the configuration object that will be initialized
   * @param defaults default values if the fitter does not contain a database configuration
   * @return whether the fitter contains a database configuration
   */
  bool getFitterDatabaseConfig(datadriven::DatabaseConfiguration &config,
                               const datadriven::DatabaseConfiguration &defaults) const;

  /**
   * Initializes the learner configuration if it exists
   * @param config the configuration instance that will be initialized
   * @param defaults default values if the fitter config does not contain a matching entry
   * @return whether the configuration contains a learner configuration
   */
  bool getFitterLearnerConfig(datadriven::LearnerConfiguration &config,
                              const datadriven::LearnerConfiguration &defaults) const;

  /**
   * Initializes the parallel configuration if it exists
   * @param config the configuration instance that will be initialized
   * @param defaults default values if the parallel config does not contain a matching entry
   * @return whether the configuration contains a parallel configuration
   */
  bool getFitterParallelConfig(datadriven::ParallelConfiguration &config,
                               const datadriven::ParallelConfiguration &defaults) const;

  /**
   * Initializes the geometry configuration if it exists
   * @param config the configuration instance that will be initialized
   * @param defaults default values if the fitter config does not contain a matching entry
   * @return whether the configuration contains a learner configuration
   */
  bool getGeometryConfig(datadriven::GeometryConfiguration &config,
                         const datadriven::GeometryConfiguration &defaults) const;

  /**
   * Initializes the visualization general configuration if it exists
   * @param config the configuration instance that will be initiazlized
   * @param defaults default values if the fitter config does not contain a matching entry
   * @return wether the configuration contains a visualization configuration
   */
  bool getVisualizationGeneralConfig(datadriven::VisualizationGeneralConfig &config,
                                     const datadriven::VisualizationGeneralConfig &defaults) const;

  /**
   * Initializes the visualization parameters configuration if it exists
   * @param config the configuration instance that will be initiazlized
   * @param defaults default values if the fitter config does not contain a matching entry
   * @return wether the configuration contains a visualization configuration
   */
  bool getVisualizationParameters(datadriven::VisualizationParameters &config,
                                  const datadriven::VisualizationParameters &defaults) const;

 private:
  std::unique_ptr<json::JSON> configFile;

  static const std::string dataSource;
  static const std::string scorer;
  static const std::string fitter;
  static const std::string visualization;

  std::string parseString(json::JSON::DictNode &dict, const std::string &key,
                          const std::string &defaultValue, const std::string &parentNode) const;
  double parseDouble(json::JSON::DictNode &dict, const std::string &key, double defaultValue,
                     const std::string &parentNode) const;
  std::vector<double> parseDoubleArray(json::JSON::DictNode &dict, const std::string &key,
                                       std::vector<double> defaultValue,
                                       const std::string &parentNode) const;

  std::vector<std::vector<int64_t>> parseArrayOfIntArrays(
      json::JSON::DictNode &dict, const std::string &key,
      std::vector<std::vector<int64_t>> defaultValue, const std::string &parentNode) const;
  std::vector<std::string> parseStringArray(json::JSON::DictNode &dict, const std::string &key,
                                            std::vector<std::string> defaultValue,
                                            const std::string &parentNode) const;
  std::vector<size_t> parseUIntArray(json::JSON::DictNode &dict, const std::string &key,
                                     std::vector<size_t> defaultValue,
                                     const std::string &parentNode) const;
  size_t parseUInt(json::JSON::DictNode &dict, const std::string &key, size_t defaultValue,
                   const std::string &parentNode) const;
  int64_t parseInt(json::JSON::DictNode &dict, const std::string &key, int64_t defaultValue,
                   const std::string &parentNode) const;
  bool parseBool(json::JSON::DictNode &dict, const std::string &key, bool defaultValue,
                 const std::string &parentNode) const;

  void parseSLESolverConfig(json::JSON::DictNode &dict, solver::SLESolverConfiguration &config,
                            const solver::SLESolverConfiguration &defaults,
                            const std::string &parentNode) const;

  void parseDataTransformationConfig(json::JSON::DictNode &dict, DataTransformationConfig &config,
                                     const DataTransformationConfig &defaults,
                                     const std::string &parentNode) const;

  void parseRosenblattTransformationConfig(json::JSON::DictNode &dict,
                                           RosenblattTransformationConfig &config,
                                           const RosenblattTransformationConfig &defaults,
                                           const std::string &parentNode) const;

  template <typename Enumeration>
  int asInteger(Enumeration const value) const {
    return static_cast<int>(value);
  }
};
} /* namespace datadriven */
} /* namespace sgpp */
