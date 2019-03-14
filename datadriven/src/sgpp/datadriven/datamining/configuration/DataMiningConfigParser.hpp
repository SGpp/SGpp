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

#pragma once

#include <sgpp/base/grid/Grid.hpp>
#include <sgpp/base/tools/json/JSON.hpp>
#include <sgpp/datadriven/datamining/modules/dataSource/DataTransformationConfig.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/FitterConfiguration.hpp>
#include <sgpp/datadriven/datamining/modules/hpo/HPOConfig.hpp>
#include <sgpp/datadriven/datamining/modules/hpo/parameters/ContinuousParameter.hpp>
#include <sgpp/datadriven/datamining/modules/hpo/parameters/DiscreteParameter.hpp>
#include <sgpp/datadriven/datamining/modules/scoring/ScorerConfig.hpp>
#include <sgpp/solver/TypesSolver.hpp>

#include <map>
#include <string>
#include <vector>

namespace sgpp {
namespace datadriven {
struct DataSourceConfig;
} /* namespace datadriven */
} /* namespace sgpp */

using json::DictNode;
using json::JSON;

using sgpp::base::AdaptivityConfiguration;
using sgpp::base::RegularGridConfiguration;
using sgpp::datadriven::CrossvalidationConfiguration;
using sgpp::datadriven::DensityEstimationConfiguration;
using sgpp::solver::SLESolverConfiguration;

namespace sgpp {
namespace datadriven {

// forward declaration to break dependency cycle
enum class FitterType;

class DataMiningConfigParser {
 public:
  explicit DataMiningConfigParser(const std::string &filepath);
  virtual ~DataMiningConfigParser();

  bool hasDataSourceConfig() const;
  bool hasDataTransformationConfig() const;
  bool hasScorerConfig() const;
  bool hasFitterConfig() const;

  void getHyperparameters(std::map<std::string, ContinuousParameter> &conpar,
                          std::map<std::string, DiscreteParameter> &dispar,
                          std::map<std::string, DiscreteParameter> &catpar,
                          std::vector<base::GridType> &basisFunctions) const;
  void getHPOConfig(HPOConfig &config);
  std::vector<int64_t> parseIntArray(DictNode &dict, const std::string &key,
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
   * Reads the configuration for the scorer
   * @param config the configuration instance to initialize
   * @param defaults a set of configurations initialized with default values
   * @return if the configuration file contained a scorer configuration
   */
  bool getScorerConfig(ScorerConfiguration &config, const ScorerConfiguration &defaults) const;

  bool getFitterConfigType(FitterType &fitter, const FitterType &defaults) const;
  bool getFitterGridConfig(RegularGridConfiguration &config,
                           const RegularGridConfiguration &defaults) const;
  bool getFitterAdaptivityConfig(AdaptivityConfiguration &config,
                                 const AdaptivityConfiguration &defaults) const;
  bool getFitterCrossvalidationConfig(CrossvalidationConfiguration &config,
                                      const CrossvalidationConfiguration &defaults) const;
  bool getFitterDensityEstimationConfig(DensityEstimationConfiguration &config,
                                        const DensityEstimationConfiguration &defaults) const;
  bool getFitterSolverRefineConfig(SLESolverConfiguration &config,
                                   const SLESolverConfiguration &defaults) const;
  bool getFitterSolverFinalConfig(SLESolverConfiguration &config,
                                  const SLESolverConfiguration &defaults) const;
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

 private:
  std::unique_ptr<JSON> configFile;

  static const std::string dataSource;
  static const std::string scorer;
  static const std::string fitter;

  std::string parseString(DictNode &dict, const std::string &key, const std::string &defaultValue,
                          const std::string &parentNode) const;
  double parseDouble(DictNode &dict, const std::string &key, double defaultValue,
                     const std::string &parentNode) const;
  std::vector<double> parseDoubleArray(DictNode &dict, const std::string &key,
                                       std::vector<double> defaultValue,
                                       const std::string &parentNode) const;
  std::vector<size_t> parseUIntArray(DictNode &dict, const std::string &key,
                                     std::vector<size_t> defaultValue,
                                     const std::string &parentNode) const;
  size_t parseUInt(DictNode &dict, const std::string &key, size_t defaultValue,
                   const std::string &parentNode) const;
  int64_t parseInt(DictNode &dict, const std::string &key, int64_t defaultValue,
                   const std::string &parentNode) const;
  bool parseBool(DictNode &dict, const std::string &key, bool defaultValue,
                 const std::string &parentNode) const;

  void parseSLESolverConfig(DictNode &dict, SLESolverConfiguration &config,
                            const SLESolverConfiguration &defaults,
                            const std::string &parentNode) const;

  void parseDataTransformationConfig(DictNode &dict, DataTransformationConfig &config,
                                     const DataTransformationConfig &defaults,
                                     const std::string &parentNode) const;
  void parseRosenblattTransformationConfig(DictNode &dict, RosenblattTransformationConfig &config,
                                           const RosenblattTransformationConfig &defaults,
                                           const std::string &parentNode) const;

  template <typename Enumeration>
  int asInteger(Enumeration const value) const {
    return static_cast<int>(value);
  }
};
} /* namespace datadriven */
} /* namespace sgpp */
