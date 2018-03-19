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

#include <sgpp/base/tools/json/JSON.hpp>
#include <sgpp/datadriven/datamining/modules/fitting/FitterConfiguration.hpp>
#include <sgpp/solver/TypesSolver.hpp>
#include <string>

namespace sgpp {
namespace datadriven {
struct CrossValidationConfiguration;
struct DataSourceConfig;
struct TestingConfiguration;
} /* namespace datadriven */
} /* namespace sgpp */

using json::JSON;
using json::DictNode;

using sgpp::solver::SLESolverConfiguration;
using sgpp::base::AdpativityConfiguration;
using sgpp::base::RegularGridConfiguration;
using sgpp::datadriven::CrossvalidationConfiguration;
using sgpp::datadriven::DensityEstimationConfiguration;

namespace sgpp {
namespace datadriven {

// forward declaration to break dependency cycle
enum class FitterType;

class DataMiningConfigParser {
 public:
  explicit DataMiningConfigParser(const std::string& filepath);
  virtual ~DataMiningConfigParser();

  bool hasDataSourceConfig() const;
  bool hasScorerConfig() const;
  bool hasScorerConfigCrossValidation() const;
  bool hasScorerConfigTesting() const;
  bool hasFitterConfig() const;
  bool hasScorerTestset() const;

  bool getScorerTestset(DataSourceConfig& config, const DataSourceConfig& defaults) const;
  bool getDataSourceConfig(DataSourceConfig& config, const DataSourceConfig& defaults) const;
  bool getScorerTestingConfig(TestingConfiguration& config,
                              const TestingConfiguration& defaults) const;
  bool getScorerCrossValidationConfig(CrossValidationConfiguration& config,
                                      const CrossValidationConfiguration& defaults) const;
  bool getFitterConfigType(FitterType& fitter, const FitterType& defaults) const;
  bool getFitterGridConfig(RegularGridConfiguration& config,
                           const RegularGridConfiguration& defaults) const;
  bool getFitterAdaptivityConfig(AdpativityConfiguration& config,
                                 const AdpativityConfiguration& defaults) const;
  bool getFitterCrossvalidationConfig(CrossvalidationConfiguration& config,
                                 const CrossvalidationConfiguration& defaults) const;
  bool getFitterDensityEstimationConfig(DensityEstimationConfiguration& config,
                                 const DensityEstimationConfiguration& defaults) const;
  bool getFitterSolverRefineConfig(SLESolverConfiguration& config,
                                   const SLESolverConfiguration& defaults) const;
  bool getFitterSolverFinalConfig(SLESolverConfiguration& config,
                                  const SLESolverConfiguration& defaults) const;
  bool getFitterRegularizationConfig(RegularizationConfiguration& config,
                                     const RegularizationConfiguration& defaults) const;
  bool getFitterLambda(double& lambda, double defaultValue) const;

 private:
  std::unique_ptr<JSON> configFile;

  static const std::string dataSource;
  static const std::string scorer;
  static const std::string fitter;

  std::string parseString(DictNode& dict, const std::string& key, const std::string& defaultValue,
                          const std::string& parentNode) const;
  double parseDouble(DictNode& dict, const std::string& key, double defaultValue,
                     const std::string& parentNode) const;
  size_t parseUInt(DictNode& dict, const std::string& key, size_t defaultValue,
                   const std::string& parentNode) const;
  int64_t parseInt(DictNode& dict, const std::string& key, int64_t defaultValue,
                   const std::string& parentNode) const;
  bool parseBool(DictNode& dict, const std::string& key, bool defaultValue,
                 const std::string& parentNode) const;

  void parseSLESolverConfig(DictNode& dict, SLESolverConfiguration& config,
                            const SLESolverConfiguration& defaults,
                            const std::string& parentNode) const;

  template <typename Enumeration>
  int asInteger(Enumeration const value) const {
    return static_cast<int>(value);
  }
};
} /* namespace datadriven */
} /* namespace sgpp */
