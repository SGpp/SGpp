/*
 * DataMiningParser.hpp
 *
 *  Created on: Aug 14, 2016
 *      Author: Michael Lettrich
 */

#pragma once

#include <sgpp/datadriven/datamining/modules/dataSource/DataSourceConfig.hpp>

#include <sgpp/base/tools/json/JSON.hpp>
#include <string>

using json::JSON;
using json::DictNode;

using sgpp::solver::SLESolverConfiguration;
using sgpp::base::AdpativityConfiguration;
using sgpp::base::RegularGridConfiguration;

namespace sgpp {
namespace datadriven {

class DataMiningConfigParser {
 public:
  DataMiningConfigParser(const std::string& filepath);
  virtual ~DataMiningConfigParser();

  bool hasDataSourceConfig() const;
  bool hasScorerConfig() const;
  bool hasFitterConfig() const;

  bool getDataSourceConfig(DataSourceConfig& config, const DataSourceConfig& defaults) const;
  //  void getScorerTestingConfig() const;
  //  void getScorerCrossValidationConfig() const;
  // bool getFitterConfigType() const;
  bool getFitterGridConfig(RegularGridConfiguration& config,
                           const RegularGridConfiguration& defaults) const;
  bool getFitterAdaptivityConfig(AdpativityConfiguration& config,
                                 const AdpativityConfiguration& defaults) const;
  bool getFitterSolverRefineConfig(SLESolverConfiguration& config,
                                   const SLESolverConfiguration& defaults) const;
  bool getFitterSolverFinalConfig(SLESolverConfiguration& config,
                                  const SLESolverConfiguration& defaults) const;
  bool getFitterRegularizationConfig(RegularizationConfiguration& config,
                                     const RegularizationConfiguration& defaults) const;
  bool getFitterLambda(double& lambda, double defaultValue) const;

 private:
  std::shared_ptr<JSON> configFile;

  static const std::string dataSource;
  static const std::string scorer;
  static const std::string fitter;

  std::string parseString(DictNode& dict, const std::string& key, const std::string& defaultValue,
                          const std::string& parentNode) const;
  double parseDouble(DictNode& dict, const std::string& key, double defaultValue,
                     const std::string& parentNode) const;
  size_t parseUint(DictNode& dict, const std::string& key, size_t defaultValue,
                   const std::string& parentNode) const;
  bool parseBool(DictNode& dict, const std::string& key, bool defaultValue,
                 const std::string& parentNode) const;
};
} /* namespace datadriven */
} /* namespace sgpp */
