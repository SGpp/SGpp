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
  //  void getFitterConfigType() const;
  //  void getFitterGridConfig() const;
  //  void getFitterAdaptivityConfig() const;
  //  void getFitterSolverRefineConfig() const;
  //  void getFitterSolverFinalConfig() const;
  //  void getFitterRegularizationConfig() const;
  //  double getFitterLambda() const;

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
