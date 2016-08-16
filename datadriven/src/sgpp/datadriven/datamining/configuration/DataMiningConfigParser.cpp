/*
 * DataMiningParser.cpp
 *
 *  Created on: Aug 14, 2016
 *      Author: michael
 */

#include "DataMiningConfigParser.hpp"
#include <sgpp/base/exception/file_exception.hpp>
#include <sgpp/base/tools/json/json_exception.hpp>

#include <string>

using json::JSON;
using json::json_exception;
using sgpp::base::file_exception;
using json::DictNode;

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
      config.fileType = parser(*dataSourceConfig["fileType"].get());
      if (config.fileType == DataSourceFileType::NONE) {
        std::cout << "# Failed to parse dataSoure[fileType]. Setting default value "
                  << DataSourceFileType::NONE << "." << std::endl;
      }
    } else {
      std::cout << "# Did Not find dataSoure[fileType]. Setting Default value " << defaults.fileType
                << "." << std::endl;
      config.fileType = defaults.fileType;
    }
  } else {
    std::cout << "# Could not find specification  of dataSource. Falling Back to default values."
              << std::endl;
    config = defaults;
  }
  return hasDataSource;
}

std::string DataMiningConfigParser::parseString(DictNode& dict, const std::string& key,
                                                const std::string& defaultValue,
                                                const std::string& parentDict) const {
  if (dict.contains(key)) {
    try {
      return dict[key].get();
    } catch (json_exception& e) {
      std::cout << "# Failed to parse string " << parentDict << "[" << key
                << "]. Setting default value " << defaultValue << "." << std::endl;
      return defaultValue;
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
      std::cout << "# Failed to parse double " << parentDict << "[" << key
                << "]. Setting default value " << defaultValue << "." << std::endl;
      return defaultValue;
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
      std::cout << "# Failed to parse unsigned integer " << parentDict << "[" << key
                << "]. Setting default value " << defaultValue << "." << std::endl;
      return defaultValue;
    }
  } else {
    std::cout << "# Did not find  " << parentDict << "[" << key << "]. Setting default value "
              << defaultValue << "." << std::endl;
    return defaultValue;
  }
}

bool DataMiningConfigParser::parseBool(DictNode& dict, const std::string& key, bool defaultValue,
                                       const std::string& parentDict) const {
  if (dict.contains(key)) {
    try {
      return dict[key].getBool();
    } catch (json_exception& e) {
      std::cout << "# Failed to parse bool " << parentDict << "[" << key
                << "]. Setting default value " << defaultValue << "." << std::endl;
      return defaultValue;
    }
  } else {
    std::cout << "# Did not find " << parentDict << "[" << key << "]. Setting default value "
              << defaultValue << "." << std::endl;
    return defaultValue;
  }
}

} /* namespace datadriven */
} /* namespace sgpp */
