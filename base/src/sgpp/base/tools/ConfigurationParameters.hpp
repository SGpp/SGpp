// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/globaldef.hpp>

#include <vector>
#include <map>
#include <string>
#include <memory>
#include <iostream>
#include <cstdint>

namespace sgpp {
namespace base {

class ConfigurationParameters {
 protected:
  std::map<std::string, std::string> parameters;

 public:
  ConfigurationParameters();

  ConfigurationParameters(std::string fileName,
                          std::map<std::string, std::string> defaultParameters =
                            std::map<std::string, std::string>());

  virtual ~ConfigurationParameters();

  void set(const std::string key, std::string value);

  std::string& get(const std::string& key);

  bool getAsBoolean(const std::string& key);

  uint64_t getAsUnsigned(const std::string& key);

  std::vector<std::string> getAsList(const std::string& key);

  void readFromMap(const std::map<std::string, std::string>& parametersMap);

  void readFromFile(std::string fileName);

  virtual std::shared_ptr<ConfigurationParameters> clone() = 0;

  std::vector<std::string> getKeys();

  void writeToFile(std::string fileName);

  void clear();

 private:
  std::vector<std::string> split(const std::string& s, char delim);
};

}  // namespace base
}  // namespace sgpp

