// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/tools/ConfigurationParameters.hpp>
#include <sgpp/base/exception/operation_exception.hpp>

#include <vector>
#include <utility>
#include <map>
#include <string>
#include <fstream>
#include <iostream>
#include <sstream>
#include <memory>

namespace sgpp {
namespace base {

ConfigurationParameters::ConfigurationParameters() {
}

ConfigurationParameters::ConfigurationParameters(std::string fileName,
    std::map<std::string, std::string> defaultParameters) {

  this->readFromMap(defaultParameters);

  this->readFromFile(fileName);
}

ConfigurationParameters::~ConfigurationParameters() {
}

void ConfigurationParameters::set(std::string key, std::string value) {
  this->parameters[key] = value;
}

std::string& ConfigurationParameters::get(const std::string& key) {
  if (this->parameters.count(key) != 1) {
    std::stringstream errorString;
    errorString << "OCL Error: parameter \"" << key << "\" used but not set" <<
                std::endl;
    throw sgpp::base::operation_exception(errorString.str());
  }

  return this->parameters[key];
}

bool ConfigurationParameters::getAsBoolean(const std::string& key) {
  if (this->parameters.count(key) != 1) {
    std::stringstream errorString;
    errorString << "OCL Error: parameter \"" << key << "\" used but not set" <<
                std::endl;
    throw sgpp::base::operation_exception(errorString.str());
  }

  bool asBool;
  std::stringstream converter(parameters[key]);
  converter >> std::boolalpha;
  converter >> asBool;
  return asBool;
}

uint64_t ConfigurationParameters::getAsUnsigned(const std::string& key) {
  if (this->parameters.count(key) != 1) {
    std::stringstream errorString;
    errorString << "OCL Error: parameter \"" << key << "\" used but not set" <<
                std::endl;
    throw sgpp::base::operation_exception(errorString.str());
  }

  uint32_t asUnsigned;
  std::istringstream(parameters[key]) >> asUnsigned;
  return asUnsigned;
}

std::vector<std::string> ConfigurationParameters::split(const std::string& s,
    char delim) {
  std::stringstream ss(s);
  std::string item;
  std::vector<std::string> splitted;

  while (std::getline(ss, item, delim)) {
    splitted.push_back(item);
  }

  return splitted;
}

void ConfigurationParameters::readFromMap(
  const std::map<std::string, std::string>& parametersMap) {
  for (auto pair : parametersMap) {
    this->parameters[pair.first] = pair.second;
  }
}

void ConfigurationParameters::readFromFile(std::string fileName) {
  std::ifstream file(fileName);

  if (file.is_open()) {
    std::string line;

    while (std::getline(file, line)) {
      if (line[0] == '#') {
        continue;
      }

      auto splitted = this->split(line, '=');

      if (splitted.size() == 2) {
        this->parameters[splitted[0]] = splitted[1];
      }
    }
  } else {
    throw;
  }

  file.close();
}

std::vector<std::string> ConfigurationParameters::getAsList(
  const std::string& key) {
  std::string& listAttribute = this->get(key);

  auto splitted = this->split(listAttribute, ',');

  return splitted;
}

std::vector<std::string> ConfigurationParameters::getKeys() {
  std::vector<std::string> keys;

  for (auto& tuple : this->parameters) {
    keys.push_back(tuple.first);
  }

  return keys;
}

void ConfigurationParameters::writeToFile(std::string fileName) {
  std::ofstream file(fileName);

  if (file.is_open()) {
    for (std::pair<std::string, std::string> tuple : this->parameters) {
      file << tuple.first << "=" << tuple.second << std::endl;
    }

  } else {
    throw;
  }

  file.close();
}

void ConfigurationParameters::clear() {
  this->parameters.clear();
}

}  // namespace base
}  // namespace sgpp
