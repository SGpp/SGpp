/*
 * ConfigurationParameterParser.cpp
 *
 *  Created on: Mar 25, 2015
 *      Author: pfandedd
 */

#include <fstream>
#include <iostream>
#include <sstream>

#include "ConfigurationParameters.hpp"

namespace SGPP {
  namespace base {

    ConfigurationParameters::ConfigurationParameters() {
    }

    ConfigurationParameters::ConfigurationParameters(std::string fileName,
        std::map<std::string, std::string> defaultParameters) {

      this->readFromMap(defaultParameters);

      this->readFromFile(fileName);
    }

    std::string ConfigurationParameters::operator[](std::string key) {
      return this->parameters[key];
    }

    bool ConfigurationParameters::getAsBoolean(std::string key) {
      bool asBool;
      std::stringstream converter(parameters[key]);
      converter >> std::boolalpha;
      converter >> asBool;
      return asBool;
    }

    uint64_t ConfigurationParameters::getAsUnsigned(std::string key) {
      uint32_t asUnsigned;
      std::istringstream(parameters[key]) >> asUnsigned;
      return asUnsigned;
    }

    std::vector<std::string> ConfigurationParameters::split(const std::string& s, char delim) {
      std::stringstream ss(s);
      std::string item;
      std::vector<std::string> splitted;

      while (std::getline(ss, item, delim)) {
        splitted.push_back(item);
      }

      return splitted;
    }

    void ConfigurationParameters::readFromMap(std::map<std::string, std::string> parametersMap) {
      for (auto pair : parametersMap) {
        this->parameters[pair.first] = pair.second;
      }
    }

    void ConfigurationParameters::readFromFile(std::string fileName) {
      std::ifstream file(fileName);

      if (file.is_open()) {
        std::string line;

        while (std::getline(file, line)) {
          auto splitted = this->split(line, '=');

          if (splitted.size() == 2) {
            this->parameters[splitted[0]] = splitted[1];
          }
        }
      }

      file.close();
    }

  }
}
