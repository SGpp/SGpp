/*
 * ConfigurationParameterParser.cpp
 *
 *  Created on: Mar 25, 2015
 *      Author: pfandedd
 */

#include <fstream>
#include <iostream>
#include <sstream>
#include <memory>

#include "ConfigurationParameters.hpp"

#include <sgpp/base/exception/operation_exception.hpp>

namespace SGPP {
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

std::string &ConfigurationParameters::get(const std::string &key) {
    if (this->parameters.count(key) != 1) {
        std::stringstream errorString;
        errorString << "OCL Error: parameter \"" << key << "\" used but not set" << std::endl;
        throw SGPP::base::operation_exception(errorString.str());
    }
    return this->parameters[key];
}

bool ConfigurationParameters::getAsBoolean(const std::string &key) {
    if (this->parameters.count(key) != 1) {
        std::stringstream errorString;
        errorString << "OCL Error: parameter \"" << key << "\" used but not set" << std::endl;
        throw SGPP::base::operation_exception(errorString.str());
    }
    bool asBool;
    std::stringstream converter(parameters[key]);
    converter >> std::boolalpha;
    converter >> asBool;
    return asBool;
}

uint64_t ConfigurationParameters::getAsUnsigned(const std::string &key) {
    if (this->parameters.count(key) != 1) {
        std::stringstream errorString;
        errorString << "OCL Error: parameter \"" << key << "\" used but not set" << std::endl;
        throw SGPP::base::operation_exception(errorString.str());
    }
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

void ConfigurationParameters::readFromMap(std::map<std::string, std::string> &parametersMap) {
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

std::vector<std::string> ConfigurationParameters::getAsList(const std::string &key) {

    std::string &listAttribute = this->get(key);

    auto splitted = this->split(listAttribute, ',');

    return splitted;
}


}
}
