/*
 * ConfigurationParser.hpp
 *
 *  Created on: Mar 25, 2015
 *      Author: pfandedd
 */

#pragma once

#include <vector>
#include <map>
#include <string>
#include <memory>
#include <iostream>

#include <sgpp/globaldef.hpp>

namespace SGPP {
namespace base {

class ConfigurationParameters {
protected:
    std::map<std::string, std::string> parameters;
    bool useDefaults;
public:
    ConfigurationParameters();

    ConfigurationParameters(std::string fileName,
            std::map<std::string, std::string> defaultParameters = std::map<std::string, std::string>());

    virtual ~ConfigurationParameters();

    bool doUseDefaults();

    void set(const std::string key, std::string value);

    std::string &get(const std::string &key);

    bool getAsBoolean(const std::string &key);

    uint64_t getAsUnsigned(const std::string &key);

    void readFromMap(std::map<std::string, std::string> &parametersMap);

    void readFromFile(std::string fileName);

    virtual std::shared_ptr<ConfigurationParameters> clone() = 0;

    void debug() {
        std::cout << "count: " << this->parameters.size() << std::endl;
    }

private:
    std::vector<std::string> split(const std::string& s, char delim);
};

}
}

