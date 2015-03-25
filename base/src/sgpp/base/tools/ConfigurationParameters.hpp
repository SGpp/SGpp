/*
 * ConfigurationParser.hpp
 *
 *  Created on: Mar 25, 2015
 *      Author: pfandedd
 */

#pragma once

#include <vector>
#include <map>

#include <sgpp/globaldef.hpp>

namespace SGPP {
namespace base {

	class ConfigurationParameters {
	private:
		std::map<std::string, std::string> parameters;
	public:
		ConfigurationParameters(std::string fileName, std::map<std::string, std::string> defaultParameters = std::map<std::string, std::string>());

		std::string operator[](std::string key);

		bool getAsBoolean(std::string key);
		uint64_t getAsUnsigned(std::string key);

	private:
		std::vector<std::string> split(const std::string &s, char delim);
	};

}
}


