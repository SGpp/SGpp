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

#include <sgpp/base/tools/ConfigurationParameters.hpp>

namespace SGPP {
namespace base {

class OCLConfigurationParameters: public ConfigurationParameters {
 public:
  OCLConfigurationParameters(std::string fileName,
                             std::map<std::string, std::string> defaultParameters);

  OCLConfigurationParameters();

  virtual ~OCLConfigurationParameters();

  virtual std::shared_ptr<ConfigurationParameters> clone() override;
};

}  // namespace base
}  // namespace SGPP

