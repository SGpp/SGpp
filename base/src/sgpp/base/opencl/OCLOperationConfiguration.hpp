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

#include <sgpp/base/tools/OperationConfiguration.hpp>

namespace SGPP {
namespace base {

class OCLOperationConfiguration: public OperationConfiguration {
 public:

  OCLOperationConfiguration();

  OCLOperationConfiguration(const std::string& fileName);

  virtual OperationConfiguration* clone() override;

};

}  // namespace base
}  // namespace SGPP

