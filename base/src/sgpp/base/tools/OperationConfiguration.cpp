/*
 * ConfigurationParser.hpp
 *
 *  Created on: Mar 25, 2015
 *      Author: pfandedd
 */

#include <sgpp/base/tools/OperationConfiguration.hpp>

namespace SGPP {
namespace base {

OperationConfiguration::OperationConfiguration(): json::JSON() {

}

OperationConfiguration::OperationConfiguration(const std::string &fileName): json::JSON(fileName) {

}

OperationConfiguration *OperationConfiguration::clone() {
    OperationConfiguration *clone = new OperationConfiguration(*this);
    return clone;
}

}
}

