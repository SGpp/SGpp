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

#include <sgpp/base/tools/json/JSON.hpp>

namespace SGPP {
namespace base {

class OperationConfiguration: public json::JSON {
public:

    OperationConfiguration();

    OperationConfiguration(const std::string &fileName);

    virtual OperationConfiguration *clone();

};

}
}

