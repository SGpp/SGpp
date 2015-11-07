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

#include "JSONToken.hpp"
#include "JSONAttributeNode.hpp"

namespace SGPP {
namespace base {

class JSON: public JSONAttributeNode {
private:
//  JSONAttributeNode root;
  std::vector<JSONToken> tokenize(std::string &input);

public:
  JSON(std::string fileName);

  void serialize(std::string outFileName);

  using JSONAttributeNode::serialize;

};

}
}
