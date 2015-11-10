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

#include "JSONDictNode.hpp"
#include "JSONToken.hpp"

namespace json {

class JSON: public JSONDictNode {
private:

  std::string fileName;

  std::vector<JSONToken> tokenize(std::string &input);

public:

  JSON(std::string fileName);

  JSON();

  void serialize(std::string outFileName);

  using JSONDictNode::serialize;

};

}
