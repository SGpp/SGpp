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

#include "DictNode.hpp"
#include "Token.hpp"

namespace json {

class JSON: public DictNode {
 private:

  std::string fileName;

  std::vector<Token> tokenize(std::string& input);

 public:

  JSON(const std::string& fileName);

  JSON();

  JSON(const JSON& original);

  virtual JSON* clone();

  void clear();

  void serialize(const std::string& outFileName);

  using DictNode::serialize;

};

}
