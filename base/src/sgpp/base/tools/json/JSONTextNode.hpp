/*
 * JSONStringNode.hpp
 *
 *  Created on: Nov 7, 2015
 *      Author: pfandedd
 */

#pragma once

#include "JSONNode.hpp"

namespace json {

class JSONTextNode: public JSONNode {
private:

  std::string value;

public:
  JSONTextNode();

  virtual void parse(std::vector<JSONToken> &stream) override;

  virtual void serialize(std::ofstream &outFile, size_t indentWidth) override;

  virtual std::string &get() override;

  virtual void set(const std::string &value) override;

  virtual size_t size() override;

  virtual JSONNode *clone() override;
};

}
