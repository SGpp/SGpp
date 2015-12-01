/*
 * JSONStringNode.hpp
 *
 *  Created on: Nov 7, 2015
 *      Author: pfandedd
 */

#pragma once

#include "Node.hpp"

namespace json {

class TextNode: public Node {
private:

  std::string value;

public:
  TextNode();

  virtual void parse(std::vector<Token> &stream) override;

  virtual void serialize(std::ofstream &outFile, size_t indentWidth) override;

  virtual std::string &get() override;

  virtual void set(const std::string &value) override;

  virtual size_t size() override;

  virtual Node *clone() override;
};

}
