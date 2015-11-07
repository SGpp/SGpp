/*
 * JSONStringNode.cpp
 *
 *  Created on: Nov 7, 2015
 *      Author: pfandedd
 */

#include "JSONStringNode.hpp"

#include <fstream>

namespace SGPP {
namespace base {

JSONStringNode::JSONStringNode() :
    value() {
}

void JSONStringNode::parse(std::vector<JSONToken> &stream) {
//create new text node
  if (stream[0].type == JSONTokenType::STRING) {
    this->value = stream[0].value;
    stream.erase(stream.begin());
  } else {
    throw; //expected colon
  }
}

void JSONStringNode::serialize(std::ofstream &outFile, size_t indentWidth) {
  outFile << "\"" << this->value << "\"";
}

std::string &JSONStringNode::getValue() {
  return this->value;
}

size_t JSONStringNode::size() {
  return 1;
}

}
}
