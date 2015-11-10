/*
 * JSONStringNode.cpp
 *
 *  Created on: Nov 7, 2015
 *      Author: pfandedd
 */

#include "JSONTextNode.hpp"

#include <fstream>

#include "json_exception.hpp"

namespace json {

JSONTextNode::JSONTextNode() :
    value() {
}

void JSONTextNode::parse(std::vector<JSONToken> &stream) {
//create new text node
  if (stream[0].type == JSONTokenType::STRING) {
    this->value = stream[0].value;
    stream.erase(stream.begin());
  } else {
    throw json_exception("expected string value");
  }
}

void JSONTextNode::serialize(std::ofstream &outFile, size_t indentWidth) {
  outFile << "\"" << this->value << "\"";
}

std::string &JSONTextNode::get() {
  return this->value;
}

void JSONTextNode::set(const std::string &value) {
  this->value = value;
}

size_t JSONTextNode::size() {
  return 1;
}

JSONNode *JSONTextNode::clone() {
  JSONTextNode *newNode = new JSONTextNode(*this);
  return newNode;
}

}
