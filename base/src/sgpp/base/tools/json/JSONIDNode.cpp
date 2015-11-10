/*
 * JSONIDNode.cpp
 *
 *  Created on: Nov 7, 2015
 *      Author: pfandedd
 */

#include <fstream>
#include <string>
#include <sstream>

#include "JSONIDNode.hpp"
#include "json_exception.hpp"

namespace json {

JSONIDNode::JSONIDNode() :
    value(), isNumber(false), numericValue(0) {
}

void JSONIDNode::parse(std::vector<JSONToken> &stream) {
//create new text node
  if (stream[0].type == JSONTokenType::ID) {
    this->value = stream[0].value;
    stream.erase(stream.begin());

    this->tryInterpretAsNumber();
  } else {
    throw json_exception(stream[0], "expected id");
  }
}

void JSONIDNode::tryInterpretAsNumber() {
  try {
    std::string::size_type size;
    double asNumber = stod(this->value, &size);
    if (this->value.size() == size) {
      this->numericValue = asNumber;
      this->isNumber = true;
    } else {
      this->isNumber = false;
    }
  } catch (std::invalid_argument &e) {
    isNumber = false;
  }
}

std::string &JSONIDNode::get() {
  return this->value;
}

void JSONIDNode::set(const std::string &value) {
  this->value = value;

  this->tryInterpretAsNumber();
}

double JSONIDNode::getNumeric() {
  if (this->isNumber) {
    return this->numericValue;
  } else {
    throw json_exception("node has not a numerical value");
  }
}

void JSONIDNode::setNumeric(double numericValue) {
  this->numericValue = numericValue;
  this->isNumber = true;
  std::stringstream stringstream;
  stringstream << numericValue;
  this->value = stringstream.str();
}

void JSONIDNode::serialize(std::ofstream &outFile, size_t indentWidth) {
  outFile << this->value;
}

size_t JSONIDNode::size() {
  return 1;
}

JSONNode *JSONIDNode::clone() {
  JSONIDNode *newNode = new JSONIDNode(*this);
  return newNode;
}

}
