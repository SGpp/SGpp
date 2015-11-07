/*
 * JSONIDNode.cpp
 *
 *  Created on: Nov 7, 2015
 *      Author: pfandedd
 */

#include "JSONIDNode.hpp"

#include <fstream>

namespace SGPP {
namespace base {

JSONIDNode::JSONIDNode() :
    value() {
}

void JSONIDNode::parse(std::vector<JSONToken> &stream) {
//create new text node
  if (stream[0].type == JSONTokenType::ID) {
    this->value = stream[0].value;
    stream.erase(stream.begin());
  } else {
    throw; //expected colon
  }
}

std::string &JSONIDNode::getValue() {
  return this->value;
}

void JSONIDNode::serialize(std::ofstream &outFile, size_t indentWidth) {
  outFile << this->value;
}

size_t JSONIDNode::size() {
  return 1;
}

}
}
