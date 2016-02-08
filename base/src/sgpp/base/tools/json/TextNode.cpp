// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/tools/json/TextNode.hpp>
#include <sgpp/base/tools/json/json_exception.hpp>

#include <fstream>
#include <string>
#include <vector>

namespace json {

TextNode::TextNode() :
  value() {
}

void TextNode::parse(std::vector<Token>& stream) {
  // create new text node
  if (stream[0].type == TokenType::STRING) {
    this->value = stream[0].value;
    stream.erase(stream.begin());
  } else {
    throw json_exception("expected string value");
  }
}

void TextNode::serialize(std::ofstream& outFile, size_t indentWidth) {
  outFile << "\"" << this->value << "\"";
}

std::string& TextNode::get() {
  return this->value;
}

void TextNode::set(const std::string& value) {
  this->value = value;
}

size_t TextNode::size() {
  return 1;
}

Node* TextNode::clone() {
  TextNode* newNode = new TextNode(*this);
  return newNode;
}

}  // namespace json
