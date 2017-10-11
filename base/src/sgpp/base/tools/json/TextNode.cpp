// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#include <sgpp/base/tools/json/TextNode.hpp>
#include <sgpp/base/tools/json/json_exception.hpp>

#include <fstream>
#include <sstream>
#include <string>
#include <vector>

namespace json {

TextNode::TextNode()
    : value(),
      isDouble(false),
      doubleValue(0.0),
      isUnsigned(false),
      unsignedValue(0),
      isSigned(false),
      signedValue(0),
      isBool(false),
      boolValue(false) {}

Node& TextNode::operator=(const Node& right) {
  const TextNode& textNode = dynamic_cast<const TextNode&>(right);
  this->operator=(textNode);
  return *this;
}

void TextNode::setupInternalType() {
  // try validating as bool
  if (this->value.compare("true") == 0) {
    this->isBool = true;
    this->boolValue = true;
  } else if (this->value.compare("false") == 0) {
    this->isBool = true;
    this->boolValue = false;
  }
  {
    std::stringstream conv_stream;
    conv_stream << this->value;
    uint64_t r;
    conv_stream >> r;
    if (!conv_stream.fail()) {
      this->isUnsigned = true;
      this->unsignedValue = r;
    }
  }
  {
    std::stringstream conv_stream;
    conv_stream << this->value;
    int64_t r;
    conv_stream >> r;
    if (!conv_stream.fail()) {
      this->isSigned = true;
      this->signedValue = r;
    }
  }
  {
    std::stringstream conv_stream;
    conv_stream << this->value;
    double r;
    conv_stream >> r;
    if (!conv_stream.fail()) {
      this->isDouble = true;
      this->doubleValue = r;
    }
  }
}

void TextNode::parse(std::vector<Token>& stream) {
  // create new text node
  if (stream[0].type == TokenType::STRING) {
    this->value = stream[0].value;
    stream.erase(stream.begin());

    this->setupInternalType();
  } else {
    throw json_exception("expected string value");
  }
}

void TextNode::serialize(std::ostream& outFile, size_t indentWidth) {
  outFile << "\"" << this->value << "\"";
}

std::string& TextNode::get() { return this->value; }

void TextNode::set(const std::string& value) {
  this->value = value;
  this->setupInternalType();
}

double TextNode::getDouble() {
  //    if (this->internalType == InternalIDType::DOUBLE) {
  if (this->isDouble) {
    return this->doubleValue;
  } else {
    throw json_exception("node has not a numerical value");
  }
}

void TextNode::setDouble(double numericValue) {
  //    this->doubleValue = numericValue;
  //    this->internalType = InternalIDType::DOUBLE;
  std::stringstream stringstream;
  stringstream << numericValue;
  this->value = stringstream.str();
  this->setupInternalType();
}

uint64_t TextNode::getUInt() {
  //    if (this->internalType == InternalIDType::UINT) {
  if (this->isUnsigned) {
    return this->unsignedValue;
  } else {
    throw json_exception("node has not an unsigned integer value");
  }
}

void TextNode::setUInt(uint64_t uintValue) {
  //    this->unsignedValue = uintValue;
  //    this->internalType = InternalIDType::UINT;
  std::stringstream stringstream;
  stringstream << uintValue;
  this->value = stringstream.str();
  this->setupInternalType();
}

int64_t TextNode::getInt() {
  //    if (this->internalType == InternalIDType::INT) {
  if (this->isSigned) {
    return this->signedValue;
  } else {
    throw json_exception("node has not an integer value");
  }
}

void TextNode::setInt(int64_t intValue) {
  //    this->unsignedValue = intValue;
  //    this->internalType = InternalIDType::INT;
  std::stringstream stringstream;
  stringstream << intValue;
  this->value = stringstream.str();
  this->setupInternalType();
}

bool TextNode::getBool() {
  //    if (this->internalType == InternalIDType::BOOL) {
  if (this->isBool) {
    return this->boolValue;
  } else {
    throw json_exception("node has not a bool value");
  }
}

void TextNode::setBool(bool boolValue) {
  //    this->boolValue = boolValue;
  //    this->internalType = InternalIDType::BOOL;
  if (boolValue) {
    this->value = std::string("true");
  } else {
    this->value = std::string("false");
  }

  this->setupInternalType();
}

size_t TextNode::size() { return 1; }

Node* TextNode::clone() {
  TextNode* newNode = new TextNode(*this);
  return newNode;
}

}  // namespace json
