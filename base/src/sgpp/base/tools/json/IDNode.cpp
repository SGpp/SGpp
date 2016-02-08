/*
 * JSONIDNode.cpp
 *
 *  Created on: Nov 7, 2015
 *      Author: pfandedd
 */

#include "IDNode.hpp"

#include <fstream>
#include <string>
#include <sstream>

#include "json_exception.hpp"

namespace json {

IDNode::IDNode() :
  value(),
  //        internalType(InternalIDType::ID),
  isDouble(false), doubleValue(0.0), isUnsigned(false), unsignedValue(0),
  isSigned(false), signedValue(0), isBool(
    false), boolValue(false) {
}

Node& IDNode::operator=(const Node& right) {
  const IDNode& idNode = dynamic_cast<const IDNode&>(right);
  this->operator =(idNode);
  return *this;
}

void IDNode::parse(std::vector<Token>& stream) {
  //create new text node
  if (stream[0].type == TokenType::ID) {
    this->value = stream[0].value;
    stream.erase(stream.begin());

    this->setupInternalType();
  } else {
    throw json_exception(stream[0], "expected id");
  }
}

void IDNode::setupInternalType() {
  //    this->internalType = InternalIDType::ID;

  //try validating as bool
  if (this->value.compare("true") == 0) {
    //        this->internalType = InternalIDType::BOOL;
    this->isBool = true;
    this->boolValue = true;
    return;
  } else if (this->value.compare("false") == 0) {
    //        this->internalType = InternalIDType::BOOL;
    this->isBool = true;
    this->boolValue = false;
    return;
  }

  //try validating as unsigned integer
  try {
    std::string::size_type size;
    uint64_t asUnsigned = stoull(this->value, &size);

    if (this->value.size() == size) {
      this->isUnsigned = true;
      this->unsignedValue = asUnsigned;
      //            this->internalType = InternalIDType::UINT;
      return;
    }
  } catch (std::invalid_argument& e) {
  }

  //try validating as signed integer
  try {
    std::string::size_type size;
    int64_t asSigned = stoll(this->value, &size);

    if (this->value.size() == size) {
      this->isSigned = true;
      this->signedValue = asSigned;
      //            this->internalType = InternalIDType::INT;
      return;
    }
  } catch (std::invalid_argument& e) {
  }

  //try validating as double
  try {
    std::string::size_type size;
    double asDouble = stod(this->value, &size);

    if (this->value.size() == size) {
      this->isDouble = true;
      this->doubleValue = asDouble;
      //            this->internalType = InternalIDType::DOUBLE;
      return;
    }
  } catch (std::invalid_argument& e) {
  }
}

std::string& IDNode::get() {
  return this->value;
}

void IDNode::set(const std::string& value) {
  this->value = value;

  this->setupInternalType();
}

double IDNode::getDouble() {
  //    if (this->internalType == InternalIDType::DOUBLE) {
  if (this->isDouble) {
    return this->doubleValue;
  } else {
    throw json_exception("node has not a numerical value");
  }
}

void IDNode::setDouble(double numericValue) {
  //    this->doubleValue = numericValue;
  //    this->internalType = InternalIDType::DOUBLE;
  std::stringstream stringstream;
  stringstream << numericValue;
  this->value = stringstream.str();
  this->setupInternalType();
}

uint64_t IDNode::getUInt() {
  //    if (this->internalType == InternalIDType::UINT) {
  if (this->isUnsigned) {
    return this->unsignedValue;
  } else {
    throw json_exception("node has not an unsigned integer value");
  }
}

void IDNode::setUInt(uint64_t uintValue) {
  //    this->unsignedValue = uintValue;
  //    this->internalType = InternalIDType::UINT;
  std::stringstream stringstream;
  stringstream << uintValue;
  this->value = stringstream.str();
  this->setupInternalType();
}

int64_t IDNode::getInt() {
  //    if (this->internalType == InternalIDType::INT) {
  if (this->isSigned) {
    return this->signedValue;
  } else {
    throw json_exception("node has not an integer value");
  }
}

void IDNode::setInt(int64_t intValue) {
  //    this->unsignedValue = intValue;
  //    this->internalType = InternalIDType::INT;
  std::stringstream stringstream;
  stringstream << intValue;
  this->value = stringstream.str();
  this->setupInternalType();
}

bool IDNode::getBool() {
  //    if (this->internalType == InternalIDType::BOOL) {
  if (this->isBool) {
    return this->boolValue;
  } else {
    throw json_exception("node has not a bool value");
  }
}

void IDNode::setBool(bool boolValue) {
  //    this->boolValue = boolValue;
  //    this->internalType = InternalIDType::BOOL;
  if (boolValue) {
    this->value = std::string("true");
  } else {
    this->value = std::string("false");
  }

  this->setupInternalType();
}

void IDNode::serialize(std::ostream& outFile, size_t indentWidth) {
  outFile << this->value;
}

size_t IDNode::size() {
  return 1;
}

Node* IDNode::clone() {
  IDNode* newNode = new IDNode(*this);
  return newNode;
}

}
