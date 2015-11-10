/*
 * JSONNode.cpp
 *
 *  Created on: Nov 7, 2015
 *      Author: pfandedd
 */

#include "JSONNode.hpp"
#include "json_exception.hpp"

namespace json {

JSONNode &JSONNode::operator[](std::string key) {
  throw json_exception("operator[] is only implemented for dict nodes");
}

JSONNode &JSONNode::operator[](size_t index) {
  throw json_exception("operator[] is only implemented for list nodes");
}

std::string &JSONNode::get() {
  throw json_exception("getValue() is only implemented for string and id nodes");
}

void JSONNode::set(const std::string &value) {
  throw json_exception("setValue() is only implemented for string and id nodes");
}

double JSONNode::getNumeric() {
  throw json_exception("getNumericValue() is only implemented for id nodes");
}

void JSONNode::setNumeric(double numericValue) {
  throw json_exception("setNumericValue() is only implemented for id nodes");
}

void JSONNode::addValue(std::unique_ptr<JSONNode> node) {
  throw json_exception("addItem() is only implemented for list nodes");
}

void JSONNode::addAttribute(const std::string &name, std::unique_ptr<JSONNode> node) {
  throw json_exception("addAttribute() is only implemented for dict nodes");
}

std::unique_ptr<JSONNode> JSONNode::removeValue(size_t index) {
  throw json_exception("removeItem() is only implemented for attribute and list nodes");
}

std::unique_ptr<JSONNode> JSONNode::removeAttribute(const std::string &name) {
  throw json_exception("removeAttribute() is only implemented for dict nodes");
}

// returns the node to which the attribute was added
JSONNode &JSONNode::addTextAttr(const std::string &name, const std::string &value) {
  throw json_exception("addTextAttr() is only implemented for dict nodes");
}

// returns the node to which the attribute was added
JSONNode &JSONNode::addIDAttr(const std::string &name, const std::string &value) {
  throw json_exception("addIDAttr() is only implemented for dict nodes");
}

// returns the node to which the attribute was added
JSONNode &JSONNode::addIDAttr(const std::string &name, const double &numericValue) {
  throw json_exception("addIDAttr() is only implemented for dict nodes");
}

// returns created dict node
JSONNode &JSONNode::addDictAttr(const std::string &name) {
  throw json_exception("addDictAttr() is only implemented for dict nodes");
}

// returns created list node
JSONNode &JSONNode::addListAttr(const std::string &name) {
  throw json_exception("addListAttr() is only implemented for dict nodes");
}

// returns created dict node
JSONNode &JSONNode::addDictValue() {
  throw json_exception("addDictValue() is only implemented for list nodes");
}

// returns created dict node
JSONNode &JSONNode::addListValue() {
  throw json_exception("addListValue() is only implemented for list nodes");
}

// returns the list node to which the value was added
JSONNode &JSONNode::addTextValue(const std::string &value) {
  throw json_exception("addTextValue() is only implemented for list nodes");
}

// returns the list node to which the value was added
JSONNode &JSONNode::addIdValue(const std::string &value) {
  throw json_exception("addIdValue() is only implemented for list nodes");
}

// returns the list node to which the value was added
JSONNode &JSONNode::addIdValue(const double &value) {
  throw json_exception("addIdValue() is only implemented for list nodes");
}

std::unique_ptr<JSONNode> JSONNode::erase(JSONNode &node) {
  throw json_exception("erase(node) is only implemented for list and dict nodes");
}

std::unique_ptr<JSONNode> JSONNode::erase() {
  if (this->parent == nullptr) {
    throw json_exception("erase(): has no parent");
  }
  std::unique_ptr<JSONNode> self = this->parent->erase(*this);
  this->parent = nullptr;
  return self;
}

}
