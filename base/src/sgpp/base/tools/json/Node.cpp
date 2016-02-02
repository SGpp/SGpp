/*
 * JSONNode.cpp
 *
 *  Created on: Nov 7, 2015
 *      Author: pfandedd
 */

#include "Node.hpp"

#include "json_exception.hpp"

namespace json {

Node::Node() :
        orderedKeyIndex(0), parent(nullptr) {

}

Node &Node::operator=(const Node& right) {
    //members of this base class are actually copied in the dervied classes, therefore nothing to do
    return *this;
}

Node &Node::operator[](const std::string &key) {
    throw json_exception("operator[] is only implemented for dict nodes");
}

Node &Node::operator[](const size_t index) {
    throw json_exception("operator[] is only implemented for list nodes");
}

std::string &Node::get() {
    throw json_exception("getValue() is only implemented for string and id nodes");
}

void Node::set(const std::string &value) {
    throw json_exception("setValue() is only implemented for string and id nodes");
}

double Node::getDouble() {
    throw json_exception("getNumericValue() is only implemented for id nodes");
}

void Node::setDouble(double doubleValue) {
    throw json_exception("setNumericValue() is only implemented for id nodes");
}

uint64_t Node::getUInt() {
    throw json_exception("getUInt() is only implemented for id nodes");
}

void Node::setUInt(uint64_t uintValue) {
    throw json_exception("setUInt() is only implemented for id nodes");
}

int64_t Node::getInt() {
    throw json_exception("getInt() is only implemented for id nodes");
}

void Node::setInt(int64_t intValue) {
    throw json_exception("setInt() is only implemented for id nodes");
}

bool Node::getBool() {
    throw json_exception("getBool() is only implemented for id nodes");
}

void Node::setBool(bool boolValue) {
    throw json_exception("setBool() is only implemented for id nodes");
}

void Node::addValue(std::unique_ptr<Node> node) {
    throw json_exception("addItem() is only implemented for list nodes");
}

void Node::addAttribute(const std::string &name, std::unique_ptr<Node> node) {
    throw json_exception("addAttribute() is only implemented for dict nodes");
}

std::unique_ptr<Node> Node::removeValue(size_t index) {
    throw json_exception("removeItem() is only implemented for attribute and list nodes");
}

std::unique_ptr<Node> Node::removeAttribute(const std::string &name) {
    throw json_exception("removeAttribute() is only implemented for dict nodes");
}

// returns the node to which the attribute was added
Node &Node::addTextAttr(const std::string &name, const std::string &value) {
    throw json_exception("addTextAttr() is only implemented for dict nodes");
}

// returns the node to which the attribute was added
Node &Node::addIDAttr(const std::string &name, const std::string &value) {
    throw json_exception("addIDAttr() is only implemented for dict nodes");
}

// returns the node to which the attribute was added
Node &Node::addIDAttr(const std::string &name, const char *value) {
    throw json_exception("addIDAttr() is only implemented for dict nodes");
}

// returns the node to which the attribute was added
Node &Node::addIDAttr(const std::string &name, const double &numericValue) {
    throw json_exception("addIDAttr() is only implemented for dict nodes");
}

// returns the node to which the attribute was added
Node &Node::addIDAttr(const std::string &name, const uint64_t &value) {
    throw json_exception("addIDAttr() is only implemented for dict nodes");
}

// returns the node to which the attribute was added
Node &Node::addIDAttr(const std::string &name, const int64_t &value) {
    throw json_exception("addIDAttr() is only implemented for dict nodes");
}

// returns the node to which the attribute was added
Node &Node::addIDAttr(const std::string &name, const bool &value) {
    throw json_exception("addIDAttr() is only implemented for dict nodes");
}

// returns created dict node
Node &Node::addDictAttr(const std::string &name) {
    throw json_exception("addDictAttr() is only implemented for dict nodes");
}

// returns created list node
Node &Node::addListAttr(const std::string &name) {
    throw json_exception("addListAttr() is only implemented for dict nodes");
}

// returns the node to which the attribute was added
// replaces a node, adds a new node, if the node does not exist, the old node is deleted
Node &Node::replaceTextAttr(const std::string &name, const std::string &value) {
    throw json_exception("replaceTextAttr() is only implemented for dict nodes");
}

// returns the node to which the attribute was added
// replaces a node, adds a new node, if the node does not exist, the old node is deleted
Node &Node::replaceIDAttr(const std::string &name, const std::string &value) {
    throw json_exception("replaceIDAttr() is only implemented for dict nodes");
}

// returns the node to which the attribute was added
// replaces a node, adds a new node, if the node does not exist, the old node is deleted
Node &Node::replaceIDAttr(const std::string &name, const double &value) {
    throw json_exception("replaceIDAttr() is only implemented for dict nodes");
}

// returns the node to which the attribute was added
// replaces a node, adds a new node, if the node does not exist, the old node is deleted
Node &Node::replaceIDAttr(const std::string &name, const uint64_t &value) {
    throw json_exception("replaceIDAttr() is only implemented for dict nodes");
}

// returns the node to which the attribute was added
// replaces a node, adds a new node, if the node does not exist, the old node is deleted
Node &Node::replaceIDAttr(const std::string &name, const int64_t &value) {
    throw json_exception("replaceIDAttr() is only implemented for dict nodes");
}

// returns the node to which the attribute was added
// replaces a node, adds a new node, if the node does not exist, the old node is deleted
Node &Node::replaceIDAttr(const std::string &name, const bool &value) {
    throw json_exception("replaceIDAttr() is only implemented for dict nodes");
}

// returns created dict node
// replaces a node, adds a new node, if the node does not exist, the old node is deleted
Node &Node::replaceDictAttr(const std::string &name) {
    throw json_exception("replaceDictAttr() is only implemented for dict nodes");
}

// returns created list node
// replaces a node, adds a new node, if the node does not exist, the old node is deleted
Node &Node::replaceListAttr(const std::string &name) {
    throw json_exception("replaceListAttr() is only implemented for dict nodes");
}

// returns created dict node
Node &Node::addDictValue() {
    throw json_exception("addDictValue() is only implemented for list nodes");
}

// returns created dict node
Node &Node::addListValue() {
    throw json_exception("addListValue() is only implemented for list nodes");
}

// returns the list node to which the value was added
Node &Node::addTextValue(const std::string &value) {
    throw json_exception("addTextValue() is only implemented for list nodes");
}

// returns the list node to which the value was added
Node &Node::addIdValue(const std::string &value) {
    throw json_exception("addIdValue() is only implemented for list nodes");
}

// returns the list node to which the value was added
Node &Node::addIdValue(const char *value) {
    throw json_exception("addIdValue() is only implemented for list nodes");
}

// returns the list node to which the value was added
Node &Node::addIdValue(const double &value) {
    throw json_exception("addIdValue() is only implemented for list nodes");
}

// returns the list node to which the value was added
Node &Node::addIdValue(const uint64_t &value) {
    throw json_exception("addIdValue() is only implemented for list nodes");
}

// returns the list node to which the value was added
Node &Node::addIdValue(const int64_t &value) {
    throw json_exception("addIdValue() is only implemented for list nodes");
}

// returns the list node to which the value was added
Node &Node::addIdValue(const bool &value) {
    throw json_exception("addIdValue() is only implemented for list nodes");
}

bool Node::contains(const std::string &key) {
    throw json_exception("contains() is only implemented for dict nodes");
}

std::unique_ptr<Node> Node::erase(Node &node) {
    throw json_exception("erase(node) is only implemented for list and dict nodes");
}

std::unique_ptr<Node> Node::erase() {
    if (this->parent == nullptr) {
        throw json_exception("erase(): has no parent");
    }
    std::unique_ptr<Node> self = this->parent->erase(*this);
    this->parent = nullptr;
    return self;
}

std::vector<std::string> &Node::keys() {
    throw json_exception("keys() is only implemented for dict nodes");
}

}

std::ostream& operator<<(std::ostream& stream, json::Node& node) {
    node.serialize(stream, 0);
    return stream;
}
