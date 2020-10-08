// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/base/tools/json/Token.hpp>

#include <memory>
#include <vector>
#include <string>

namespace json {

class Node {
 protected:
  static const int SERIALIZE_INDENT = 3;

 public:
  // only relevant if the parent is a DictNode, managed by the DictNode
  size_t orderedKeyIndex;

  // managed by the parent, i.e. DictNode or ListNode
  Node* parent;

  Node();

  Node(const Node& right) = default;

  virtual ~Node() = default;

  virtual Node& operator=(const Node& right);

  virtual void parse(std::vector<Token>& stream) = 0;

  virtual void serialize(std::ostream& outFile, size_t indentWidth) = 0;

  virtual Node& operator[](const std::string& key);

  virtual Node& operator[](const size_t index);

  virtual std::string& get();

  virtual void set(const std::string& value);

  virtual double getDouble();

  virtual void setDouble(double doubleValue);

  virtual uint64_t getUInt();

  virtual void setUInt(uint64_t uintValue);

  virtual int64_t getInt();

  virtual void setInt(int64_t intValue);

  virtual bool getBool();

  virtual void setBool(bool boolValue);

  //  virtual JSONNode &getItem(size_t index);

  virtual size_t size() = 0;

  virtual void addValue(std::unique_ptr<Node> node);

  virtual void addAttribute(const std::string& name, std::unique_ptr<Node> node);

  virtual std::unique_ptr<Node> removeValue(size_t index);

  virtual std::unique_ptr<Node> removeAttribute(const std::string name);

  virtual Node* clone() = 0;

  // returns the node to which the attribute was added
  virtual Node& addTextAttr(const std::string& name, const std::string& value);

  // returns the node to which the attribute was added
  virtual Node& addIDAttr(const std::string& name, const std::string& value);

  // returns the node to which the attribute was added
  // cast internally to string, prevents the boolean overload from being used, if the value is a
  // string literal
  virtual Node& addIDAttr(const std::string& name, const char* value);

  // returns the node to which the attribute was added
  virtual Node& addIDAttr(const std::string& name, const double& value);

  // returns the node to which the attribute was added
  virtual Node& addIDAttr(const std::string& name, const uint64_t& value);

  // returns the node to which the attribute was added
  virtual Node& addIDAttr(const std::string& name, const int64_t& value);

  // returns the node to which the attribute was added
  virtual Node& addIDAttr(const std::string& name, const bool& value);

  // returns created dict node
  virtual Node& addDictAttr(const std::string& name);

  // returns created list node
  virtual Node& addListAttr(const std::string& name);

  // returns the node to which the attribute was added
  // replaces a node, adds a new node, if the node does not exist, the old node is deleted
  virtual Node& replaceTextAttr(const std::string& name, const std::string& value);

  // returns the node to which the attribute was added
  // replaces a node, adds a new node, if the node does not exist, the old node is deleted
  virtual Node& replaceIDAttr(const std::string& name, const std::string& value);

  // returns the node to which the attribute was added
  // replaces a node, adds a new node, if the node does not exist, the old node is deleted
  // cast internally to string, prevents the boolean overload from being used, if the value is a
  // string literal
  virtual Node& replaceIDAttr(const std::string& name, const char* value);

  // returns the node to which the attribute was added
  // replaces a node, adds a new node, if the node does not exist, the old node is deleted
  virtual Node& replaceIDAttr(const std::string& name, const double& value);

  // returns the node to which the attribute was added
  // replaces a node, adds a new node, if the node does not exist, the old node is deleted
  virtual Node& replaceIDAttr(const std::string& name, const uint64_t& value);

  // returns the node to which the attribute was added
  // replaces a node, adds a new node, if the node does not exist, the old node is deleted
  virtual Node& replaceIDAttr(const std::string& name, const int64_t& value);

  // returns the node to which the attribute was added
  // replaces a node, adds a new node, if the node does not exist, the old node is deleted
  virtual Node& replaceIDAttr(const std::string& name, const bool& value);

  // returns created dict node
  // replaces a node, adds a new node, if the node does not exist, the old node is deleted
  virtual Node& replaceDictAttr(const std::string& name);

  // returns created list node
  // replaces a node, adds a new node, if the node does not exist, the old node is deleted
  virtual Node& replaceListAttr(const std::string& name);

  // returns created dict node
  virtual Node& addDictValue();

  // returns created dict node
  virtual Node& addListValue();

  // returns the list node to which the value was added
  virtual Node& addTextValue(const std::string& value);

  // returns the list node to which the value was added
  virtual Node& addIdValue(const std::string& value);

  // returns the list node to which the value was added
  // cast internally to string, prevents the boolean overload from being used, if the value is a
  // string literal
  virtual Node& addIdValue(const char* value);

  // returns the list node to which the value was added
  virtual Node& addIdValue(const double& value);

  // returns the list node to which the value was added
  virtual Node& addIdValue(const uint64_t& value);

  // returns the list node to which the value was added
  virtual Node& addIdValue(const int64_t& value);

  // returns the list node to which the value was added
  virtual Node& addIdValue(const bool& value);

  virtual bool contains(const std::string& key);

  virtual std::unique_ptr<Node> erase(Node& node);

  virtual std::unique_ptr<Node> erase();

  virtual std::vector<std::string>& keys();
};

}  // namespace json

std::ostream& operator<<(std::ostream& stream, json::Node& node);
