// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/base/tools/json/Node.hpp>

#include <map>
#include <memory>
#include <string>
#include <vector>

namespace json {

class DictNode : public Node {
 protected:
  std::map<std::string, std::unique_ptr<Node>> attributes;

  std::vector<std::string> keyOrder;

 public:
  DictNode();

  DictNode(const DictNode& original);

  DictNode& operator=(const DictNode& right);

  Node& operator=(const Node& right) override;

  void parse(std::vector<Token>& stream) override;

  void parseAttributes(std::vector<Token>& stream);

  void serialize(std::ostream& outFile, size_t indentWidth) override;

  Node& operator[](const std::string& key)override;

  size_t size() override;

  Node* clone() override;

  void addAttribute(const std::string& name, std::unique_ptr<Node> node) override;

  std::unique_ptr<Node> removeAttribute(const std::string& name) override;

  // returns the node to which the attribute was added
  Node& addTextAttr(const std::string& name, const std::string& value) override;

  // returns the node to which the attribute was added
  Node& addIDAttr(const std::string& name, const std::string& value) override;

  // returns the node to which the attribute was added
  // cast internally to string, prevents the boolean overload from being used, if the value is a
  // string literal
  Node& addIDAttr(const std::string& name, const char* value) override;

  // returns the node to which the attribute was added
  Node& addIDAttr(const std::string& name, const double& value) override;

  // returns the node to which the attribute was added
  Node& addIDAttr(const std::string& name, const uint64_t& value) override;

  // returns the node to which the attribute was added
  Node& addIDAttr(const std::string& name, const int64_t& value) override;

  // returns the node to which the attribute was added
  Node& addIDAttr(const std::string& name, const bool& value) override;

  // returns created dict node
  Node& addDictAttr(const std::string& name) override;

  // returns created list node
  Node& addListAttr(const std::string& name) override;

  // returns the node to which the attribute was added
  // replaces a node, adds a new node, if the node does not exist,
  // the old node is deleted
  Node& replaceTextAttr(const std::string& name, const std::string& value) override;

  // returns the node to which the attribute was added
  // replaces a node, adds a new node, if the node does not exist,
  // the old node is deleted
  Node& replaceIDAttr(const std::string& name, const std::string& value) override;

  // returns the node to which the attribute was added
  // replaces a node, adds a new node, if the node does not exist,
  // the old node is deleted
  // cast internally to string, prevents the boolean overload from being used, if the value is a
  // string literal
  Node& replaceIDAttr(const std::string& name, const char* value) override;

  // returns the node to which the attribute was added
  // replaces a node, adds a new node, if the node does not exist,
  // the old node is deleted
  Node& replaceIDAttr(const std::string& name, const double& value) override;

  // returns the node to which the attribute was added
  // replaces a node, adds a new node, if the node does not exist,
  // the old node is deleted
  Node& replaceIDAttr(const std::string& name, const uint64_t& value) override;

  // returns the node to which the attribute was added
  // replaces a node, adds a new node, if the node does not exist,
  // the old node is deleted
  Node& replaceIDAttr(const std::string& name, const int64_t& value) override;

  // returns the node to which the attribute was added
  // replaces a node, adds a new node, if the node does not exist,
  // the old node is deleted
  Node& replaceIDAttr(const std::string& name, const bool& value) override;

  // returns created dict node
  // replaces a node, adds a new node, if the node does not exist,
  // the old node is deleted
  Node& replaceDictAttr(const std::string& name) override;

  // returns created list node
  // replaces a node, adds a new node, if the node does not exist,
  // the old node is deleted
  Node& replaceListAttr(const std::string& name) override;

  bool contains(const std::string& key) override;

  std::unique_ptr<Node> erase(Node& node) override;

  std::vector<std::string>& keys() override;
};

}  // namespace json
