// Copyright (C) 2008-today The SG++ project
// This file is part of the SG++ project. For conditions of distribution and
// use, please see the copyright notice provided with SG++ or at
// sgpp.sparsegrids.org

#pragma once

#include <sgpp/base/tools/json/Node.hpp>

#include <memory>
#include <vector>
#include <string>

namespace json {

class ListNode: public Node {
 private:
  std::vector<std::unique_ptr<Node>> list;

 public:
  ListNode();

  ListNode(const ListNode& original);

  ListNode& operator=(const ListNode& right);

  Node& operator=(const Node& right) override;

  void parse(std::vector<Token>& stream) override;

  void serialize(std::ostream& outFile, size_t indentWidth) override;

  Node& operator[](const size_t index) override;

  size_t size() override;

  void addValue(std::unique_ptr<Node> node) override;

  std::unique_ptr<Node> removeValue(size_t index) override;

  Node* clone() override;

  // returns created dict node
  Node& addDictValue() override;

  // returns created dict node
  Node& addListValue() override;

  // returns the list node to which the value was added
  Node& addTextValue(const std::string& value) override;

  // returns the list node to which the value was added
  Node& addIdValue(const std::string& value) override;

  // returns the list node to which the value was added
  Node& addIdValue(const char* value) override;

  // returns the list node to which the value was added
  Node& addIdValue(const double& value) override;

  // returns the list node to which the value was added
  Node& addIdValue(const uint64_t& value) override;

  // returns the list node to which the value was added
  Node& addIdValue(const int64_t& value) override;

  // returns the list node to which the value was added
  Node& addIdValue(const bool& value) override;

  std::unique_ptr<Node> erase(Node& node) override;
};

}  // namespace json
