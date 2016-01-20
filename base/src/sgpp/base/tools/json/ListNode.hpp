/*
 * JSONListNode.hpp
 *
 *  Created on: Nov 7, 2015
 *      Author: pfandedd
 */

#pragma once

#include <memory>
#include <vector>

#include <sgpp/base/tools/json/Node.hpp>

namespace json {

class ListNode: public Node {
private:

  std::vector<std::unique_ptr<Node>> list;

public:
  ListNode();

  ListNode(const ListNode &original);

  ListNode &operator=(const ListNode &right);

  void parse(std::vector<Token> &stream) override;

  virtual void serialize(std::ofstream &outFile, size_t indentWidth) override;

  virtual Node &operator[](size_t index) override;

  virtual size_t size() override;

  virtual void addValue(std::unique_ptr<Node> node) override;

  virtual std::unique_ptr<Node> removeValue(size_t index) override;

  virtual Node *clone() override;

  // returns created dict node
  virtual Node &addDictValue() override;

  // returns created dict node
  virtual Node &addListValue() override;

  // returns the list node to which the value was added
  virtual Node &addTextValue(const std::string &value) override;

  // returns the list node to which the value was added
  virtual Node &addIdValue(const std::string &value) override;

  // returns the list node to which the value was added
  virtual Node &addIdValue(const char *value) override;

  // returns the list node to which the value was added
  virtual Node &addIdValue(const double &value) override;

  // returns the list node to which the value was added
  virtual Node &addIdValue(const uint64_t &value) override;

  // returns the list node to which the value was added
  virtual Node &addIdValue(const int64_t &value) override;

  // returns the list node to which the value was added
  virtual Node &addIdValue(const bool &value) override;

  virtual std::unique_ptr<Node> erase(Node &node) override;

};

}
